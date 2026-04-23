mod kmer;
mod pattern;
mod wordlist;
mod normalization;
mod validator;
mod metrics;
mod output;

use std::path::{Path, PathBuf};
use clap::Parser;
use anyhow::{Context, Result};
use rayon::prelude::*;

use validator::{Validator, TaskSpec};
use metrics::compute_value;
use output::write_outputs;

// ---------------------------------------------------------------------------
// CLI
// ---------------------------------------------------------------------------

#[derive(Parser, Debug)]
#[command(
    name = "oligowords",
    version = env!("CARGO_PKG_VERSION"),
    about = "Sliding-window oligonucleotide composition analysis of DNA sequences.\nReads FASTA input and writes tab-separated output files."
)]
struct Cli {
    /// Input FASTA file (may contain multiple records)
    #[arg(short = 'i', long = "input", value_name = "FASTA", conflicts_with = "path")]
    input: Option<PathBuf>,

    /// Directory to scan for *.fasta / *.fa / *.fna files
    #[arg(long = "path", value_name = "DIR")]
    path: Option<PathBuf>,

    /// Base name or path for output files (appends _oligowords.tsv / .bedgraph). If omitted, uses input filename stem.
    #[arg(long = "output", value_name = "PATH")]
    output: Option<PathBuf>,

    /// Semicolon-separated list of tasks (e.g. "n0_4mer:D;n1_4mer:V")
    #[arg(long, default_value = "n0_4mer:D;n0_4mer:PS;n1_4mer:V;n1_4mer:GV", value_name = "TASKS")]
    task: String,

    /// Sliding-window size in bp
    #[arg(long, value_name = "BP")]
    frame: Option<usize>,

    /// Step size in bp (default: frame / 2)
    #[arg(long, value_name = "BP")]
    step: Option<usize>,

    /// Suppress per-window progress output
    #[arg(short = 'q', long)]
    quiet: bool,

    /// Number of rayon threads (default: 8)
    #[arg(short = 't', long, value_name = "N", default_value_t = 8)]
    threads: usize,

    /// Allow circular genome wrapping (windows extending beyond sequence end wrap to start)
    #[arg(long, default_value_t = false)]
    circular: bool,
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

fn main() -> Result<()> {
    let cli = Cli::parse();

    // Require either -i/--input or --path
    if cli.input.is_none() && cli.path.is_none() {
        eprintln!("Error: either -i/--input or --path must be specified.");
        eprintln!("\nFor help, run: oligowords -h");
        std::process::exit(1);
    }

    // Configure rayon thread pool FIRST — before any par_iter() call, which would
    // otherwise lazily initialise the global pool with the default (all CPUs).
    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.threads)
        .build_global()
        .context("Failed to build rayon thread pool")?;

    // Resolve the list of FASTA files to process
    let file_paths: Vec<PathBuf> = if let Some(ref input) = cli.input {
        vec![input.clone()]
    } else if let Some(ref path) = cli.path {
        fasta_files_in_dir(path)
    } else {
        unreachable!("Checked above that at least one is provided")
    };

    if file_paths.is_empty() {
        eprintln!("No FASTA files found. Use -i/--input to specify a file or --path for a directory.");
        return Ok(());
    }

    // Validate tasks / frame / step
    let mut validator = Validator::new();
    let (tasklist, frame, step) = validator.validate(&cli.task, cli.frame, cli.step)
        .context("Argument validation failed")?;

    // Process files in parallel when there are multiple; each writes its own output file.
    file_paths.par_iter().try_for_each(|fpath| {
        process_file(fpath, &tasklist, frame, step, cli.quiet, cli.circular, cli.output.as_deref())
    })?;

    Ok(())
}

// ---------------------------------------------------------------------------
// File-level processing
// ---------------------------------------------------------------------------

fn process_file(
    fpath: &Path,
    tasklist: &[TaskSpec],
    frame: usize,
    step: usize,
    quiet: bool,
    circular: bool,
    output_base: Option<&Path>,
) -> Result<()> {
    use needletail::parse_fastx_file;

    let mut reader = parse_fastx_file(fpath)
        .with_context(|| format!("Cannot open {}", fpath.display()))?;

    let stem = if let Some(output_path) = output_base {
        output_path.to_string_lossy().to_string()
    } else {
        fpath.file_stem().and_then(|s| s.to_str()).unwrap_or("out").to_string()
    };
    let dir = Path::new(".");
    let tsv_path = dir.join(format!("{}_oligowords.tsv", stem));
    let bg_path = dir.join(format!("{}_oligowords.bedgraph", stem));

    let mut tsv_buf = String::new();
    let mut bg_buf = String::new();

    while let Some(record) = reader.next() {
        let rec = record.with_context(|| format!("Error reading {}", fpath.display()))?;
        let raw_id = std::str::from_utf8(rec.id()).unwrap_or("unknown");
        let name = raw_id.split_whitespace().next().unwrap_or(raw_id);
        // Normalise to uppercase ASCII, stripping non-ACGT characters
        let seq: Vec<u8> = rec.seq()
            .iter()
            .map(|&b| b.to_ascii_uppercase())
            .filter(|&b| matches!(b, b'A' | b'C' | b'G' | b'T'))
            .collect();

        if seq.is_empty() {
            continue;
        }

        if !quiet {
            println!("{}", name);
        }

        process_sequence(name, &seq, tasklist, frame, step, quiet, &mut tsv_buf, &mut bg_buf, circular);
    }

    if !tsv_buf.is_empty() {
        write_outputs(&tsv_path, &bg_path, &tsv_buf, &bg_buf)?;
    }

    Ok(())
}

fn process_sequence(
    name: &str,
    seq: &[u8],
    tasklist: &[TaskSpec],
    frame: usize,
    step: usize,
    quiet: bool,
    tsv_buf: &mut String,
    bg_buf: &mut String,
    circular: bool,
) {
    use std::fmt::Write;
    use pattern::Pattern;

    let seq_len = seq.len();
    let now = chrono_now();

    // TSV header
    tsv_buf.push_str(&format!(
        "# {}\n# {}\n# Sequence length: {} bp.\n# Window: {}\tStep: {}\n#\n# Chrom\tStart\tStop\t",
        now, name, seq_len, frame, step
    ));

    let mut header = String::new();
    for task in tasklist {
        let str_subtr = match &task.subtr {
            Some(s) => format!("-{}_{} ", s, task.wlength_str()),
            None => " ".to_string(),
        };
        if let Some(norm) = task.norm {
            header.push_str(&format!("n{}_{}mer{}", norm, task.wlength, str_subtr));
        }
        header.push_str(task.task_label());
    }
    tsv_buf.push_str(&header);
    tsv_buf.push('\n');

    // bedGraph track header
    let first = &tasklist[0];
    let first_subtr = match &first.subtr {
        Some(s) => format!("-{}_{} ", s, first.wlength_str()),
        None => " ".to_string(),
    };
    let mut metric_name = String::new();
    if let Some(norm) = first.norm {
        metric_name.push_str(&format!("n{}_{}mer{}", norm, first.wlength, first_subtr));
    }
    metric_name.push_str(first.task_label().trim());
    bg_buf.push_str(&format!(
        "track type=bedGraph name=\"{}_{metric_name}\" description=\"{metric_name} for {}\"\n",
        name, name
    ));

    // Pre-compute global / standard patterns for tasks that need them (sequential —
    // they depend on the full sequence and are small).
    let std_patterns: Vec<Option<Pattern>> = tasklist.iter().map(|task| {
        if matches!(task.task.as_str(), "D" | "GD" | "GPS" | "GRPS" | "GV" | "GRV") {
            let mut p = Pattern::new(task.wlength);
            p.set_pattern(seq, task.norm);
            Some(p)
        } else {
            None
        }
    }).collect();

    // Enumerate all window positions upfront so rayon can index into them.
    // Using half-open intervals [start, stop) for proper bedgraph format.
    let windows: Vec<(usize, usize)> = {
        let mut w = Vec::new();
        let mut start = 0usize;
        let mut stop = frame;  // stop is exclusive (half-open interval)

        if circular {
            // Allow windows to wrap around (circular genome)
            while stop <= seq_len + frame {
                w.push((start, stop));
                start += step;
                stop += step;
            }
        } else {
            // Only include windows that don't extend beyond sequence end (linear genome)
            while stop <= seq_len {
                w.push((start, stop));
                start += step;
                stop += step;
            }
        }
        w
    };

    // Compute all windows in parallel, preserving order via indexed collect.
    // Each element: (start, coord_stop, values_string, first_value_opt)
    let rows: Vec<(usize, usize, String, Option<f64>)> = windows
        .par_iter()
        .map(|&(start, stop)| {
            if !quiet {
                println!("[{}, {})", start, stop);
            }

            // Build locus (may wrap around circularly).
            let locus: Vec<u8> = if stop > seq_len {
                let tail = &seq[start..];
                let wrap_end = (stop - seq_len).min(seq_len);
                [tail, &seq[..wrap_end]].concat()
            } else {
                seq[start..stop].to_vec()
            };

            // coord_stop is already exclusive (stop is the exclusive end in half-open interval)
            let coord_stop = if stop > seq_len { stop - seq_len } else { stop };

            // Compute all tasks for this window (optionally in parallel for large task lists).
            let values: Vec<String> = if tasklist.len() >= 4 {
                tasklist.par_iter().enumerate().map(|(i, task)| {
                    compute_value(
                        &locus,
                        &task.task,
                        task.norm,
                        task.wlength,
                        &task.subtr,
                        std_patterns[i].as_ref(),
                    )
                }).collect()
            } else {
                tasklist.iter().enumerate().map(|(i, task)| {
                    compute_value(
                        &locus,
                        &task.task,
                        task.norm,
                        task.wlength,
                        &task.subtr,
                        std_patterns[i].as_ref(),
                    )
                }).collect()
            };

            let first_value = values.first().and_then(|v| v.parse::<f64>().ok());

            let mut row = format!("{}\t{}\t{}", name, start, coord_stop);
            for v in &values {
                let _ = write!(row, "\t{}", v);
            }

            (start, coord_stop, row, first_value)
        })
        .collect();

    // Write rows sequentially in order (already sorted by window position).
    for (start, coord_stop, row, first_value) in rows {
        tsv_buf.push_str(&row);
        tsv_buf.push('\n');
        if let Some(fv) = first_value {
            let _ = writeln!(bg_buf, "{}\t{}\t{}\t{}", name, start, coord_stop, fv);
        }
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn fasta_files_in_dir(dir: &Path) -> Vec<PathBuf> {
    let Ok(entries) = std::fs::read_dir(dir) else { return vec![] };
    let mut paths: Vec<PathBuf> = entries
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| {
            p.extension()
                .and_then(|e| e.to_str())
                .map(|e| matches!(e.to_lowercase().as_str(), "fasta" | "fa" | "fna"))
                .unwrap_or(false)
        })
        .collect();
    paths.sort();
    paths
}

fn chrono_now() -> String {
    // Simple wall-clock timestamp without extra dependencies
    use std::time::{SystemTime, UNIX_EPOCH};
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    let (h, m, s) = (secs / 3600 % 24, secs / 60 % 60, secs % 60);
    format!("{}:{}:{}", h, m, s)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::TempDir;

    fn make_fasta(dir: &TempDir, name: &str, content: &str) -> PathBuf {
        let p = dir.path().join(name);
        let mut f = std::fs::File::create(&p).unwrap();
        f.write_all(content.as_bytes()).unwrap();
        p
    }

    #[test]
    fn test_fasta_files_in_dir_finds_fasta() {
        let tmp = TempDir::new().unwrap();
        make_fasta(&tmp, "a.fasta", ">seq1\nACGT\n");
        make_fasta(&tmp, "b.txt", "not fasta");
        let files = fasta_files_in_dir(tmp.path());
        assert_eq!(files.len(), 1);
        assert!(files[0].to_str().unwrap().ends_with("a.fasta"));
    }

    #[test]
    fn test_gc_integration() {
        let tmp = TempDir::new().unwrap();
        let fpath = make_fasta(&tmp, "test.fasta", &format!(">test\n{}\n", "ACGT".repeat(2000)));
        let mut v = Validator::new();
        let (tasks, frame, step) = v.validate("GC", Some(4600), Some(2300)).unwrap();
        // Run from the tmp dir so output lands there (binary writes to cwd = tmp).
        let orig_dir = std::env::current_dir().unwrap();
        std::env::set_current_dir(tmp.path()).unwrap();
        process_file(&fpath, &tasks, frame, step, true, false).unwrap();
        std::env::set_current_dir(&orig_dir).unwrap();
        assert!(tmp.path().join("test_oligowords.tsv").exists());
        assert!(tmp.path().join("test_oligowords.bedgraph").exists());
    }
}
