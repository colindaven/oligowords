# OligoWords

Sliding-window oligonucleotide composition analysis of DNA sequences. For each window along a chromosome or contig, OligoWords computes one or more metrics that describe how the local k-mer frequency distribution compares to an expected or global background. The output is suitable for visualisation in a genome browser (bedGraph) or downstream statistical analysis (TSV).

# Source of original code

This was originally written by Reva and Tuemmler 2005 [link](https://link.springer.com/article/10.1186/1471-2105-6-251) but was in outdated Python2. It was rewritten with Claude AI is rust in April 2026. New features like bedgraph and bigwig functionality were added. 

# Status - alpha
Functionality is being tested in April 2026, mainly focused on the distance parameter and its use in predicting genomic islands or aberrant sequences.



---

## Table of contents

1. [Quick start (Rust)](#quick-start-rust)
2. [Installation](#installation)
   - [Rust (recommended)](#rust-recommended)
   - [Python](#python)
3. [Input](#input)
4. [Parameters](#parameters)
5. [Task syntax](#task-syntax)
6. [Available metrics](#available-metrics)
7. [Output files](#output-files)
8. [Interpretation guide](#interpretation-guide)
9. [Examples](#examples)
10. [Implementation notes](#implementation-notes)

---

## Quick start (Rust)

### Easiest - 
* Download a release using wget.

### Compile with cargo
```bash
cd rust
cargo build --release

# Single file, default tasks, default window/step
./target/release/oligowords -i genome.fasta

# GC-content only, 5 kb windows, 1 kb steps
./target/release/oligowords -i genome.fasta --task GC --frame 5000 --step 1000

# Multiple metrics, suppress progress output
./target/release/oligowords -i genome.fasta --task "n0_4mer:D;n1_4mer:V" --quiet

# Batch-process all FASTA files in a directory
./target/release/oligowords --path /data/fastas/ --task "n1_4mer:V"
```

Output files are written to the **current working directory**.

---

## Installation

### Rust (recommended)

Requires [Rust ≥ 1.75](https://rustup.rs).

```bash
cd rust
cargo build --release
# Binary at: rust/target/release/oligowords
```

The Rust implementation uses [needletail](https://github.com/onecodex/needletail) for fast FASTA parsing and [rayon](https://github.com/rayon-rs/rayon) for multi-threaded window computation. It is significantly faster than the Python version and handles large genomes with low memory overhead.

### Python

Requires Python ≥ 3.8. No external dependencies beyond the standard library.

```bash
cd python3
python3 oligowords.py -i genome.fasta
```

---

## Input

OligoWords accepts standard FASTA files (`.fasta`, `.fa`, `.fna`). Multi-record files are supported; each record is processed independently as a separate chromosome or contig. Sequences are upper-cased automatically; any characters that are not `A`, `C`, `G`, or `T` are silently discarded.

| Option | Description |
|--------|-------------|
| `-i / --input FILE` | A single FASTA file (one or more records) |
| `--path DIR` | Directory scanned recursively for `*.fasta / *.fa / *.fna` files (default: current directory) |

`-i` and `--path` are mutually exclusive.

---

## Parameters

| Option | Default | Description |
|--------|---------|-------------|
| `--task TASKS` | `n0_4mer:D;n0_4mer:PS;n1_4mer:V;n1_4mer:GV` | Semicolon-separated list of analyses (see [Task syntax](#task-syntax)) |
| `--frame BP` | Minimum valid size for chosen k-mer | Sliding-window size in base pairs |
| `--step BP` | `frame / 2` | Step between successive windows (50 % overlap by default) |
| `-q / --quiet` | off | Suppress per-window progress messages |

### Minimum window sizes

The window must contain enough k-mers for the statistics to be meaningful. The minimum enforced frame sizes are:

| k-mer | Minimum frame |
|-------|--------------|
| 2-mer | 300 bp |
| 3-mer | 1 200 bp |
| 4-mer | 4 600 bp |
| 5-mer | 18 500 bp |
| 6-mer | 74 000 bp |
| 7-mer | 295 000 bp |

---

## Task syntax

Each task follows the pattern:

```
[nN_Kmer:]TASKCODE[-nN_Kmer]
```

| Component | Meaning | Example |
|-----------|---------|---------|
| `nN` | Normalisation order (0 = none, 1 = mononucleotide, 2 = dinucleotide, …) | `n0`, `n1` |
| `Kmer` | k-mer size (`2mer`–`7mer`) | `4mer` |
| `TASKCODE` | Metric to compute (see table below) | `D`, `V`, `GC` |
| `-nN_Kmer` | Optional subtraction: subtract this order from the result | `n0_4mer:D-n1_4mer` |

The prefix `nN_Kmer:` is optional for simple metrics like `GC`, `GCS`, `ATS`.

Multiple tasks are separated by `;`:

```
n0_4mer:D;n1_4mer:V;GC
```

---

## Available metrics

### Nucleotide composition

| Code | Name | Description |
|------|------|-------------|
| `GC` | GC-content | Fraction of G+C bases in the window. Range 0–1. |
| `GCS` | G/C-skew | `(G − C) / (G + C)`. Positive on the leading strand during replication. Range −1 to +1. |
| `ATS` | A/T-skew | `(A − T) / (A + T)`. Complementary to GCS. Range −1 to +1. |

### k-mer pattern metrics

All pattern metrics compare the k-mer frequency distribution of the local window to a reference distribution. The reference is either the **whole-sequence** distribution (global tasks: `GD`, `GPS`, `GV`, `GRV`, `GRPS`) or an **independent per-window expectation** (local tasks: `D`, `PS`, `V`, `RV`, `RPS`).

| Code | Name | Description |
|------|------|-------------|
| `D` | Distance | Rank-order distance between the local k-mer distribution and the local expected distribution. High values indicate unusual composition. |
| `GD` | Global distance | As `D`, but the reference is the whole-sequence k-mer distribution. Highlights regions that differ from the genome average. |
| `PS` | Pattern skew | Rank-order asymmetry between a word and its reverse complement across the window. Detects strand-biased gene content or replication asymmetry. |
| `GPS` | Global pattern skew | As `PS` using the global reference. |
| `RPS` | Relative pattern skew | `PS` normalised for sequence length. |
| `GRPS` | Global relative pattern skew | `GPS` normalised for sequence length. |
| `V` | Variance | Variance of normalised k-mer deviations within the window. High variance indicates a heterogeneous or unusual compositional signature. |
| `GV` | Global variance | As `V` using global frequencies. |
| `RV` | Relative variance | `V` normalised for window length. |
| `GRV` | Global relative variance | `GV` normalised for window length. |

### Normalisation orders (`nN`)

| Order | Expectation model |
|-------|-------------------|
| `n0` | Uniform (all k-mers equally likely) |
| `n1` | Mononucleotide frequencies of the window |
| `n2` | Dinucleotide Markov chain |
| `n3`–`n7` | Higher-order Markov chains |

Higher normalisation orders progressively remove background compositional bias, making the metric more sensitive to *relative* departures (e.g. CpG suppression, oligonucleotide avoidance).

---

## Output files

Two files are written per input FASTA file to the **current working directory**:

### `<stem>_oligowords.tsv`

Tab-separated values, human-readable.

```
# 10:44:37
# pKLC102
# Sequence length: 103609 bp.
# Window: 4600	Step: 2300
#
# Chrom	Start	Stop	n0_4mer Distance  n1_4mer Variance  ...
pKLC102	0	4600	16.46	0.502
pKLC102	2300	6900	16.41	0.491
...
```

- Header lines begin with `#`.
- Coordinates are 0-based, half-open `[Start, Stop)`.
- One column per task, in the order specified on the command line.
- The final windows may wrap around to the start of the sequence (circular chromosomes).

### `<stem>_oligowords.bedgraph`

[UCSC bedGraph](https://genome.ucsc.edu/goldenPath/help/bedgraph.html) format, containing the **first task only**. Load directly into IGV, UCSC, or any bedGraph-compatible genome browser.

```
track type=bedGraph name="pKLC102_n1_4mer Variance" description="..."
pKLC102	0	4600	0.2893
pKLC102	2300	6900	0.2574
...
```

---

## Interpretation guide

### GC-content (`GC`)
Values near 0.5 are typical for most bacteria. Sustained deviations often mark horizontally acquired islands, prophages, or AT-rich regulatory regions.

### GC- and AT-skew (`GCS`, `ATS`)
In circular chromosomes with a single replication origin, GCS and ATS show a characteristic sinusoidal pattern that switches sign at the *oriC* (replication origin) and *ter* (replication terminus). Sharp inflection points identify these landmarks.

### Distance (`D`, `GD`)
Peaks identify windows whose k-mer composition is atypical. `GD` peaks that co-localise with `D` valleys often indicate alien DNA (genomic islands, phage) whose local composition is internally consistent but globally divergent.

### Pattern skew (`PS`, `GPS`)
Strand asymmetry in k-mer usage. In prokaryotes, elevated PS reflects the excess of highly expressed genes on the leading strand. Combined with GCS, PS asymmetry can help partition the chromosome into replichores.

### Variance (`V`, `GV`)
Measures how heterogeneously k-mers are distributed within a window. Low-complexity regions (tandem repeats, homopolymers) give very high variance. Horizontally transferred DNA often shows elevated variance under `n1` normalisation because its k-mer usage does not match the host mononucleotide background.

---

## Examples

```bash
# Default analysis (4 metrics) on a bacterial genome
./target/release/oligowords -i genome.fasta --quiet

# GC-content and GC-skew with 10 kb windows, 2 kb step
./target/release/oligowords -i genome.fasta --task "GC;GCS" --frame 10000 --step 2000

# Detect genomic islands: distance from global composition
./target/release/oligowords -i genome.fasta --task "n1_4mer:GD" --frame 4600 --quiet

# High-order normalisation to remove mononucleotide bias
./target/release/oligowords -i genome.fasta --task "n2_4mer:V" --frame 4600

# Batch process a directory of assemblies (multi-threaded)
./target/release/oligowords --path /data/genomes/ --task "n1_4mer:V;GCS" --quiet

# Python equivalent
python3 python3/oligowords.py -i genome.fasta --task "n1_4mer:V" --frame 4600 --quiet
```

---

## Implementation notes

| | Rust | Python |
|-|------|--------|
| FASTA parsing | [needletail](https://github.com/onecodex/needletail) (streaming, handles gzip/bgzip) | Custom line-by-line parser |
| Parallelism | [rayon](https://github.com/rayon-rs/rayon): windows × tasks | None |
| Dependencies | `needletail`, `clap`, `rayon`, `anyhow` | Standard library only |
| Output directory | Current working directory | Directory of the input file |

The Rust and Python implementations produce numerically identical results (differences ≤ 10⁻¹⁶, sub-machine-epsilon floating-point noise).
