use crate::pattern::Pattern;
use crate::wordlist::WordList;

/// Compute a single metric value for a sequence window.
///
/// Returns the result as a string (matching Python behaviour).
pub fn compute_value(
    seq: &[u8],
    task: &str,
    normalization: Option<usize>,
    wlength: usize,
    subtr: &Option<String>,
    std_pattern: Option<&Pattern>,
) -> String {
    match task {
        "GC" => gc_content(seq),
        "GCS" => gc_skew(seq),
        "ATS" => at_skew(seq),
        "D" | "GD" => {
            let curr_pat = make_pat(seq, wlength, normalization, std_pattern, task);
            let subtr_pat = make_subtr_pat(seq, wlength, subtr);
            let dist = pattern_distance(&curr_pat, std_pattern.unwrap());
            let s_dist = subtr_pat.as_ref()
                .map(|sp| pattern_distance(sp, std_pattern.unwrap()))
                .unwrap_or(0.0);
            format!("{:.2}", dist - s_dist)
        }
        "PS" | "GPS" | "RPS" | "GRPS" => {
            let curr_pat = make_pat(seq, wlength, normalization, std_pattern, task);
            let subtr_pat = make_subtr_pat(seq, wlength, subtr);
            let ps = oligo_skewness(&curr_pat);
            let s_ps = subtr_pat.as_ref().map(|sp| oligo_skewness(sp)).unwrap_or(0.0);
            if task == "PS" || task == "GPS" {
                format!("{:.2}", ps - s_ps)
            } else {
                let exp = expected_ps(seq.len());
                format!("{:.2}", (ps - exp[0] - 2.0 * s_ps) / exp[1])
            }
        }
        "V" | "GV" | "RV" | "GRV" => {
            let curr_pat = make_pat(seq, wlength, normalization, std_pattern, task);
            let subtr_pat = make_subtr_pat(seq, wlength, subtr);
            let mut wl = WordList::new(curr_pat.word_length());
            wl.set_value(&curr_pat);
            let variance = if task == "V" || task == "GV" {
                wl.get_variance(None)
            } else {
                wl.get_variance(Some(&curr_pat))
            };
            let s_variance = if let Some(sp) = &subtr_pat {
                let mut swl = WordList::new(sp.word_length());
                swl.set_value(sp);
                if task == "V" || task == "GV" { swl.get_variance(None) } else { swl.get_variance(Some(sp)) }
            } else {
                0.0
            };
            format!("{:.2}", variance - s_variance)
        }
        other => {
            eprintln!("Error: unknown task '{}'", other);
            "0".to_string()
        }
    }
}

// ---------------------------------------------------------------------------
// Nucleotide metrics
// ---------------------------------------------------------------------------

fn gc_content(seq: &[u8]) -> String {
    if seq.is_empty() { return "0".to_string(); }
    let gc = seq.iter().filter(|&&b| b == b'G' || b == b'C').count();
    format!("{:.2}", gc as f64 / seq.len() as f64)
}

fn gc_skew(seq: &[u8]) -> String {
    if seq.is_empty() { return "0".to_string(); }
    let g = seq.iter().filter(|&&b| b == b'G').count() as f64;
    let c = seq.iter().filter(|&&b| b == b'C').count() as f64;
    if g + c == 0.0 { return "0".to_string(); }
    format!("{:.2}", (g - c) / (g + c))
}

fn at_skew(seq: &[u8]) -> String {
    if seq.is_empty() { return "0".to_string(); }
    let a = seq.iter().filter(|&&b| b == b'A').count() as f64;
    let t = seq.iter().filter(|&&b| b == b'T').count() as f64;
    if a + t == 0.0 { return "0".to_string(); }
    format!("{:.2}", (a - t) / (a + t))
}

// ---------------------------------------------------------------------------
// Pattern helpers
// ---------------------------------------------------------------------------

fn make_pat(
    seq: &[u8],
    wlength: usize,
    normalization: Option<usize>,
    std_pattern: Option<&Pattern>,
    task: &str,
) -> Pattern {
    let mut p = Pattern::new(wlength);
    // For global tasks with normalisation, borrow the normalization table from std_pattern
    if matches!(task, "GD" | "GPS" | "GRPS" | "GV" | "GRV") && normalization.is_some() {
        p.set_pattern(seq, normalization);
        if let Some(sp) = std_pattern {
            if let Some(nt) = sp.normalization_table() {
                p.set_normalization_table(nt.clone());
            }
        }
    } else {
        p.set_pattern(seq, normalization);
    }
    p
}

fn make_subtr_pat(seq: &[u8], wlength: usize, subtr: &Option<String>) -> Option<Pattern> {
    subtr.as_ref().map(|s| {
        // subtr is like "n1" or "n1_4mer" — extract the normalization order
        let norm: usize = s.trim_start_matches('n')
            .split('_').next()
            .and_then(|n| n.parse().ok())
            .unwrap_or(0);
        let mut sp = Pattern::new(wlength);
        sp.set_pattern(seq, if norm == 0 { None } else { Some(norm) });
        sp
    })
}

// ---------------------------------------------------------------------------
// Distance / skewness
// ---------------------------------------------------------------------------

fn pattern_distance(curr_pat: &Pattern, std_pat: &Pattern) -> f64 {
    let mut swl = make_ranked_word_list(std_pat);
    let mut cwl = make_ranked_word_list(curr_pat);
    rank_word_list(&mut swl);
    rank_word_list(&mut cwl);
    get_distance(&swl, &cwl, std_pat.word_length(), true, 0.0)
}

fn oligo_skewness(curr_pat: &Pattern) -> f64 {
    let mut cwl = make_ranked_word_list(curr_pat);
    let mut cwl_compl = get_complement(&cwl);
    rank_word_list(&mut cwl);
    rank_word_list(&mut cwl_compl);
    let min_dist = strand_minimal_distance(curr_pat.word_length()) as f64;
    get_distance(&cwl, &cwl_compl, curr_pat.word_length(), false, min_dist)
}

/// A word-list entry: [word, original_index, deviation, rank_or_freq]
type WlEntry = (String, usize, f64, f64);

fn make_ranked_word_list(pat: &Pattern) -> Vec<WlEntry> {
    use crate::wordlist::WordList;
    let mut wl = WordList::new(pat.word_length());
    wl.set_value(pat);
    wl.entries().iter().map(|e| {
        (e.word.clone(), e.index, e.deviation, e.rank)
    }).collect()
}

fn rank_word_list(wl: &mut Vec<WlEntry>) {
    // Sort by deviation descending, breaking ties by word descending
    wl.sort_by(|a, b| {
        b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| b.0.cmp(&a.0))
    });
    for (j, entry) in wl.iter_mut().enumerate() {
        entry.3 = j as f64;
    }
    // Restore alphabetical order by word
    wl.sort_by(|a, b| a.0.cmp(&b.0));
}

fn get_distance(
    wl1: &[WlEntry],
    wl2: &[WlEntry],
    wlength: usize,
    best_hit: bool,
    min_dist: f64,
) -> f64 {
    if wl1.len() != wl2.len() || wl1.is_empty() {
        return 0.0;
    }
    let word_number = 4usize.pow(wlength as u32);
    let sum1: f64 = (0..word_number)
        .map(|i| (wl1[i].3 - wl2[i].3).abs())
        .sum();

    let dev = if best_hit {
        let mut wl2_compl = get_complement(wl2);
        rank_word_list(&mut wl2_compl);
        let sum2: f64 = (0..word_number)
            .map(|i| (wl1[i].3 - wl2_compl[i].3).abs())
            .sum();
        sum1.min(sum2)
    } else {
        sum1
    };

    let denom = word_number as f64 * (word_number as f64 + 1.0) / 2.0 - min_dist;
    100.0 * (dev - min_dist) / denom
}

fn strand_minimal_distance(wlength: usize) -> usize {
    if wlength % 2 == 0 {
        4usize.pow(wlength as u32) - 2usize.pow(wlength as u32)
    } else {
        4usize.pow(wlength as u32)
    }
}

/// Build the reverse-complement mapping of a word list.
fn get_complement(wl: &[WlEntry]) -> Vec<WlEntry> {
    let _n = wl.len();
    // Sort by original index
    let mut sorted = wl.to_vec();
    sorted.sort_by(|a, b| b.1.cmp(&a.1)); // reverse order by index
    let word_map: std::collections::HashMap<String, usize> =
        wl.iter().map(|e| (e.0.clone(), e.1)).collect();

    let mut result = wl.to_vec();
    for (i, entry) in result.iter_mut().enumerate() {
        let rev_word: String = sorted[i].0.chars().rev().collect();
        if let Some(&_orig_idx) = word_map.get(&rev_word) {
            entry.2 = sorted[i].2;
            entry.3 = f64::NAN; // will be overwritten by rank_word_list
        }
    }
    result
}

fn expected_ps(seq_len: usize) -> [f64; 2] {
    let n = seq_len as f64;
    [
        100.0 - 95.24 * (-2796.6 / n).exp(),
        100.0 - 98.66 * (-652.27 / n).exp(),
    ]
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gc_content_balanced() {
        let seq = b"AACCGGTT";
        let v: f64 = gc_content(seq).parse().unwrap();
        assert!((v - 0.5).abs() < 0.01);
    }

    #[test]
    fn test_gc_skew_all_g() {
        let seq = b"GGGGGGGG";
        let v: f64 = gc_skew(seq).parse().unwrap();
        assert!((v - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_at_skew_balanced() {
        let seq = b"AAAATTTT";
        let v: f64 = at_skew(seq).parse().unwrap();
        assert!(v.abs() < 0.01);
    }

    #[test]
    fn test_compute_gc_value() {
        let seq = b"AACCGGTT";
        let v = compute_value(seq, "GC", None, 4, &None, None);
        let fv: f64 = v.parse().unwrap();
        assert!((fv - 0.5).abs() < 0.01);
    }

    #[test]
    fn test_compute_gcs_value() {
        let seq = b"GGGGCCCC";
        let v = compute_value(seq, "GCS", None, 4, &None, None);
        let fv: f64 = v.parse().unwrap();
        assert!(fv.abs() < 0.01);
    }

    #[test]
    fn test_compute_ats_value_all_a() {
        let seq = b"AAAAAAAA";
        let v = compute_value(seq, "ATS", None, 4, &None, None);
        let fv: f64 = v.parse().unwrap();
        assert!((fv - 1.0).abs() < 0.01);
    }
}
