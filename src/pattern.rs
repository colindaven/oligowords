use std::collections::HashMap;
use crate::kmer::generate_word_list;
use crate::normalization::NormalizationTable;

// ---------------------------------------------------------------------------
// Pattern
// ---------------------------------------------------------------------------

/// Holds k-mer counts and statistics for a DNA sequence window.
///
/// Mirrors the Python `Pattern` class.
#[derive(Clone)]
pub struct Pattern {
    pub wlength: usize,
    pub total_word_number: usize,
    /// word → [real_count, expected_count]
    pub words: HashMap<String, [f64; 2]>,
    /// nucleotide frequencies A/C/G/T
    pub boxagct: HashMap<char, f64>,
    pub pat_type: usize,
    pub norm_table: Option<NormalizationTable>,
}

impl Pattern {
    pub fn new(wlength: usize) -> Self {
        let word_list = generate_word_list(wlength);
        let mut words = HashMap::with_capacity(word_list.len());
        for w in word_list {
            words.insert(w, [0.0, 0.0]);
        }
        let mut boxagct = HashMap::new();
        for c in ['A', 'C', 'G', 'T'] {
            boxagct.insert(c, 0.0);
        }
        Pattern {
            wlength,
            total_word_number: 0,
            words,
            boxagct,
            pat_type: 0,
            norm_table: None,
        }
    }

    // ------------------------------------------------------------------
    // Setters
    // ------------------------------------------------------------------

    /// Compute counts and statistics from `seq`, optionally with normalisation.
    ///
    /// * `normalization = None` → uniform expectation
    /// * `normalization = Some(0)` → same as None (no normalisation)
    /// * `normalization = Some(1)` → 1st-order Markov (mono-nuc frequencies)
    /// * `normalization = Some(n)` → (n-1)-mer frequency table
    pub fn set_pattern(&mut self, seq: &[u8], normalization: Option<usize>) {
        self.set_matrix(seq);
        self.set_word_statistics(seq);

        match normalization {
            None | Some(0) => {
                self.norm_table = None;
                self.set_expectation();
            }
            Some(1) => {
                // Use mononucleotide frequencies
                let word_freqs: HashMap<String, f64> = self.boxagct
                    .iter()
                    .map(|(&c, &f)| (c.to_string(), f))
                    .collect();
                let nt = NormalizationTable::new(word_freqs);
                self.norm_table = Some(nt);
                self.pat_type = 1;
                self.set_expectation();
            }
            Some(n) => {
                let mut norm_pat = Pattern::new(n);
                norm_pat.set_pattern(seq, None);
                let word_freqs: HashMap<String, f64> = norm_pat.words.keys()
                    .map(|w| (w.clone(), norm_pat.word_frequency(w)))
                    .collect();
                let nt = NormalizationTable::new(word_freqs);
                self.norm_table = Some(nt);
                self.pat_type = n;
                self.set_expectation();
            }
        }
    }

    fn set_matrix(&mut self, seq: &[u8]) {
        let n = seq.len() as f64;
        if n == 0.0 {
            return;
        }
        for c in ['A', 'C', 'G', 'T'] {
            let count = seq.iter().filter(|&&b| b == c as u8).count();
            self.boxagct.insert(c, count as f64 / n);
        }
    }

    fn set_word_statistics(&mut self, seq: &[u8]) {
        self.total_word_number = seq.len().saturating_sub(self.wlength);
        // Reset counts
        for v in self.words.values_mut() {
            v[0] = 0.0;
        }
        for i in 0..self.total_word_number {
            let slice = &seq[i..i + self.wlength];
            if let Ok(word) = std::str::from_utf8(slice) {
                if let Some(v) = self.words.get_mut(word) {
                    v[0] += 1.0;
                }
            }
        }
    }

    fn set_expectation(&mut self) {
        let n_words = self.words.len() as f64;
        let total = self.total_word_number as f64;
        let uniform_p = if total > 0.0 { n_words / total } else { 0.0 };

        let keys: Vec<String> = self.words.keys().cloned().collect();
        for word in keys {
            let p_val = if let Some(ref nt) = self.norm_table {
                nt.word_likelihood(&word)
            } else {
                uniform_p
            };
            self.words.get_mut(&word).unwrap()[1] = p_val * total;
        }
    }

    pub fn set_normalization_table(&mut self, table: NormalizationTable) {
        self.norm_table = Some(table);
        self.set_expectation();
    }

    // ------------------------------------------------------------------
    // Getters
    // ------------------------------------------------------------------

    pub fn word_length(&self) -> usize { self.wlength }
    pub fn seq_length(&self) -> usize { self.total_word_number + self.wlength }
    #[allow(dead_code)]
    pub fn word_list(&self) -> Vec<String> { self.words.keys().cloned().collect() }
    pub fn real_number(&self, word: &str) -> f64 { self.words.get(word).map(|v| v[0]).unwrap_or(0.0) }
    pub fn normalized_value(&self, word: &str) -> f64 { self.words.get(word).map(|v| v[1]).unwrap_or(0.0) }
    pub fn normalization_table(&self) -> Option<&NormalizationTable> { self.norm_table.as_ref() }
    pub fn mean_word_number(&self) -> f64 {
        if self.words.is_empty() { 0.0 } else { self.total_word_number as f64 / self.words.len() as f64 }
    }
    #[allow(dead_code)]
    pub fn nucleotide_frequency(&self, letter: char) -> f64 {
        *self.boxagct.get(&letter).unwrap_or(&0.0)
    }
    pub fn word_frequency(&self, word: &str) -> f64 {
        if self.total_word_number == 0 { return 0.0; }
        self.words.get(word).map(|v| v[0] / self.total_word_number as f64).unwrap_or(0.0)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_word_length() {
        let p = Pattern::new(3);
        assert_eq!(p.word_length(), 3);
    }

    #[test]
    fn test_word_list_length_2mer() {
        let mut p = Pattern::new(2);
        p.set_pattern(b"ACGTACGT", None);
        assert_eq!(p.word_list().len(), 16);
    }

    #[test]
    fn test_total_word_number() {
        let mut p = Pattern::new(2);
        p.set_pattern(b"ACGTACGT", None);
        // len=8, wlength=2 → 8-2=6
        assert_eq!(p.total_word_number, 6);
    }

    #[test]
    fn test_aa_count() {
        let mut p = Pattern::new(2);
        p.set_pattern(b"AAAA", None);
        // AAAA has 2 2-mers of length 2 (total_word_number = 4-2=2), both "AA"
        assert_eq!(p.real_number("AA"), 2.0);
    }

    #[test]
    fn test_nucleotide_frequency_all_a() {
        let mut p = Pattern::new(3);
        p.set_pattern(b"AAAA", None);
        assert!((p.nucleotide_frequency('A') - 1.0).abs() < 1e-9);
        assert_eq!(p.nucleotide_frequency('C'), 0.0);
    }

    #[test]
    fn test_word_freq_sum_approx_1() {
        let seq = b"ACGTACGT".repeat(10);
        let mut p = Pattern::new(2);
        p.set_pattern(&seq, None);
        let total: f64 = p.word_list().iter().map(|w| p.word_frequency(w)).sum();
        assert!((total - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_pat_type_set_with_normalization() {
        let mut p = Pattern::new(2);
        p.set_pattern(b"ACGTACGT", Some(0));
        assert_eq!(p.pat_type, 0);
    }

    #[test]
    fn test_4mer_acgt_appears() {
        let mut p = Pattern::new(4);
        p.set_pattern(b"ACGTACGT", None);
        assert!(p.real_number("ACGT") >= 1.0);
    }

    #[test]
    fn test_copy_via_clone() {
        let mut p1 = Pattern::new(2);
        p1.set_pattern(b"ACGT", None);
        let p2 = p1.clone();
        assert_eq!(p2.word_length(), p1.word_length());
    }
}
