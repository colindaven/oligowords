use std::collections::HashMap;

/// Markov-chain normalization table built from a word-frequency dictionary.
///
/// Mirrors the Python `NormalizationTable` class.
#[derive(Clone)]
pub struct NormalizationTable {
    words: HashMap<String, f64>,
    frame_len: usize,
}

impl NormalizationTable {
    pub fn new(words: HashMap<String, f64>) -> Self {
        let frame_len = words.keys().next().map(|k| k.len()).unwrap_or(0);
        NormalizationTable { words, frame_len }
    }

    /// Compute the Markov likelihood of `word` using the stored frequency table.
    pub fn word_likelihood(&self, word: &str) -> f64 {
        if self.frame_len == 0 || word.is_empty() {
            return 0.0;
        }
        // Probability of the first frame-length prefix
        let prefix = &word[..self.frame_len.min(word.len())];
        let mut p_val = *self.words.get(prefix).unwrap_or(&0.0);

        if word.len() <= self.frame_len {
            return p_val;
        }

        let mut start = 1usize;
        let mut stop = self.frame_len + 1;
        while stop <= word.len() {
            let subword = &word[start..stop - 1]; // context = word[start..stop-1], matches Python word[start:stop-1]
            let subset = self.set_subset(subword);
            let kmer = &word[start..stop];
            p_val *= subset.get(kmer).copied().unwrap_or(0.0);
            start += 1;
            stop += 1;
        }
        p_val
    }

    fn set_subset(&self, subword: &str) -> HashMap<String, f64> {
        let sw_len = subword.len();
        let mut total = 0.0f64;
        let mut subset = HashMap::new();
        for (word, &freq) in &self.words {
            if word.starts_with(subword) && word.len() == sw_len + 1 {
                subset.insert(word.clone(), freq);
                total += freq;
            }
        }
        if total > 0.0 {
            for v in subset.values_mut() {
                *v /= total;
            }
        }
        subset
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use crate::pattern::Pattern;

    #[test]
    fn test_word_likelihood_in_range() {
        let seq = b"ACGTACGT".repeat(10);
        let mut _p = Pattern::new(3);
        _p.set_pattern(&seq, Some(0));
        // With normalization=2 (2-mer table):
        let mut p2 = Pattern::new(3);
        p2.set_pattern(&seq, Some(2));
        for word in p2.word_list() {
            let nt = p2.normalization_table().unwrap();
            let l = nt.word_likelihood(&word);
            assert!(l >= 0.0 && l <= 1.0, "likelihood out of range: {}", l);
        }
    }
}
