use crate::kmer::generate_word_list;
use crate::pattern::Pattern;

/// Entry in the ranked word list: [word, original_index, deviation, rank_or_freq]
#[derive(Clone, Debug)]
pub struct WordEntry {
    pub word: String,
    pub index: usize,
    pub deviation: f64,
    pub rank: f64, // used as rank (usize cast to f64) after ranking, or word_frequency before
}

/// Ranked word-frequency list derived from a Pattern.
pub struct WordList {
    #[allow(dead_code)]
    pub wlength: usize,
    entries: Vec<WordEntry>,
}

#[allow(dead_code)]

impl WordList {
    pub fn new(wlength: usize) -> Self {
        let word_list = generate_word_list(wlength);
        let entries = word_list
            .into_iter()
            .enumerate()
            .map(|(i, w)| WordEntry { word: w, index: i, deviation: 0.0, rank: 0.0 })
            .collect();
        WordList { wlength, entries }
    }

    /// Compute deviation values from the pattern.
    pub fn set_value(&mut self, pattern: &Pattern) {
        let normal_expect = pattern.mean_word_number();
        for entry in &mut self.entries {
            let value = pattern.real_number(&entry.word);
            entry.rank = pattern.word_frequency(&entry.word);
            if pattern.pat_type != 0 {
                let expect = pattern.normalized_value(&entry.word).max(1.0);
                entry.deviation = (value - expect) / normal_expect;
            } else if normal_expect != 0.0 {
                entry.deviation = (value - normal_expect) / normal_expect;
            } else {
                entry.deviation = 0.0;
            }
        }
    }

    pub fn entries(&self) -> &[WordEntry] {
        &self.entries
    }

    pub fn entries_mut(&mut self) -> &mut Vec<WordEntry> {
        &mut self.entries
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Compute variance of deviations, optionally normalised.
    pub fn get_variance(&self, pattern: Option<&Pattern>) -> f64 {
        let n = self.entries.len() as f64;
        if n <= 1.0 {
            return 0.0;
        }
        let total: f64 = self.entries.iter().map(|e| e.deviation).sum();
        let total_sq: f64 = self.entries.iter().map(|e| e.deviation * e.deviation).sum();
        let var = (total_sq - total * total / n) / (n - 1.0);
        let norm = if let Some(p) = pattern {
            get_norm_stdev(p.word_length(), p.seq_length())
        } else {
            1.0
        };
        var / norm
    }
}

fn get_norm_stdev(_word_length: usize, _seq_length: usize) -> f64 {
    1.0
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_length_2mer() {
        let wl = WordList::new(2);
        assert_eq!(wl.len(), 16);
    }

    #[test]
    fn test_wlength_set() {
        let wl = WordList::new(2);
        assert_eq!(wl.wlength, 2);
    }

    #[test]
    fn test_variance_non_negative() {
        let mut wl = WordList::new(2);
        let mut p = Pattern::new(2);
        p.set_pattern(b"ACGTACGT", None);
        wl.set_value(&p);
        let v = wl.get_variance(None);
        assert!(v >= 0.0);
    }
}
