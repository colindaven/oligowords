/// Generates all 4^k k-mers in the same canonical order as the Python PatternMethod (BOXAGCT).
///
/// The BOXAGCT enumeration is a stateful cursor algorithm.  The cursor starts
/// at the rightmost position and moves leftward.  When it goes negative it
/// wraps around using Python-style negative indexing (cursor -1 → last char,
/// -2 → second-to-last, etc.).  The character replacement also uses 'a' for
/// non-negative positions and 'b' for negative ones, mirroring Python's
/// `_changeLetters` implementation.

pub fn generate_word_list(wlength: usize) -> Vec<String> {
    let total = 4usize.pow(wlength as u32);
    let mut word = vec![b'A'; wlength];
    let mut list = Vec::with_capacity(total);
    list.push(bytes_to_string(&word));
    for _ in 1..total {
        advance_boxagct(&mut word);
        list.push(bytes_to_string(&word));
    }
    list
}

fn bytes_to_string(word: &[u8]) -> String {
    String::from_utf8(word.to_vec()).unwrap()
}

/// Set word[cursor] using Python negative-index semantics.
/// Positive cursor → use char `a`; negative cursor → use char `b`.
fn change_letter(word: &mut [u8], a: u8, b: u8, cursor: isize) {
    let k = word.len() as isize;
    if cursor >= 0 {
        word[cursor as usize] = a;
    } else {
        word[(k + cursor) as usize] = b;
    }
}

/// Advance the word by one step using the BOXAGCT rule.
/// Mirrors Python's `getNextWord` cursor logic exactly.
fn advance_boxagct(word: &mut [u8]) {
    let k = word.len() as isize;
    let mut cursor: isize = k - 1;
    let mut running = true;

    while running {
        let idx = if cursor >= 0 { cursor as usize } else { (k + cursor) as usize };

        match word[idx] {
            b'A' => {
                change_letter(word, b'G', b'C', cursor);
                running = false;
            }
            b'G' => {
                change_letter(word, b'A', b'C', cursor);
                if cursor >= 0 {
                    cursor -= 1;
                    // flg_state stays 1 (keep iterating)
                } else {
                    running = false;
                }
            }
            b'C' => {
                change_letter(word, b'T', b'A', cursor);
                if cursor >= 0 {
                    running = false;
                } else {
                    cursor -= 1;
                }
            }
            b'T' => {
                change_letter(word, b'C', b'A', cursor);
                cursor -= 1;
            }
            _ => {
                running = false;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_word_list_length_2mer() {
        let wl = generate_word_list(2);
        assert_eq!(wl.len(), 16);
    }

    #[test]
    fn test_word_list_length_4mer() {
        let wl = generate_word_list(4);
        assert_eq!(wl.len(), 256);
    }

    #[test]
    fn test_first_word_is_all_a() {
        let wl = generate_word_list(3);
        assert_eq!(wl[0], "AAA");
    }

    #[test]
    fn test_all_words_are_unique() {
        let wl = generate_word_list(3);
        let set: std::collections::HashSet<_> = wl.iter().collect();
        assert_eq!(set.len(), wl.len());
    }

    #[test]
    fn test_words_only_acgt() {
        for word in generate_word_list(2) {
            for ch in word.chars() {
                assert!(matches!(ch, 'A' | 'C' | 'G' | 'T'));
            }
        }
    }

    #[test]
    fn test_contains_aaa() {
        assert!(generate_word_list(3).contains(&"AAA".to_string()));
    }
}
