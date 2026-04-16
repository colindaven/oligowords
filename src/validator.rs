use std::collections::HashMap;
use anyhow::{anyhow, Result};

// ---------------------------------------------------------------------------
// Minimal word-length to sliding-window size constraints
// ---------------------------------------------------------------------------

fn minimal_length(wlength: usize) -> Option<usize> {
    let table: HashMap<usize, usize> = [
        (2, 300), (3, 1200), (4, 4600), (5, 18500), (6, 74000), (7, 295000),
    ].iter().cloned().collect();
    table.get(&wlength).copied()
}

// ---------------------------------------------------------------------------
// TaskSpec
// ---------------------------------------------------------------------------

/// A parsed, validated task specification.
#[derive(Debug, Clone)]
pub struct TaskSpec {
    pub task: String,
    pub norm: Option<usize>,
    pub wlength: usize,
    pub subtr: Option<String>,
    #[allow(dead_code)]
    pub id: String,
}

impl TaskSpec {
    /// Human-readable label matching the Python `Validator.tasks` dict values.
    pub fn task_label(&self) -> &'static str {
        match self.task.as_str() {
            "GC"  => "GC-content ",
            "GCS" => "G/C-skew ",
            "ATS" => "A/T-skew ",
            "D"   => "Distance ",
            "GD"  => "Global distance ",
            "PS"  => "Pattern skew ",
            "RPS" => "Relative pattern skew ",
            "GPS" => "Global pattern skew ",
            "GRPS"=> "Global relative pattern skew ",
            "V"   => "Variance ",
            "GV"  => "Global variance ",
            "RV"  => "Relative variance ",
            "GRV" => "Global relative variance ",
            _     => "Unknown ",
        }
    }

    pub fn wlength_str(&self) -> String {
        format!("{}mer", self.wlength)
    }
}

// ---------------------------------------------------------------------------
// Validator
// ---------------------------------------------------------------------------

pub struct Validator;

impl Validator {
    pub fn new() -> Self { Validator }

    /// Parse the task string and validate frame/step.
    /// Returns (tasklist, frame, step).
    pub fn validate(
        &mut self,
        task_str: &str,
        frame_opt: Option<usize>,
        step_opt: Option<usize>,
    ) -> Result<(Vec<TaskSpec>, usize, usize)> {
        let valid_tasks = [
            "GC", "GCS", "ATS", "D", "GD", "PS", "RPS", "GPS", "GRPS",
            "V", "GV", "RV", "GRV",
        ];

        let mut tasklist = Vec::new();
        let mut length_limits = Vec::new();

        for task_id in task_str.split(';') {
            let task_id = task_id.trim();
            if task_id.is_empty() {
                continue;
            }

            let (type_part, task_code, subtr) = parse_task_id(task_id)?;

            if !valid_tasks.contains(&task_code.as_str()) {
                return Err(anyhow!(
                    "Wrong task: {task_code}. Must be one of {:?}",
                    valid_tasks
                ));
            }

            let (norm, wlength) = if let Some(ref tp) = type_part {
                parse_type_part(tp)?
            } else {
                (None, 4)
            };

            let min_len = minimal_length(wlength)
                .ok_or_else(|| anyhow!("Unsupported k-mer length: {}", wlength))?;
            length_limits.push(min_len);

            tasklist.push(TaskSpec {
                task: task_code,
                norm,
                wlength,
                subtr,
                id: task_id.to_string(),
            });
        }

        if tasklist.is_empty() {
            return Err(anyhow!("No valid tasks specified"));
        }

        let min_required = *length_limits.iter().max().unwrap();

        let frame = match frame_opt {
            Some(f) => {
                if f < min_required {
                    return Err(anyhow!(
                        "The sliding window size must be at least {} bp.", min_required
                    ));
                }
                f
            }
            None => min_required,
        };

        let step = match step_opt {
            Some(s) => {
                if s == 0 {
                    return Err(anyhow!("The step must be a positive integer."));
                }
                if s > frame {
                    return Err(anyhow!("The step {} bp is longer than the window.", s));
                }
                s
            }
            None => frame / 2,
        };

        Ok((tasklist, frame, step))
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Splits "n0_4mer:D" into (Some("n0_4mer"), "D", None).
/// Splits "GC" into (None, "GC", None).
/// Also handles subtraction like "n0_4mer:D-n1_4mer".
fn parse_task_id(s: &str) -> Result<(Option<String>, String, Option<String>)> {
    let parts: Vec<&str> = s.splitn(2, ':').collect();
    let (type_part, rest) = if parts.len() == 2 {
        (Some(parts[0].to_string()), parts[1])
    } else {
        (None, parts[0])
    };

    // rest may be "D-n1_4mer" (subtraction)
    let super_parts: Vec<&str> = rest.splitn(2, '-').collect();
    let task_code = super_parts[0].to_string();
    let subtr = if super_parts.len() == 2 {
        Some(super_parts[1].to_string())
    } else {
        None
    };

    Ok((type_part, task_code, subtr))
}

/// Parses "n0_4mer" → (Some(0), 4).
fn parse_type_part(s: &str) -> Result<(Option<usize>, usize)> {
    let parts: Vec<&str> = s.splitn(2, '_').collect();
    if parts.len() != 2 {
        return Err(anyhow!("Wrong type: {s}. Must be like 'n0_4mer'"));
    }
    let norm_str = parts[0];
    let wlen_str = parts[1];

    if !norm_str.starts_with('n') || norm_str.len() < 2 {
        return Err(anyhow!("Wrong type: {s}. Must be like 'n0_4mer'"));
    }
    let norm: usize = norm_str[1..].parse()
        .map_err(|_| anyhow!("Wrong type: {s}"))?;

    if !wlen_str.ends_with("mer") || wlen_str.len() != 4 {
        return Err(anyhow!("Wrong type: {s}. Must be like 'n0_4mer'"));
    }
    let wlength: usize = wlen_str[..1].parse()
        .map_err(|_| anyhow!("Wrong k-mer length in: {s}"))?;

    Ok((Some(norm), wlength))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_valid_task() {
        let mut v = Validator::new();
        let (tasks, frame, step) = v.validate("n0_4mer:D", Some(8000), Some(2000)).unwrap();
        assert_eq!(tasks.len(), 1);
        assert_eq!(tasks[0].task, "D");
        assert_eq!(tasks[0].wlength, 4);
        assert_eq!(frame, 8000);
        assert_eq!(step, 2000);
    }

    #[test]
    fn test_invalid_task_errors() {
        let mut v = Validator::new();
        assert!(v.validate("INVALID", Some(8000), Some(2000)).is_err());
    }

    #[test]
    fn test_frame_too_small_errors() {
        let mut v = Validator::new();
        assert!(v.validate("n0_4mer:D", Some(10), Some(5)).is_err());
    }

    #[test]
    fn test_step_zero_errors() {
        let mut v = Validator::new();
        assert!(v.validate("n0_4mer:D", Some(8000), Some(0)).is_err());
    }

    #[test]
    fn test_step_larger_than_frame_errors() {
        let mut v = Validator::new();
        assert!(v.validate("n0_4mer:D", Some(5000), Some(8000)).is_err());
    }

    #[test]
    fn test_bare_gc_task() {
        let mut v = Validator::new();
        let (tasks, _frame, _step) = v.validate("GC", Some(4600), Some(2300)).unwrap();
        assert_eq!(tasks[0].task, "GC");
        assert_eq!(tasks[0].norm, None);
        assert_eq!(tasks[0].wlength, 4);
    }

    #[test]
    fn test_multiple_tasks() {
        let mut v = Validator::new();
        let (tasks, _, _) = v.validate("n0_4mer:D;n1_4mer:V", Some(8000), Some(2000)).unwrap();
        assert_eq!(tasks.len(), 2);
    }

    #[test]
    fn test_default_frame_and_step() {
        let mut v = Validator::new();
        let (_, frame, step) = v.validate("n0_4mer:D", None, None).unwrap();
        assert_eq!(frame, 4600);
        assert_eq!(step, 2300);
    }
}
