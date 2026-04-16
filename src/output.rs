use std::path::Path;
use anyhow::Result;

/// Write the accumulated TSV and bedGraph buffers to disk.
pub fn write_outputs(tsv_path: &Path, bg_path: &Path, tsv: &str, bg: &str) -> Result<()> {
    std::fs::write(tsv_path, tsv)?;
    std::fs::write(bg_path, bg)?;
    Ok(())
}
