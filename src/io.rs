//! Shared I/O and numeric utilities.
//!
//! Provides [`open_reader`] which returns a buffered reader that transparently
//! handles `.gz` (gzip) compressed files by detecting the gzip magic bytes.
//! Also provides [`median`] for computing the median of an `f64` slice.

use anyhow::{Context, Result};
use flate2::read::GzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek};
use std::path::Path;

/// Gzip magic bytes: the first two bytes of any gzip-compressed file.
const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

/// Open a file and return a buffered reader, transparently decompressing gzip.
///
/// Detection is based on the file's magic bytes (first two bytes), not the
/// file extension. This means renamed `.gz` files without the extension are
/// still handled correctly.
///
/// # Arguments
/// * `path` - Path to the file (plain-text or gzip-compressed)
///
/// # Returns
/// A boxed [`BufRead`] implementation suitable for line-by-line iteration.
pub fn open_reader<P: AsRef<Path>>(path: P) -> Result<Box<dyn BufRead>> {
    let path = path.as_ref();
    let mut file =
        File::open(path).with_context(|| format!("Failed to open file: {}", path.display()))?;

    // Read the first two bytes to check for gzip magic number
    let mut magic = [0u8; 2];
    let bytes_read = file
        .read(&mut magic)
        .with_context(|| format!("Failed to read from file: {}", path.display()))?;

    // Seek back to the beginning so the reader starts from byte 0
    file.seek(std::io::SeekFrom::Start(0))
        .with_context(|| format!("Failed to seek in file: {}", path.display()))?;

    if bytes_read >= 2 && magic == GZIP_MAGIC {
        log::debug!("Detected gzip compression: {}", path.display());
        let decoder = GzDecoder::new(file);
        Ok(Box::new(BufReader::new(decoder)))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

// ============================================================
// Formatting helpers
// ============================================================

/// Format a `u64` with comma-separated thousands (e.g. `1,234,567`).
pub fn format_with_commas(n: u64) -> String {
    let s = n.to_string();
    let bytes = s.as_bytes();
    let len = bytes.len();
    if len <= 3 {
        return s;
    }
    let mut result = String::with_capacity(len + (len - 1) / 3);
    for (i, &b) in bytes.iter().enumerate() {
        if i > 0 && (len - i).is_multiple_of(3) {
            result.push(',');
        }
        result.push(b as char);
    }
    result
}

// ============================================================
// Numeric helpers
// ============================================================

/// Compute the median of an `f64` slice.
///
/// Returns `0.0` for an empty slice.  A sorted copy is made internally so
/// the input is not modified.
///
/// # Arguments
/// * `values` - Slice of `f64` values.
pub fn median(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mut sorted: Vec<f64> = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = sorted.len();
    if n.is_multiple_of(2) {
        (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
    } else {
        sorted[n / 2]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::io::Write;

    /// Helper to create a temporary plain-text file.
    fn write_temp_plain(name: &str, content: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir();
        let path = dir.join(name);
        std::fs::write(&path, content).unwrap();
        path
    }

    /// Helper to create a temporary gzip-compressed file.
    fn write_temp_gz(name: &str, content: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir();
        let path = dir.join(name);
        let file = File::create(&path).unwrap();
        let mut encoder = GzEncoder::new(file, Compression::default());
        encoder.write_all(content.as_bytes()).unwrap();
        encoder.finish().unwrap();
        path
    }

    #[test]
    fn test_format_with_commas() {
        assert_eq!(format_with_commas(0), "0");
        assert_eq!(format_with_commas(999), "999");
        assert_eq!(format_with_commas(1000), "1,000");
        assert_eq!(format_with_commas(1234567), "1,234,567");
    }

    #[test]
    fn test_open_reader_plain() {
        let content = "line1\nline2\nline3\n";
        let path = write_temp_plain("rustqc_test_io_plain.txt", content);
        let reader = open_reader(&path).unwrap();
        let lines: Vec<String> = reader.lines().map(|l| l.unwrap()).collect();
        assert_eq!(lines, vec!["line1", "line2", "line3"]);
    }

    #[test]
    fn test_open_reader_gzip() {
        let content = "line1\nline2\nline3\n";
        let path = write_temp_gz("rustqc_test_io_gzip.txt.gz", content);
        let reader = open_reader(&path).unwrap();
        let lines: Vec<String> = reader.lines().map(|l| l.unwrap()).collect();
        assert_eq!(lines, vec!["line1", "line2", "line3"]);
    }

    #[test]
    fn test_open_reader_empty_file() {
        let path = write_temp_plain("rustqc_test_io_empty.txt", "");
        let reader = open_reader(&path).unwrap();
        let lines: Vec<String> = reader.lines().map(|l| l.unwrap()).collect();
        assert!(lines.is_empty());
    }

    #[test]
    fn test_open_reader_nonexistent() {
        let result = open_reader("/nonexistent/path/file.txt");
        assert!(result.is_err());
    }
}
