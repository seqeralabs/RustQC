//! Shared I/O and numeric utilities.
//!
//! Provides [`open_reader`] which returns a buffered reader that transparently
//! handles `.gz` (gzip) compressed files by detecting the gzip magic bytes.
//! Also provides [`median`] for computing the median of an `f64` slice.

use anyhow::{Context, Result};
use flate2::read::GzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
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

    // Re-open the file so the reader starts from byte 0
    let file =
        File::open(path).with_context(|| format!("Failed to open file: {}", path.display()))?;

    if bytes_read >= 2 && magic == GZIP_MAGIC {
        log::debug!("Detected gzip compression: {}", path.display());
        let decoder = GzDecoder::new(file);
        Ok(Box::new(BufReader::new(decoder)))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Read the entire contents of a file into a `String`, transparently
/// decompressing gzip if the magic bytes are present.
///
/// This is the gzip-aware equivalent of [`std::fs::read_to_string`].
///
/// # Arguments
/// * `path` - Path to the file (plain-text or gzip-compressed)
///
/// # Returns
/// The full file contents as a `String`.
pub fn read_to_string<P: AsRef<Path>>(path: P) -> Result<String> {
    let path = path.as_ref();
    let mut reader =
        open_reader(path).with_context(|| format!("Failed to open file: {}", path.display()))?;
    let mut contents = String::new();
    reader
        .read_to_string(&mut contents)
        .with_context(|| format!("Failed to read file: {}", path.display()))?;
    Ok(contents)
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
    fn test_read_to_string_plain() {
        let content = "hello world\n";
        let path = write_temp_plain("rustqc_test_io_rts_plain.txt", content);
        let result = read_to_string(&path).unwrap();
        assert_eq!(result, content);
    }

    #[test]
    fn test_read_to_string_gzip() {
        let content = "hello world\n";
        let path = write_temp_gz("rustqc_test_io_rts_gzip.txt.gz", content);
        let result = read_to_string(&path).unwrap();
        assert_eq!(result, content);
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
