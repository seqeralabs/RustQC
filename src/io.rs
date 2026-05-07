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
use std::time::Duration;

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
// Hashing helpers
// ============================================================

/// FNV-1a offset basis constant.
pub const FNV1A_OFFSET: u64 = 0xcbf2_9ce4_8422_2325;

/// FNV-1a prime constant.
const FNV1A_PRIME: u64 = 0x0000_0100_0000_01b3;

/// FNV-1a hash of a byte slice.
///
/// A fast, non-cryptographic hash used throughout RustQC for hashing
/// read names, position keys, and fragment identifiers without heap
/// allocation.
#[inline(always)]
pub fn fnv1a(bytes: &[u8]) -> u64 {
    let mut hash = FNV1A_OFFSET;
    fnv1a_update(&mut hash, bytes);
    hash
}

/// FNV-1a streaming update: mix additional bytes into an existing hash state.
#[inline(always)]
pub fn fnv1a_update(hash: &mut u64, bytes: &[u8]) {
    for &b in bytes {
        *hash ^= b as u64;
        *hash = hash.wrapping_mul(FNV1A_PRIME);
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

/// Format a count with SI suffixes (e.g. "1.5K", "48.2M", "2.3G").
///
/// Used for compact human-readable counts in progress messages and summaries.
pub fn format_count(n: u64) -> String {
    use number_prefix::NumberPrefix;
    match NumberPrefix::decimal(n as f64) {
        NumberPrefix::Standalone(n) => format!("{n}"),
        NumberPrefix::Prefixed(prefix, n) => {
            // Map SI prefixes to short single-char suffixes
            let suffix = match prefix {
                number_prefix::Prefix::Kilo => "K",
                number_prefix::Prefix::Mega => "M",
                number_prefix::Prefix::Giga => "G",
                number_prefix::Prefix::Tera => "T",
                _ => return format!("{:.1}{prefix:?}", n),
            };
            format!("{n:.1}{suffix}")
        }
    }
}

/// Format a percentage string (e.g. "(83.3%)").
pub fn format_pct(n: u64, total: u64) -> String {
    if total == 0 {
        return "(0.0%)".to_string();
    }
    format!("({:.1}%)", n as f64 / total as f64 * 100.0)
}

/// Format a duration as human-friendly mm:ss or h:mm:ss.
///
/// - Under 60s: `"45.2s"`
/// - Under 1h: `"1:23"`
/// - Over 1h: `"1:02:34"`
pub fn format_duration(d: Duration) -> String {
    let total_secs = d.as_secs_f64();
    if total_secs < 60.0 {
        return format!("{total_secs:.1}s");
    }
    let total_secs = d.as_secs();
    let hours = total_secs / 3600;
    let minutes = (total_secs % 3600) / 60;
    let seconds = total_secs % 60;
    if hours > 0 {
        format!("{hours}:{minutes:02}:{seconds:02}")
    } else {
        format!("{minutes}:{seconds:02}")
    }
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
    fn test_fnv1a_deterministic() {
        // Same input always produces the same hash
        let hash1 = fnv1a(b"test_read_name");
        let hash2 = fnv1a(b"test_read_name");
        assert_eq!(hash1, hash2);
        // Different inputs produce different hashes
        let hash3 = fnv1a(b"other_read_name");
        assert_ne!(hash1, hash3);
    }

    #[test]
    fn test_fnv1a_streaming_matches_oneshot() {
        // Streaming update should match one-shot hash
        let oneshot = fnv1a(b"hello world");
        let mut streaming = FNV1A_OFFSET;
        fnv1a_update(&mut streaming, b"hello ");
        fnv1a_update(&mut streaming, b"world");
        assert_eq!(oneshot, streaming);
    }

    #[test]
    fn test_format_with_commas() {
        assert_eq!(format_with_commas(0), "0");
        assert_eq!(format_with_commas(999), "999");
        assert_eq!(format_with_commas(1000), "1,000");
        assert_eq!(format_with_commas(1234567), "1,234,567");
    }

    #[test]
    fn test_format_count_small() {
        assert_eq!(format_count(0), "0");
        assert_eq!(format_count(42), "42");
        assert_eq!(format_count(999), "999");
    }

    #[test]
    fn test_format_count_thousands() {
        assert_eq!(format_count(1000), "1.0K");
        assert_eq!(format_count(1500), "1.5K");
        assert_eq!(format_count(50000), "50.0K");
    }

    #[test]
    fn test_format_count_millions() {
        assert_eq!(format_count(1_000_000), "1.0M");
        assert_eq!(format_count(48_200_000), "48.2M");
        assert_eq!(format_count(50_000_000), "50.0M");
    }

    #[test]
    fn test_format_count_billions() {
        assert_eq!(format_count(1_000_000_000), "1.0G");
        assert_eq!(format_count(5_000_000_000), "5.0G");
    }

    #[test]
    fn test_format_pct() {
        assert_eq!(format_pct(833, 1000), "(83.3%)");
        assert_eq!(format_pct(0, 0), "(0.0%)");
        assert_eq!(format_pct(1000, 1000), "(100.0%)");
    }

    #[test]
    fn test_format_duration_seconds() {
        assert_eq!(format_duration(Duration::from_secs_f64(0.5)), "0.5s");
        assert_eq!(format_duration(Duration::from_secs_f64(45.2)), "45.2s");
        assert_eq!(format_duration(Duration::from_secs_f64(59.9)), "59.9s");
    }

    #[test]
    fn test_format_duration_minutes() {
        assert_eq!(format_duration(Duration::from_secs(60)), "1:00");
        assert_eq!(format_duration(Duration::from_secs(83)), "1:23");
        assert_eq!(format_duration(Duration::from_secs(3599)), "59:59");
    }

    #[test]
    fn test_format_duration_hours() {
        assert_eq!(format_duration(Duration::from_secs(3600)), "1:00:00");
        assert_eq!(format_duration(Duration::from_secs(3754)), "1:02:34");
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
