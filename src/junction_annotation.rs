//! Splice junction annotation and classification.
//!
//! Reimplementation of RSeQC's `junction_annotation.py`: extracts splice junctions
//! from BAM files (CIGAR N-operations) and classifies them as known (annotated),
//! partial novel, or complete novel by comparing against a BED12 gene model.

use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};
use log::{debug, info};
use rust_htslib::bam::{self, Read as BamRead};

// ===================================================================
// Data types
// ===================================================================

/// Classification of a splice junction.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum JunctionClass {
    /// Both splice sites match known annotation positions.
    Annotated,
    /// Exactly one splice site matches a known position.
    PartialNovel,
    /// Neither splice site matches any known position.
    CompleteNovel,
}

impl JunctionClass {
    /// Label string matching RSeQC output.
    pub fn label(self) -> &'static str {
        match self {
            JunctionClass::Annotated => "annotated",
            JunctionClass::PartialNovel => "partial_novel",
            JunctionClass::CompleteNovel => "complete_novel",
        }
    }

    /// RGB color for BED12 output.
    pub fn color(self) -> &'static str {
        match self {
            JunctionClass::Annotated => "205,0,0",
            JunctionClass::PartialNovel => "0,205,0",
            JunctionClass::CompleteNovel => "0,0,205",
        }
    }
}

/// A unique splice junction identified by chromosome and intron coordinates.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Junction {
    /// Chromosome name (uppercased internally).
    pub chrom: String,
    /// 0-based intron start (first intronic base).
    pub intron_start: u64,
    /// 0-based intron end (exclusive — first base after intron).
    pub intron_end: u64,
}

/// Reference annotation: sets of known intron start and end positions per chromosome.
#[derive(Debug)]
pub struct ReferenceJunctions {
    /// Known intron start positions per chromosome (uppercased).
    pub intron_starts: HashMap<String, HashSet<u64>>,
    /// Known intron end positions per chromosome (uppercased).
    pub intron_ends: HashMap<String, HashSet<u64>>,
}

/// Aggregated results from junction annotation.
#[derive(Debug)]
pub struct JunctionResults {
    /// Per-junction read counts and classification.
    pub junctions: HashMap<Junction, (u64, JunctionClass)>,
    /// Total splicing events (per-read N-operations, includes short introns).
    pub total_events: u64,
    /// Events classified as known.
    pub known_events: u64,
    /// Events classified as partial novel.
    pub partial_novel_events: u64,
    /// Events classified as complete novel.
    pub complete_novel_events: u64,
    /// Events filtered by min_intron threshold.
    pub filtered_events: u64,
}

// ===================================================================
// BED12 reference parsing
// ===================================================================

/// Parse a BED12 file and extract reference intron start/end positions.
///
/// Single-exon transcripts are skipped. Chromosome names are uppercased.
///
/// # Arguments
/// * `bed_path` - Path to BED12 gene model file.
///
/// # Returns
/// A `ReferenceJunctions` with intron start and end position sets per chromosome.
pub fn parse_reference_bed(bed_path: &str) -> Result<ReferenceJunctions> {
    let content = std::fs::read_to_string(bed_path)
        .with_context(|| format!("Failed to read BED file: {}", bed_path))?;

    let mut intron_starts: HashMap<String, HashSet<u64>> = HashMap::new();
    let mut intron_ends: HashMap<String, HashSet<u64>> = HashMap::new();
    let mut transcript_count = 0u64;
    let mut skipped = 0u64;

    for line in content.lines() {
        let line = line.trim();
        if line.is_empty()
            || line.starts_with('#')
            || line.starts_with("track")
            || line.starts_with("browser")
        {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            debug!("Skipping BED line with {} fields (need 12)", fields.len());
            skipped += 1;
            continue;
        }

        let block_count: usize = fields[9].parse().unwrap_or(0);
        // Skip single-exon transcripts (matching RSeQC behavior)
        if block_count <= 1 {
            continue;
        }

        let chrom = fields[0].to_uppercase();
        let tx_start: u64 = fields[1].parse().unwrap_or(0);

        let block_sizes: Vec<u64> = fields[10]
            .trim_end_matches(',')
            .split(',')
            .filter_map(|s| s.parse().ok())
            .collect();
        let block_starts: Vec<u64> = fields[11]
            .trim_end_matches(',')
            .split(',')
            .filter_map(|s| s.parse().ok())
            .collect();

        if block_sizes.len() != block_count || block_starts.len() != block_count {
            debug!("Skipping BED line: block count mismatch");
            skipped += 1;
            continue;
        }

        // Compute exon coordinates
        let exon_starts: Vec<u64> = block_starts.iter().map(|&bs| tx_start + bs).collect();
        let exon_ends: Vec<u64> = exon_starts
            .iter()
            .zip(block_sizes.iter())
            .map(|(&s, &sz)| s + sz)
            .collect();

        // Introns are the gaps between consecutive exons
        let starts_set = intron_starts.entry(chrom.clone()).or_default();
        let ends_set = intron_ends.entry(chrom).or_default();

        for i in 0..block_count - 1 {
            let intron_start = exon_ends[i]; // first intronic base (0-based)
            let intron_end = exon_starts[i + 1]; // first base after intron (0-based, exclusive)
            starts_set.insert(intron_start);
            ends_set.insert(intron_end);
        }

        transcript_count += 1;
    }

    info!(
        "Parsed {} multi-exon transcripts from BED12 ({} lines skipped)",
        transcript_count, skipped
    );

    Ok(ReferenceJunctions {
        intron_starts,
        intron_ends,
    })
}

// ===================================================================
// CIGAR intron extraction
// ===================================================================

/// Extract intron intervals from a CIGAR string.
///
/// Matches RSeQC's `fetch_intron()` behavior: only CIGAR `N` operations produce
/// intron intervals. Soft clips do NOT advance position (unlike `fetch_exon()`).
/// `=` and `X` operations are ignored (do not advance position) — matching the
/// original Python bug.
///
/// # Arguments
/// * `start_pos` - Alignment start position (0-based, from BAM record).
/// * `cigar` - CIGAR operations from the BAM record.
///
/// # Returns
/// Vector of `(intron_start, intron_end)` tuples (0-based coordinates).
fn fetch_introns(start_pos: u64, cigar: &[rust_htslib::bam::record::Cigar]) -> Vec<(u64, u64)> {
    use rust_htslib::bam::record::Cigar::*;

    let mut pos = start_pos;
    let mut introns = Vec::new();

    for op in cigar {
        match op {
            Match(len) => pos += *len as u64, // M: advance position
            Ins(_) => {}                      // I: no position change
            Del(len) => pos += *len as u64,   // D: advance position
            RefSkip(len) => {
                // N: intron!
                let intron_start = pos;
                let intron_end = pos + *len as u64;
                introns.push((intron_start, intron_end));
                pos = intron_end;
            }
            SoftClip(_) => {} // S: no position change (unlike fetch_exon)
            HardClip(_) | Pad(_) | Equal(_) | Diff(_) => {} // ignored entirely
        }
    }

    introns
}

// ===================================================================
// Junction classification
// ===================================================================

/// Classify a junction as annotated, partial novel, or complete novel.
fn classify_junction(
    chrom: &str,
    intron_start: u64,
    intron_end: u64,
    reference: &ReferenceJunctions,
) -> JunctionClass {
    let start_known = reference
        .intron_starts
        .get(chrom)
        .is_some_and(|s| s.contains(&intron_start));
    let end_known = reference
        .intron_ends
        .get(chrom)
        .is_some_and(|s| s.contains(&intron_end));

    match (start_known, end_known) {
        (true, true) => JunctionClass::Annotated,
        (false, false) => JunctionClass::CompleteNovel,
        _ => JunctionClass::PartialNovel,
    }
}

// ===================================================================
// Main analysis
// ===================================================================

/// Run junction annotation analysis on a BAM file.
///
/// # Arguments
/// * `bam_path` - Path to input BAM file.
/// * `bed_path` - Path to BED12 reference gene model.
/// * `min_intron` - Minimum intron size to include (default 50).
/// * `mapq_cut` - Minimum mapping quality (default 30).
/// * `reference` - Optional path to reference FASTA (required for CRAM).
///
/// # Returns
/// `JunctionResults` containing all junction classifications and summary counts.
pub fn junction_annotation(
    bam_path: &str,
    bed_path: &str,
    min_intron: u64,
    mapq_cut: u8,
    reference: Option<&str>,
) -> Result<JunctionResults> {
    info!("Parsing reference gene model: {}", bed_path);
    let ref_junctions = parse_reference_bed(bed_path)?;

    info!("Processing BAM file: {}", bam_path);
    let mut bam = if let Some(ref_path) = reference {
        let mut r = bam::Reader::from_path(bam_path)
            .with_context(|| format!("Failed to open BAM file: {}", bam_path))?;
        r.set_reference(ref_path)
            .with_context(|| format!("Failed to set reference: {}", ref_path))?;
        r
    } else {
        bam::Reader::from_path(bam_path)
            .with_context(|| format!("Failed to open BAM file: {}", bam_path))?
    };
    // Build tid -> uppercase chrom name mapping
    let tid_to_chrom: HashMap<i32, String> = {
        let header_view = bam.header();
        let mut map = HashMap::new();
        for tid in 0..header_view.target_count() {
            if let Ok(name) = std::str::from_utf8(header_view.tid2name(tid)) {
                map.insert(tid as i32, name.to_uppercase());
            }
        }
        map
    };

    // Counters
    let mut total_events: u64 = 0;
    let mut known_events: u64 = 0;
    let mut partial_novel_events: u64 = 0;
    let mut complete_novel_events: u64 = 0;
    let mut filtered_events: u64 = 0;

    // Splicing events per unique junction
    let mut junction_counts: HashMap<Junction, u64> = HashMap::new();

    let mut total_reads: u64 = 0;
    let mut skipped_qcfail: u64 = 0;
    let mut skipped_dup: u64 = 0;
    let mut skipped_secondary: u64 = 0;
    let mut skipped_unmapped: u64 = 0;
    let mut skipped_mapq: u64 = 0;

    for result in bam.records() {
        let record = result.context("Failed to read BAM record")?;
        total_reads += 1;

        // Flag filtering — same order as RSeQC
        let flags = record.flags();
        if flags & 0x200 != 0 {
            skipped_qcfail += 1;
            continue;
        }
        if flags & 0x400 != 0 {
            skipped_dup += 1;
            continue;
        }
        if flags & 0x100 != 0 {
            skipped_secondary += 1;
            continue;
        }
        if record.is_unmapped() {
            skipped_unmapped += 1;
            continue;
        }
        if record.mapq() < mapq_cut {
            skipped_mapq += 1;
            continue;
        }

        let tid = record.tid();
        if tid < 0 {
            continue;
        }

        let chrom = match tid_to_chrom.get(&tid) {
            Some(c) => c,
            None => continue,
        };

        let start_pos = record.pos() as u64;
        let cigar = record.cigar();
        let introns = fetch_introns(start_pos, cigar.as_ref());

        for (intron_start, intron_end) in introns {
            total_events += 1;

            // Filter by minimum intron size (after counting toward total_events)
            let intron_size = intron_end.saturating_sub(intron_start);
            if intron_size < min_intron {
                filtered_events += 1;
                continue;
            }

            let class = classify_junction(chrom, intron_start, intron_end, &ref_junctions);

            match class {
                JunctionClass::Annotated => known_events += 1,
                JunctionClass::PartialNovel => partial_novel_events += 1,
                JunctionClass::CompleteNovel => complete_novel_events += 1,
            }

            let junction = Junction {
                chrom: chrom.clone(),
                intron_start,
                intron_end,
            };
            *junction_counts.entry(junction).or_insert(0) += 1;
        }
    }

    debug!(
        "Processed {} reads: {} qcfail, {} dup, {} secondary, {} unmapped, {} low mapq",
        total_reads, skipped_qcfail, skipped_dup, skipped_secondary, skipped_unmapped, skipped_mapq
    );

    // Classify each unique junction
    let mut junctions: HashMap<Junction, (u64, JunctionClass)> = HashMap::new();
    for (junction, count) in &junction_counts {
        let class = classify_junction(
            &junction.chrom,
            junction.intron_start,
            junction.intron_end,
            &ref_junctions,
        );
        junctions.insert(junction.clone(), (*count, class));
    }

    Ok(JunctionResults {
        junctions,
        total_events,
        known_events,
        partial_novel_events,
        complete_novel_events,
        filtered_events,
    })
}

// ===================================================================
// Output formatting
// ===================================================================

/// Convert an uppercased chromosome name back to lowercase for output.
///
/// Matches RSeQC's `chrom.replace("CHR","chr")` behavior.
fn format_chrom(chrom: &str) -> String {
    chrom.replace("CHR", "chr")
}

/// Write the junction annotation XLS (TSV) file.
///
/// Format matches RSeQC's `.junction.xls` output.
pub fn write_junction_xls(results: &JunctionResults, path: &Path) -> Result<()> {
    let mut f = std::fs::File::create(path)
        .with_context(|| format!("Failed to create file: {}", path.display()))?;

    writeln!(
        f,
        "chrom\tintron_st(0-based)\tintron_end(1-based)\tread_count\tannotation"
    )?;

    // Sort junctions for deterministic output (RSeQC output order is non-deterministic)
    let mut sorted_junctions: Vec<_> = results.junctions.iter().collect();
    sorted_junctions.sort_by(|(a, _), (b, _)| {
        a.chrom
            .cmp(&b.chrom)
            .then(a.intron_start.cmp(&b.intron_start))
            .then(a.intron_end.cmp(&b.intron_end))
    });

    for (junction, (count, class)) in &sorted_junctions {
        // RSeQC output has a tab then space before the annotation label
        writeln!(
            f,
            "{}\t{}\t{}\t{}\t {}",
            format_chrom(&junction.chrom),
            junction.intron_start,
            junction.intron_end,
            count,
            class.label()
        )?;
    }

    Ok(())
}

/// Write the junction BED12 file.
///
/// Each junction is represented with 1bp exon blocks flanking the intron,
/// matching the reference output format.
pub fn write_junction_bed(results: &JunctionResults, path: &Path) -> Result<()> {
    let mut f = std::fs::File::create(path)
        .with_context(|| format!("Failed to create file: {}", path.display()))?;

    // Sort for deterministic output
    let mut sorted_junctions: Vec<_> = results.junctions.iter().collect();
    sorted_junctions.sort_by(|(a, _), (b, _)| {
        a.chrom
            .cmp(&b.chrom)
            .then(a.intron_start.cmp(&b.intron_start))
            .then(a.intron_end.cmp(&b.intron_end))
    });

    for (junction, (count, class)) in &sorted_junctions {
        let chrom = format_chrom(&junction.chrom);
        // BED coordinates: 1bp block before intron and 1bp block after
        let bed_start = junction.intron_start.saturating_sub(1); // 1bp before intron
        let bed_end = junction.intron_end + 1; // 1bp after intron
        let span = bed_end - bed_start;

        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}\t.\t{}\t{}\t{}\t2\t1,1\t0,{}",
            chrom,
            bed_start,
            bed_end,
            class.label(),
            count,
            bed_start,
            bed_end,
            class.color(),
            span - 1
        )?;
    }

    Ok(())
}

/// Write the R script for junction pie charts.
///
/// Generates two pie charts: splice events and splice junctions.
pub fn write_junction_plot_r(results: &JunctionResults, prefix: &str, path: &Path) -> Result<()> {
    let mut f = std::fs::File::create(path)
        .with_context(|| format!("Failed to create file: {}", path.display()))?;

    // Event-level percentages
    let total_classified_events =
        results.known_events + results.partial_novel_events + results.complete_novel_events;
    let (e_known_pct, e_partial_pct, e_novel_pct) = if total_classified_events > 0 {
        (
            results.known_events as f64 * 100.0 / total_classified_events as f64,
            results.partial_novel_events as f64 * 100.0 / total_classified_events as f64,
            results.complete_novel_events as f64 * 100.0 / total_classified_events as f64,
        )
    } else {
        (0.0, 0.0, 0.0)
    };

    // Junction-level counts
    let mut junc_known: u64 = 0;
    let mut junc_partial: u64 = 0;
    let mut junc_novel: u64 = 0;
    for (_, class) in results.junctions.values() {
        match class {
            JunctionClass::Annotated => junc_known += 1,
            JunctionClass::PartialNovel => junc_partial += 1,
            JunctionClass::CompleteNovel => junc_novel += 1,
        }
    }
    let total_junctions = junc_known + junc_partial + junc_novel;
    let (j_known_pct, j_partial_pct, j_novel_pct) = if total_junctions > 0 {
        (
            junc_known as f64 * 100.0 / total_junctions as f64,
            junc_partial as f64 * 100.0 / total_junctions as f64,
            junc_novel as f64 * 100.0 / total_junctions as f64,
        )
    } else {
        (0.0, 0.0, 0.0)
    };

    // Events pie chart
    writeln!(f, "pdf(\"{}.splice_events.pdf\")", prefix)?;
    writeln!(
        f,
        "events=c({:.3},{:.3},{:.3})",
        e_partial_pct, e_novel_pct, e_known_pct
    )?;
    writeln!(
        f,
        "pie(events,col=c(2,3,4),init.angle=30,angle=c(60,120,150),density=c(70,70,70),main=\"splicing events\",labels=c(\"partial_novel {}%\",\"complete_novel {}%\",\"known {}%\"))",
        e_partial_pct.round() as u64,
        e_novel_pct.round() as u64,
        e_known_pct.round() as u64
    )?;
    writeln!(f, "dev.off()")?;

    // Junctions pie chart
    writeln!(f, "pdf(\"{}.splice_junction.pdf\")", prefix)?;
    writeln!(
        f,
        "junction=c({:.3},{:.3},{:.3})",
        j_partial_pct, j_novel_pct, j_known_pct
    )?;
    writeln!(
        f,
        "pie(junction,col=c(2,3,4),init.angle=30,angle=c(60,120,150),density=c(70,70,70),main=\"splicing junctions\",labels=c(\"partial_novel {}%\",\"complete_novel {}%\",\"known {}%\"))",
        j_partial_pct.round() as u64,
        j_novel_pct.round() as u64,
        j_known_pct.round() as u64
    )?;
    writeln!(f, "dev.off()")?;

    Ok(())
}

/// Write the summary statistics to a file.
///
/// Format matches RSeQC's text summary output.
pub fn write_summary(results: &JunctionResults, path: &Path) -> Result<()> {
    let mut f = std::fs::File::create(path)
        .with_context(|| format!("Failed to create file: {}", path.display()))?;

    // Junction-level counts
    let mut junc_known: u64 = 0;
    let mut junc_partial: u64 = 0;
    let mut junc_novel: u64 = 0;
    for (_, class) in results.junctions.values() {
        match class {
            JunctionClass::Annotated => junc_known += 1,
            JunctionClass::PartialNovel => junc_partial += 1,
            JunctionClass::CompleteNovel => junc_novel += 1,
        }
    }
    let total_junctions = junc_known + junc_partial + junc_novel;

    writeln!(
        f,
        "==================================================================="
    )?;

    // Events section
    writeln!(f, "Total splicing  Events:\t{}", results.total_events)?;
    writeln!(f, "Known Splicing Events:\t{}", results.known_events)?;
    writeln!(
        f,
        "Partial Novel Splicing Events:\t{}",
        results.partial_novel_events
    )?;
    writeln!(
        f,
        "Novel Splicing Events:\t{}",
        results.complete_novel_events
    )?;
    writeln!(f, "Filtered Splicing Events:\t{}", results.filtered_events)?;

    // Junctions section
    writeln!(f)?;
    writeln!(f, "Total splicing  Junctions:\t{}", total_junctions)?;
    writeln!(f, "Known Splicing Junctions:\t{}", junc_known)?;
    writeln!(f, "Partial Novel Splicing Junctions:\t{}", junc_partial)?;
    writeln!(f, "Novel Splicing Junctions:\t{}", junc_novel)?;

    writeln!(f)?;
    writeln!(
        f,
        "==================================================================="
    )?;

    Ok(())
}

/// Print the summary to stderr (matching RSeQC behavior).
pub fn print_summary(results: &JunctionResults) {
    // Junction-level counts
    let mut junc_known: u64 = 0;
    let mut junc_partial: u64 = 0;
    let mut junc_novel: u64 = 0;
    for (_, class) in results.junctions.values() {
        match class {
            JunctionClass::Annotated => junc_known += 1,
            JunctionClass::PartialNovel => junc_partial += 1,
            JunctionClass::CompleteNovel => junc_novel += 1,
        }
    }
    let total_junctions = junc_known + junc_partial + junc_novel;

    eprintln!("===================================================================");
    eprintln!("Total splicing  Events:\t{}", results.total_events);
    eprintln!("Known Splicing Events:\t{}", results.known_events);
    eprintln!(
        "Partial Novel Splicing Events:\t{}",
        results.partial_novel_events
    );
    eprintln!("Novel Splicing Events:\t{}", results.complete_novel_events);
    eprintln!("Filtered Splicing Events:\t{}", results.filtered_events);
    eprintln!();
    eprintln!("Total splicing  Junctions:\t{}", total_junctions);
    eprintln!("Known Splicing Junctions:\t{}", junc_known);
    eprintln!("Partial Novel Splicing Junctions:\t{}", junc_partial);
    eprintln!("Novel Splicing Junctions:\t{}", junc_novel);
    eprintln!();
    eprintln!("===================================================================");
}

// ===================================================================
// Unit tests
// ===================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fetch_introns_simple() {
        use rust_htslib::bam::record::Cigar::*;
        // 50M500N50M — one intron at position 100+50=150 to 150+500=650
        let cigar = vec![Match(50), RefSkip(500), Match(50)];
        let introns = fetch_introns(100, &cigar);
        assert_eq!(introns.len(), 1);
        assert_eq!(introns[0], (150, 650));
    }

    #[test]
    fn test_fetch_introns_multiple() {
        use rust_htslib::bam::record::Cigar::*;
        // 10M500N20M300N10M — two introns
        let cigar = vec![Match(10), RefSkip(500), Match(20), RefSkip(300), Match(10)];
        let introns = fetch_introns(100, &cigar);
        assert_eq!(introns.len(), 2);
        assert_eq!(introns[0], (110, 610)); // first intron
        assert_eq!(introns[1], (630, 930)); // second intron
    }

    #[test]
    fn test_fetch_introns_with_deletions() {
        use rust_htslib::bam::record::Cigar::*;
        // 10M5D10M500N10M
        let cigar = vec![Match(10), Del(5), Match(10), RefSkip(500), Match(10)];
        let introns = fetch_introns(100, &cigar);
        assert_eq!(introns.len(), 1);
        assert_eq!(introns[0], (125, 625)); // 100+10+5+10=125
    }

    #[test]
    fn test_fetch_introns_no_introns() {
        use rust_htslib::bam::record::Cigar::*;
        let cigar = vec![Match(100)];
        let introns = fetch_introns(100, &cigar);
        assert!(introns.is_empty());
    }

    #[test]
    fn test_fetch_introns_soft_clip_no_advance() {
        use rust_htslib::bam::record::Cigar::*;
        // 5S50M500N50M — soft clip should NOT advance position
        let cigar = vec![SoftClip(5), Match(50), RefSkip(500), Match(50)];
        let introns = fetch_introns(100, &cigar);
        assert_eq!(introns.len(), 1);
        assert_eq!(introns[0], (150, 650)); // same as without soft clip
    }

    #[test]
    fn test_classify_junction_annotated() {
        let mut starts = HashMap::new();
        let mut ends = HashMap::new();
        starts.insert("CHR1".to_string(), HashSet::from([100, 500]));
        ends.insert("CHR1".to_string(), HashSet::from([400, 900]));
        let reference = ReferenceJunctions {
            intron_starts: starts,
            intron_ends: ends,
        };

        assert_eq!(
            classify_junction("CHR1", 100, 400, &reference),
            JunctionClass::Annotated
        );
    }

    #[test]
    fn test_classify_junction_partial_novel() {
        let mut starts = HashMap::new();
        let mut ends = HashMap::new();
        starts.insert("CHR1".to_string(), HashSet::from([100]));
        ends.insert("CHR1".to_string(), HashSet::from([400]));
        let reference = ReferenceJunctions {
            intron_starts: starts,
            intron_ends: ends,
        };

        // Start known, end unknown
        assert_eq!(
            classify_junction("CHR1", 100, 999, &reference),
            JunctionClass::PartialNovel
        );
        // Start unknown, end known
        assert_eq!(
            classify_junction("CHR1", 999, 400, &reference),
            JunctionClass::PartialNovel
        );
    }

    #[test]
    fn test_classify_junction_complete_novel() {
        let mut starts = HashMap::new();
        let mut ends = HashMap::new();
        starts.insert("CHR1".to_string(), HashSet::from([100]));
        ends.insert("CHR1".to_string(), HashSet::from([400]));
        let reference = ReferenceJunctions {
            intron_starts: starts,
            intron_ends: ends,
        };

        assert_eq!(
            classify_junction("CHR1", 999, 888, &reference),
            JunctionClass::CompleteNovel
        );
    }

    #[test]
    fn test_classify_junction_unknown_chrom() {
        let reference = ReferenceJunctions {
            intron_starts: HashMap::new(),
            intron_ends: HashMap::new(),
        };

        assert_eq!(
            classify_junction("CHRX", 100, 400, &reference),
            JunctionClass::CompleteNovel
        );
    }

    #[test]
    fn test_junction_class_labels() {
        assert_eq!(JunctionClass::Annotated.label(), "annotated");
        assert_eq!(JunctionClass::PartialNovel.label(), "partial_novel");
        assert_eq!(JunctionClass::CompleteNovel.label(), "complete_novel");
    }

    #[test]
    fn test_junction_class_colors() {
        assert_eq!(JunctionClass::Annotated.color(), "205,0,0");
        assert_eq!(JunctionClass::PartialNovel.color(), "0,205,0");
        assert_eq!(JunctionClass::CompleteNovel.color(), "0,0,205");
    }

    #[test]
    fn test_format_chrom() {
        assert_eq!(format_chrom("CHR1"), "chr1");
        assert_eq!(format_chrom("CHRX"), "chrX");
        assert_eq!(format_chrom("6"), "6"); // numeric chroms unchanged
    }
}
