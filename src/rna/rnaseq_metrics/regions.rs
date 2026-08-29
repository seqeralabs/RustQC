//! Genomic region sets built from a gene annotation.
//!
//! # Upstream semantics
//!
//! Follows Picard `CollectRnaSeqMetrics`. Every aligned base of a qualifying
//! read is assigned to exactly one class, and the classes are checked in a
//! fixed order because annotations overlap: a base inside a coding sequence of
//! one transcript and an intron of another counts as coding.
//!
//! "UTR" means exonic but not coding, which is what Picard reports under that
//! name. A transcript with no annotated CDS therefore contributes only UTR,
//! never coding, and on a non-coding annotation `CODING_BASES` is legitimately
//! zero rather than missing.

use std::collections::HashMap;

use crate::gtf::Gene;

/// What a base was assigned to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Region {
    /// Inside an annotated coding sequence.
    Coding,
    /// Exonic but outside any coding sequence.
    Utr,
    /// Within a transcript's span but not exonic.
    Intronic,
    /// Outside every transcript.
    Intergenic,
}

/// Merged, sorted interval lists per contig.
#[derive(Debug, Default)]
pub struct RegionSets {
    coding: HashMap<String, Vec<(u64, u64)>>,
    exonic: HashMap<String, Vec<(u64, u64)>>,
    transcribed: HashMap<String, Vec<(u64, u64)>>,
}

/// Sort and merge overlapping half-open intervals.
fn merge(mut intervals: Vec<(u64, u64)>) -> Vec<(u64, u64)> {
    if intervals.is_empty() {
        return intervals;
    }
    intervals.sort_unstable();
    let mut merged: Vec<(u64, u64)> = Vec::with_capacity(intervals.len());
    for (start, end) in intervals {
        match merged.last_mut() {
            Some(last) if start <= last.1 => last.1 = last.1.max(end),
            _ => merged.push((start, end)),
        }
    }
    merged
}

/// Whether a sorted, merged interval list contains a position.
fn contains(intervals: &[(u64, u64)], position: u64) -> bool {
    intervals
        .binary_search_by(|(start, end)| {
            if position < *start {
                std::cmp::Ordering::Greater
            } else if position >= *end {
                std::cmp::Ordering::Less
            } else {
                std::cmp::Ordering::Equal
            }
        })
        .is_ok()
}

impl RegionSets {
    /// Build the region sets from parsed genes.
    ///
    /// The GTF parser keeps coordinates one-based and inclusive, as the file
    /// has them, while alignment positions are zero-based. Everything is
    /// converted here to zero-based half-open, once, so the classifier
    /// compares like with like. Skipping the conversion shifts every boundary
    /// by one base, which leaves interior bases right and quietly misassigns
    /// the edges: on the project fixture exactly one base moved from UTR to
    /// intergenic.
    pub fn from_genes<'a, I>(genes: I) -> Self
    where
        I: IntoIterator<Item = &'a Gene>,
    {
        let mut coding: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
        let mut exonic: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
        let mut transcribed: HashMap<String, Vec<(u64, u64)>> = HashMap::new();

        for gene in genes {
            for transcript in &gene.transcripts {
                let chrom = transcript.chrom.clone();
                transcribed
                    .entry(chrom.clone())
                    .or_default()
                    .push((transcript.start.saturating_sub(1), transcript.end));
                for (start, end) in &transcript.exons {
                    exonic
                        .entry(chrom.clone())
                        .or_default()
                        .push((start.saturating_sub(1), *end));
                }
                // A transcript without a CDS contributes no coding bases; its
                // exons are all UTR, which is what Picard reports.
                if let (Some(cds_start), Some(cds_end)) = (transcript.cds_start, transcript.cds_end)
                {
                    if cds_end > cds_start {
                        let cds_start = cds_start.saturating_sub(1);
                        // Only the exonic part of the CDS span is coding; the
                        // span itself can cross introns.
                        for (start, end) in &transcript.exons {
                            let overlap_start = start.saturating_sub(1).max(cds_start);
                            let overlap_end = (*end).min(cds_end);
                            if overlap_end > overlap_start {
                                coding
                                    .entry(chrom.clone())
                                    .or_default()
                                    .push((overlap_start, overlap_end));
                            }
                        }
                    }
                }
            }
            // A gene whose exons were parsed without transcripts still
            // contributes; otherwise a GTF with no transcript rows would look
            // entirely intergenic.
            if gene.transcripts.is_empty() {
                transcribed
                    .entry(gene.chrom.clone())
                    .or_default()
                    .push((gene.start.saturating_sub(1), gene.end));
                for exon in &gene.exons {
                    exonic
                        .entry(exon.chrom.clone())
                        .or_default()
                        .push((exon.start.saturating_sub(1), exon.end));
                }
            }
        }

        Self {
            coding: coding.into_iter().map(|(k, v)| (k, merge(v))).collect(),
            exonic: exonic.into_iter().map(|(k, v)| (k, merge(v))).collect(),
            transcribed: transcribed
                .into_iter()
                .map(|(k, v)| (k, merge(v)))
                .collect(),
        }
    }

    /// Classify one zero-based position.
    ///
    /// The order is fixed and matters: annotations overlap, and a base that is
    /// coding in one transcript and intronic in another is coding.
    pub fn classify(&self, chrom: &str, position: u64) -> Region {
        if self
            .coding
            .get(chrom)
            .is_some_and(|set| contains(set, position))
        {
            return Region::Coding;
        }
        if self
            .exonic
            .get(chrom)
            .is_some_and(|set| contains(set, position))
        {
            return Region::Utr;
        }
        if self
            .transcribed
            .get(chrom)
            .is_some_and(|set| contains(set, position))
        {
            return Region::Intronic;
        }
        Region::Intergenic
    }

    /// Total bases of each class in the annotation itself, which is the
    /// territory the metrics are reported against.
    pub fn territories(&self) -> (u64, u64, u64) {
        let span = |sets: &HashMap<String, Vec<(u64, u64)>>| -> u64 {
            sets.values()
                .flat_map(|v| v.iter())
                .map(|(start, end)| end - start)
                .sum()
        };
        let coding = span(&self.coding);
        let exonic = span(&self.exonic);
        let transcribed = span(&self.transcribed);
        (
            coding,
            exonic.saturating_sub(coding),
            transcribed.saturating_sub(exonic),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn intervals_merge_when_they_touch_or_overlap() {
        assert_eq!(merge(vec![(0, 10), (5, 20)]), vec![(0, 20)]);
        assert_eq!(merge(vec![(0, 10), (10, 20)]), vec![(0, 20)], "touching");
        assert_eq!(merge(vec![(0, 10), (11, 20)]), vec![(0, 10), (11, 20)]);
        assert_eq!(
            merge(vec![(5, 20), (0, 10)]),
            vec![(0, 20)],
            "unsorted input"
        );
    }

    #[test]
    fn containment_is_half_open() {
        let set = vec![(10, 20)];
        assert!(!contains(&set, 9));
        assert!(contains(&set, 10), "the start is inside");
        assert!(contains(&set, 19));
        assert!(!contains(&set, 20), "the end is not");
    }

    #[test]
    fn classification_prefers_coding_then_utr_then_intronic() {
        let mut sets = RegionSets::default();
        sets.coding.insert("chr1".into(), vec![(100, 200)]);
        sets.exonic.insert("chr1".into(), vec![(50, 300)]);
        sets.transcribed.insert("chr1".into(), vec![(0, 1000)]);

        assert_eq!(sets.classify("chr1", 150), Region::Coding);
        assert_eq!(sets.classify("chr1", 60), Region::Utr, "exonic, not coding");
        assert_eq!(sets.classify("chr1", 500), Region::Intronic);
        assert_eq!(sets.classify("chr1", 2000), Region::Intergenic);
        assert_eq!(
            sets.classify("chr2", 150),
            Region::Intergenic,
            "an unannotated contig is entirely intergenic"
        );
    }

    #[test]
    fn territories_subtract_the_nested_classes() {
        let mut sets = RegionSets::default();
        sets.coding.insert("chr1".into(), vec![(100, 200)]);
        sets.exonic.insert("chr1".into(), vec![(50, 300)]);
        sets.transcribed.insert("chr1".into(), vec![(0, 1000)]);
        let (coding, utr, intronic) = sets.territories();
        assert_eq!(coding, 100);
        assert_eq!(utr, 150, "exonic 250 less the 100 coding");
        assert_eq!(intronic, 750, "transcribed 1000 less the 250 exonic");
    }
}
