//! Defects a protein FASTA can carry that seqkit does not report.
//!
//! These are the things that make a proteome unusable downstream: a residue
//! that is not an amino acid, a stop codon in the middle of a sequence, a
//! sequence that is a byte-for-byte duplicate of another, or an identifier
//! used twice for different sequences.

use std::collections::HashMap;

use super::Record;

/// The twenty standard amino acids, plus the two that appear in translated
/// proteomes: `U` for selenocysteine and `O` for pyrrolysine.
pub const STANDARD_RESIDUES: &[u8] = b"ACDEFGHIKLMNPQRSTVWYUO";

/// Ambiguity codes that are valid in a protein FASTA but carry no identity:
/// `B` is D or N, `Z` is E or Q, `J` is I or L, and `X` is anything.
pub const AMBIGUOUS_RESIDUES: &[u8] = b"BZJX";

/// What is wrong with a proteome, if anything.
#[derive(Debug, Clone, Default)]
pub struct Defects {
    /// Sequences holding a `*` anywhere but the final position.
    pub internal_stops: Vec<String>,
    /// Sequences holding a residue that is neither standard, ambiguous, a gap
    /// nor a terminal stop, with the offending residue.
    pub non_standard: Vec<(String, u8)>,
    /// Sequences whose residues are identical to an earlier one, paired with
    /// the identifier they duplicate.
    pub duplicates: Vec<(String, String)>,
    /// Identifiers used more than once.
    pub duplicate_ids: Vec<String>,
    /// Sequences with no residues at all.
    pub empty: Vec<String>,
    /// Sequences not ending in `*`, reported only when asked for.
    pub missing_terminal_stop: Vec<String>,
    /// Residues that are ambiguity codes.
    pub ambiguous_residues: u64,
}

impl Defects {
    /// Whether anything at all was found.
    pub fn is_clean(&self) -> bool {
        self.internal_stops.is_empty()
            && self.non_standard.is_empty()
            && self.duplicates.is_empty()
            && self.duplicate_ids.is_empty()
            && self.empty.is_empty()
            && self.missing_terminal_stop.is_empty()
    }

    /// Total number of defective sequences, counting a sequence once per
    /// category it falls into.
    pub fn count(&self) -> usize {
        self.internal_stops.len()
            + self.non_standard.len()
            + self.duplicates.len()
            + self.duplicate_ids.len()
            + self.empty.len()
            + self.missing_terminal_stop.len()
    }
}

/// Inspect every record.
///
/// `expect_stop` turns a missing terminal `*` into a reported defect. It is off
/// by default because most reference proteomes do not carry one.
pub fn inspect(records: &[Record], expect_stop: bool) -> Defects {
    let mut defects = Defects::default();
    let mut seen_sequences: HashMap<&[u8], &str> = HashMap::new();
    let mut seen_ids: HashMap<&str, usize> = HashMap::new();

    for record in records {
        *seen_ids.entry(record.id.as_str()).or_insert(0) += 1;

        if record.is_empty() {
            defects.empty.push(record.id.clone());
            continue;
        }

        // A trailing stop is expected; anywhere else it truncates the protein.
        let body = record
            .residues
            .strip_suffix(b"*")
            .unwrap_or(&record.residues);
        if body.contains(&b'*') {
            defects.internal_stops.push(record.id.clone());
        }
        if expect_stop && !record.residues.ends_with(b"*") {
            defects.missing_terminal_stop.push(record.id.clone());
        }

        for residue in body {
            if AMBIGUOUS_RESIDUES.contains(residue) {
                defects.ambiguous_residues += 1;
            } else if !STANDARD_RESIDUES.contains(residue)
                && *residue != b'-'
                && *residue != b'.'
                && *residue != b'*'
            {
                defects.non_standard.push((record.id.clone(), *residue));
                break;
            }
        }

        match seen_sequences.get(record.residues.as_slice()) {
            Some(first) => defects
                .duplicates
                .push((record.id.clone(), (*first).to_string())),
            None => {
                seen_sequences.insert(record.residues.as_slice(), record.id.as_str());
            }
        }
    }

    for (id, count) in seen_ids {
        if count > 1 {
            defects.duplicate_ids.push(id.to_string());
        }
    }
    defects.duplicate_ids.sort();

    defects
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(id: &str, residues: &str) -> Record {
        Record {
            id: id.to_string(),
            residues: residues.as_bytes().to_vec(),
        }
    }

    #[test]
    fn a_terminal_stop_is_fine_and_an_internal_one_is_not() {
        let defects = inspect(&[record("ok", "MKTAYI*"), record("bad", "MKT*AYI")], false);
        assert_eq!(defects.internal_stops, vec!["bad".to_string()]);
    }

    #[test]
    fn ambiguity_codes_are_counted_rather_than_flagged() {
        let defects = inspect(&[record("a", "MKTXBZJ")], false);
        assert!(defects.non_standard.is_empty(), "X, B, Z and J are valid");
        assert_eq!(defects.ambiguous_residues, 4);
    }

    #[test]
    fn a_residue_that_is_not_an_amino_acid_is_flagged_once() {
        let defects = inspect(&[record("a", "MK1T2AYI")], false);
        assert_eq!(defects.non_standard.len(), 1, "one report per sequence");
        assert_eq!(defects.non_standard[0].1, b'1', "the first offender");
    }

    #[test]
    fn selenocysteine_and_pyrrolysine_are_standard_enough() {
        let defects = inspect(&[record("a", "MKUTAYIO")], false);
        assert!(defects.non_standard.is_empty(), "U and O are real residues");
    }

    #[test]
    fn identical_sequences_are_paired_with_the_first_that_carried_them() {
        let defects = inspect(
            &[
                record("first", "MKTAYI"),
                record("second", "MSEQ"),
                record("third", "MKTAYI"),
            ],
            false,
        );
        assert_eq!(
            defects.duplicates,
            vec![("third".to_string(), "first".to_string())]
        );
    }

    #[test]
    fn a_repeated_identifier_is_reported_even_with_different_sequences() {
        let defects = inspect(&[record("a", "MKT"), record("a", "MSEQ")], false);
        assert_eq!(defects.duplicate_ids, vec!["a".to_string()]);
        assert!(defects.duplicates.is_empty(), "the sequences differ");
    }

    #[test]
    fn a_missing_terminal_stop_is_reported_only_when_asked_for() {
        let records = [record("a", "MKTAYI")];
        assert!(inspect(&records, false).missing_terminal_stop.is_empty());
        assert_eq!(
            inspect(&records, true).missing_terminal_stop,
            vec!["a".to_string()]
        );
    }

    #[test]
    fn a_clean_proteome_reports_clean() {
        let defects = inspect(&[record("a", "MKTAYI"), record("b", "MSEQ")], false);
        assert!(defects.is_clean());
        assert_eq!(defects.count(), 0);
    }
}
