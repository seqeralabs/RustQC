//! Coverage distribution (`COV` section) for samtools stats output.
//!
//! `samtools stats` builds the coverage histogram with a round buffer that is
//! flushed as reads stream past in coordinate order (`round_buffer_flush` and
//! `round_buffer_insert_read` in samtools `stats.c`). Two quirks of that buffer
//! show up in the output, so they are reproduced here:
//!
//! - A flush to "the end" (`round_buffer_flush(stats, -1)`, called on every
//!   contig change and once at end of input) only visits `size - 1` slots. One
//!   slot keeps its depth: it is carried into the next contig, where it is
//!   counted on that contig's first visit of the slot (on its own, or added to
//!   the depth already there), and it is never counted after the last contig.
//! - Reference positions more than `size` bases past the read start (long
//!   introns) wrap around and add to the depth of an earlier slot.
//!
//! RustQC splits the BAM by contig across threads, so each contig is simulated
//! on its own from the state samtools starts it with (`start = 0`, `pos = -1`,
//! empty buffer). The slot carried in from the previous contig is resolved in
//! [`resolve_carries`] once all contigs are known: only its first visit in the
//! next contig depends on it. This keeps the result independent of `--threads`.

use std::collections::HashMap;

/// Initial round buffer size in samtools: `nbases * 5` with `nbases = 300`.
pub(crate) const INITIAL_BUFFER_SIZE: usize = 1500;

/// Sentinel for a slot that no flush has visited yet.
const UNVISITED: u32 = u32::MAX;

/// `round_buffer_lidx2ridx`: buffer slot holding reference position `pos`.
fn lidx2ridx(offset: usize, size: usize, refpos: i64, pos: i64) -> usize {
    let size = size as i64;
    (offset as i64 + (pos - refpos) % size).rem_euclid(size) as usize
}

/// Round buffer for one contig, mirroring samtools' `round_buffer_t`.
#[derive(Debug)]
pub(crate) struct ContigCoverage {
    tid: i32,
    buffer: Vec<u32>,
    start: usize,
    pos: i64,
    /// Depth of each slot when a flush first visited it (`UNVISITED` if not yet).
    /// Dropped once every slot has been visited.
    first_visit: Option<Vec<u32>>,
    /// Non-zero first-visit depths, `(slot, depth)`, filled from `first_visit`.
    first_nonzero: Vec<(u32, u32)>,
    unvisited: usize,
    /// The buffer was enlarged for a long read (slot numbers no longer line up
    /// with the contig's initial layout).
    grown: bool,
}

impl ContigCoverage {
    /// Start a contig in the state samtools leaves after `round_buffer_flush(-1)`.
    pub(crate) fn new(tid: i32) -> Self {
        Self {
            tid,
            buffer: vec![0; INITIAL_BUFFER_SIZE],
            start: 0,
            pos: -1,
            first_visit: Some(vec![UNVISITED; INITIAL_BUFFER_SIZE]),
            first_nonzero: Vec::new(),
            unvisited: INITIAL_BUFFER_SIZE,
            grown: false,
        }
    }

    pub(crate) fn tid(&self) -> i32 {
        self.tid
    }

    /// Make room for a read of `read_len` bases (samtools sizes the buffer at
    /// five read lengths). Keeps pending depths, linearised from `start`.
    pub(crate) fn ensure_read_len(&mut self, read_len: usize) {
        let need = read_len * 5;
        if need <= self.buffer.len() {
            return;
        }
        let mut new_buf = vec![0u32; need];
        let head = self.buffer.len() - self.start;
        new_buf[..head].copy_from_slice(&self.buffer[self.start..]);
        new_buf[head..head + self.start].copy_from_slice(&self.buffer[..self.start]);
        self.buffer = new_buf;
        self.start = 0;
        // Slot numbers changed; stop tracking first visits for the carry.
        self.take_first_visits();
        self.grown = true;
    }

    /// Record the first visit of `slot`, then count and clear its depth.
    #[inline]
    fn visit(&mut self, slot: usize, hist: &mut HashMap<u32, u64>) {
        let depth = self.buffer[slot];
        if let Some(first) = self.first_visit.as_mut() {
            if first[slot] == UNVISITED {
                first[slot] = depth;
                self.unvisited -= 1;
                if self.unvisited == 0 {
                    self.take_first_visits();
                }
            }
        }
        if depth != 0 {
            *hist.entry(depth).or_insert(0) += 1;
            self.buffer[slot] = 0;
        }
    }

    /// Move the non-zero first visits into `first_nonzero` and stop tracking.
    fn take_first_visits(&mut self) {
        if let Some(first) = self.first_visit.take() {
            self.first_nonzero.extend(
                first
                    .iter()
                    .enumerate()
                    .filter(|&(_, &d)| d != UNVISITED && d != 0)
                    .map(|(slot, &d)| (slot as u32, d)),
            );
        }
    }

    /// `round_buffer_flush`: count every position before `pos` (`-1` = flush
    /// the whole buffer, which visits only `size - 1` slots).
    pub(crate) fn flush(&mut self, pos: i64, hist: &mut HashMap<u32, u64>) {
        if pos == self.pos {
            return;
        }
        let size = self.buffer.len();
        let new_pos = pos;
        let mut pos = pos;
        if pos == -1 || pos - self.pos >= size as i64 {
            pos = self.pos + size as i64 - 1;
        }
        let mut ifrom = self.start;
        let ito = lidx2ridx(self.start, size, self.pos, pos - 1);
        if ifrom > ito {
            for slot in ifrom..size {
                self.visit(slot, hist);
            }
            ifrom = 0;
        }
        for slot in ifrom..=ito {
            self.visit(slot, hist);
        }
        self.start = if new_pos == -1 {
            0
        } else {
            lidx2ridx(self.start, size, self.pos, pos)
        };
        self.pos = new_pos;
    }

    /// `round_buffer_insert_read`: add one to the depth of `[from, to)`.
    ///
    /// Like samtools, a block exactly `size` long adds nothing and positions
    /// past the end of the buffer wrap around onto earlier slots.
    pub(crate) fn insert(&mut self, from: i64, to: i64) {
        let size = self.buffer.len();
        let mut ifrom = lidx2ridx(self.start, size, self.pos, from);
        let ito = lidx2ridx(self.start, size, self.pos, to);
        if ifrom > ito {
            for d in &mut self.buffer[ifrom..] {
                *d += 1;
            }
            ifrom = 0;
        }
        for d in &mut self.buffer[ifrom..ito] {
            *d += 1;
        }
    }

    /// End the contig the way samtools does on a contig change or at end of
    /// input: `round_buffer_flush(-1)`, which leaves one slot unflushed.
    pub(crate) fn finish(mut self, hist: &mut HashMap<u32, u64>) -> ContigCovSummary {
        let size = self.buffer.len();
        let end_slot = (self.start + size - 1) % size;
        let end_slot_visited = match &self.first_visit {
            Some(first) => first[end_slot] != UNVISITED,
            None => true,
        };
        self.flush(-1, hist);
        let end_depth = self.buffer[end_slot];
        self.take_first_visits();
        self.first_nonzero.sort_unstable();
        ContigCovSummary {
            tid: self.tid,
            first_nonzero: self.first_nonzero,
            end_slot: end_slot as u32,
            end_depth,
            end_slot_visited,
            grown: self.grown,
        }
    }
}

/// What the next contig needs from a finished contig to resolve the slot
/// samtools carries across the contig change.
#[derive(Debug, Clone)]
pub(crate) struct ContigCovSummary {
    tid: i32,
    /// Slots whose first visit saw a non-zero depth, sorted by slot.
    first_nonzero: Vec<(u32, u32)>,
    /// The slot left unflushed at the end of the contig, and its depth.
    end_slot: u32,
    end_depth: u32,
    /// Whether `end_slot` was visited before the final flush.
    end_slot_visited: bool,
    grown: bool,
}

/// Apply the carried slots between consecutive contigs (in BAM order, i.e.
/// ascending tid for a coordinate-sorted BAM) to `hist`. The depth left in the
/// buffer after the last contig is dropped, as in samtools.
pub(crate) fn resolve_carries(mut contigs: Vec<ContigCovSummary>, hist: &mut HashMap<u32, u64>) {
    contigs.sort_by_key(|c| c.tid);
    let mut carry: Option<(u32, u32, bool)> = None; // (slot, depth, grown)
    for contig in &mut contigs {
        if let Some((slot, depth, prev_grown)) = carry {
            // A grown buffer has a different slot layout; the carry cannot be
            // mapped (and samtools' own buffer resize loses pending depths).
            if depth > 0 && !prev_grown && !contig.grown {
                if slot == contig.end_slot && !contig.end_slot_visited {
                    contig.end_depth += depth;
                } else {
                    match contig
                        .first_nonzero
                        .binary_search_by_key(&slot, |&(s, _)| s)
                    {
                        Ok(i) => {
                            let first = contig.first_nonzero[i].1;
                            if let Some(n) = hist.get_mut(&first) {
                                *n -= 1;
                                if *n == 0 {
                                    hist.remove(&first);
                                }
                            }
                            *hist.entry(first + depth).or_insert(0) += 1;
                        }
                        Err(_) => *hist.entry(depth).or_insert(0) += 1,
                    }
                }
            }
        }
        carry = Some((contig.end_slot, contig.end_depth, contig.grown));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Straight port of samtools' round buffer over a whole BAM (one global
    /// buffer, contig changes flush with -1). Reads are `(tid, [(from, to)])`.
    fn reference_hist(reads: &[(i32, Vec<(i64, i64)>)]) -> HashMap<u32, u64> {
        let mut hist = HashMap::new();
        let mut buf = ContigCoverage::new(-1);
        buf.first_visit = None;
        buf.pos = 0; // calloc'd stats struct
        let mut tid = -1;
        for (t, blocks) in reads {
            if tid == -1 || tid != *t {
                buf.flush(-1, &mut hist);
                tid = *t;
            }
            buf.flush(blocks[0].0, &mut hist);
            for &(from, to) in blocks {
                buf.insert(from, to);
            }
        }
        buf.flush(-1, &mut hist);
        hist
    }

    /// The per-contig simulation plus carry resolution, as the accumulator uses it.
    fn split_hist(reads: &[(i32, Vec<(i64, i64)>)]) -> HashMap<u32, u64> {
        let mut hist = HashMap::new();
        let mut done = Vec::new();
        let mut cur: Option<ContigCoverage> = None;
        for (t, blocks) in reads {
            if cur.as_ref().is_some_and(|c| c.tid() != *t) {
                done.push(cur.take().unwrap().finish(&mut hist));
            }
            let c = cur.get_or_insert_with(|| ContigCoverage::new(*t));
            c.flush(blocks[0].0, &mut hist);
            for &(from, to) in blocks {
                c.insert(from, to);
            }
        }
        if let Some(c) = cur {
            done.push(c.finish(&mut hist));
        }
        resolve_carries(done, &mut hist);
        hist
    }

    #[test]
    fn last_slot_is_not_counted_at_end_of_input() {
        // A spliced read whose second block ends exactly on the last buffer
        // slot: samtools never counts that slot's position.
        let reads = vec![(0, vec![(100, 150), (100 + 1499, 100 + 1500)])];
        let hist = reference_hist(&reads);
        assert_eq!(hist.get(&1), Some(&50));
        assert_eq!(split_hist(&reads), hist);
    }

    #[test]
    fn carried_slot_is_counted_in_next_contig() {
        // Same read, then a second contig: the unflushed slot is carried and
        // counted while flushing the next contig.
        let reads = vec![
            (0, vec![(100, 150), (1599, 1600)]),
            (1, vec![(10, 60)]),
            (1, vec![(2000, 2050)]),
        ];
        let hist = reference_hist(&reads);
        assert_eq!(split_hist(&reads), hist);
        // 50 + 1 (carried) + 50 + 50 positions at depth 1.
        assert_eq!(hist.get(&1), Some(&151));
        assert_eq!(hist.len(), 1);
    }

    #[test]
    fn carry_adds_to_depth_in_next_contig() {
        // The carried slot lands on a position the next contig also covers.
        for shift in [0i64, 1, 5, 700, 1497, 1498, 1499, 1500, 3000] {
            let reads = vec![
                (0, vec![(0, 10), (1499, 1500)]),
                (2, vec![(shift, shift + 75)]),
                (2, vec![(shift + 3, shift + 78)]),
            ];
            assert_eq!(split_hist(&reads), reference_hist(&reads), "shift {shift}");
        }
    }

    #[test]
    fn wrapped_blocks_and_gaps_match_reference() {
        // Pseudo-random spliced reads on three contigs, including intron
        // lengths that wrap the buffer and gaps larger than the buffer.
        let mut state = 12345u64;
        let mut next = |m: u64| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 33) % m
        };
        let mut reads = Vec::new();
        for tid in 0..3 {
            let mut pos = next(3000) as i64;
            for _ in 0..2000 {
                let max_gap = if next(10) == 0 { 5000 } else { 40 };
                pos += next(max_gap) as i64;
                let len = 1 + next(75) as i64;
                let mut blocks = vec![(pos, pos + len)];
                if next(3) == 0 {
                    let intron = next(4000) as i64;
                    let s = pos + len + intron;
                    blocks.push((s, s + 1 + next(75) as i64));
                }
                reads.push((tid, blocks));
            }
        }
        assert_eq!(split_hist(&reads), reference_hist(&reads));
    }

    #[test]
    fn block_of_exactly_buffer_size_adds_nothing() {
        let mut c = ContigCoverage::new(0);
        let mut hist = HashMap::new();
        c.flush(0, &mut hist);
        c.insert(0, INITIAL_BUFFER_SIZE as i64);
        c.finish(&mut hist);
        assert!(hist.is_empty());
    }
}
