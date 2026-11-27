//! SAM/BAM header wrapper.

/// SAM header wrapper with rust-htslib-compatible helpers.
#[derive(Debug, Clone)]
pub struct Header {
    pub(crate) inner: noodles::sam::Header,
    ref_names: Vec<String>,
    ref_lengths: Vec<u64>,
}

impl Header {
    /// Wrap a noodles SAM header.
    pub fn from_noodles(inner: noodles::sam::Header) -> Self {
        let mut ref_names = Vec::new();
        let mut ref_lengths = Vec::new();
        for (name, map) in inner.reference_sequences() {
            ref_names.push(String::from_utf8_lossy(name.as_ref()).into_owned());
            let len = u64::try_from(usize::from(map.length())).unwrap_or(0);
            ref_lengths.push(len);
        }
        Self {
            inner,
            ref_names,
            ref_lengths,
        }
    }

    /// Create an empty header (rust-htslib-compatible constructor).
    pub fn new() -> Self {
        Self::empty()
    }

    /// Create an empty header.
    pub fn empty() -> Self {
        Self::from_noodles(noodles::sam::Header::default())
    }

    /// Number of reference sequences.
    pub fn target_count(&self) -> u32 {
        self.ref_names.len() as u32
    }

    /// Reference name for a target ID.
    pub fn tid2name(&self, tid: u32) -> &[u8] {
        self.ref_names
            .get(tid as usize)
            .map(|s| s.as_bytes())
            .unwrap_or(b"*")
    }

    /// Reference length for a target ID.
    pub fn target_len(&self, tid: u32) -> Option<u64> {
        self.ref_lengths.get(tid as usize).copied()
    }

    /// Borrow the underlying noodles header.
    pub fn noodles_header(&self) -> &noodles::sam::Header {
        &self.inner
    }
}

impl Default for Header {
    fn default() -> Self {
        Self::empty()
    }
}
