//! FFI bindings to the C++ RNG shim for preseq-compatible random sampling.
//!
//! This module wraps `std::mt19937` and `std::binomial_distribution` from the
//! host C++ standard library. By using the exact same RNG and distribution
//! implementations as upstream preseq, the bootstrap resampling produces
//! byte-identical results when compiled on the same platform.
//!
//! ## Why FFI instead of a pure-Rust implementation?
//!
//! Different C++ standard libraries use different algorithms for
//! `std::binomial_distribution`:
//!
//! - **libstdc++ (GCC/Linux)**: Devroye's 4-region rejection method with a
//!   normal distribution proposal, falling back to a waiting-time method for
//!   small `t*p`. Uses `std::normal_distribution` internally (Marsaglia polar
//!   method with cached state).
//!
//! - **libc++ (Apple Clang/macOS)**: Kemp's modal method — a single uniform
//!   draw followed by a walk outward from the mode.
//!
//! These algorithms consume different numbers of RNG values per sample, so even
//! with identical seeds and identical MT19937 state, the bootstrap diverges
//! immediately. Rather than reimplementing a specific algorithm (which would
//! only match one platform), we link against the host's C++ stdlib so that
//! RustQC automatically uses whichever algorithm preseq would use on that
//! platform.

// Raw FFI declarations matching cpp/rng_shim.cpp
extern "C" {
    fn mt19937_new(seed: u32) -> *mut std::ffi::c_void;
    fn mt19937_free(handle: *mut std::ffi::c_void);
    fn binomial_sample(handle: *mut std::ffi::c_void, trials: u64, p: f64) -> u64;
}

/// A Mersenne Twister 19937 (32-bit) RNG backed by the C++ standard library.
///
/// This is a thin wrapper around `std::mt19937`. The RNG state lives on the
/// C++ heap and is freed on drop.
pub struct CppMt19937 {
    handle: *mut std::ffi::c_void,
}

// SAFETY: The C++ mt19937 handle is only accessed through &mut self methods,
// and CppMt19937 is never shared across threads (preseq bootstrap is sequential).
unsafe impl Send for CppMt19937 {}

impl CppMt19937 {
    /// Create a new MT19937 RNG with the given seed.
    ///
    /// This matches `std::mt19937 rng(seed)` in C++.
    pub fn new(seed: u32) -> Self {
        let handle = unsafe { mt19937_new(seed) };
        assert!(!handle.is_null(), "mt19937_new returned null");
        Self { handle }
    }

    /// Draw a single sample from `Binomial(trials, p)`.
    ///
    /// Uses `std::binomial_distribution<uint32_t>` (or `uint64_t` for large
    /// trial counts), matching upstream preseq's `multinomial()` function
    /// which constructs a temporary distribution per call.
    pub fn binomial(&mut self, trials: u64, p: f64) -> u64 {
        unsafe { binomial_sample(self.handle, trials, p) }
    }
}

impl Drop for CppMt19937 {
    fn drop(&mut self) {
        unsafe { mt19937_free(self.handle) };
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cpp_mt19937_creates_and_drops() {
        let _rng = CppMt19937::new(42);
        // Should not crash or leak
    }

    #[test]
    fn test_binomial_deterministic() {
        // Same seed should produce same sequence
        let mut rng1 = CppMt19937::new(1);
        let mut rng2 = CppMt19937::new(1);

        for _ in 0..20 {
            let a = rng1.binomial(1000, 0.3);
            let b = rng2.binomial(1000, 0.3);
            assert_eq!(a, b, "Same seed must produce identical samples");
        }
    }

    #[test]
    fn test_binomial_edge_cases() {
        let mut rng = CppMt19937::new(1);

        // p = 0 should always return 0
        for _ in 0..10 {
            assert_eq!(rng.binomial(1000, 0.0), 0);
        }

        // p = 1 should always return trials
        for _ in 0..10 {
            assert_eq!(rng.binomial(1000, 1.0), 1000);
        }

        // trials = 0 should always return 0
        for _ in 0..10 {
            assert_eq!(rng.binomial(0, 0.5), 0);
        }
    }

    #[test]
    fn test_binomial_range() {
        let mut rng = CppMt19937::new(42);

        for _ in 0..100 {
            let trials = 500;
            let sample = rng.binomial(trials, 0.5);
            assert!(sample <= trials, "Sample {} > trials {}", sample, trials);
        }
    }
}
