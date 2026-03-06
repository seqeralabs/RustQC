// C++ shim providing std::mt19937 + std::binomial_distribution via C API.
//
// This ensures RustQC's preseq bootstrap produces byte-identical results to
// upstream preseq, which uses the host C++ standard library's implementations.
// Different C++ standard libraries use different algorithms:
//   - libstdc++ (GCC/Linux): Devroye's rejection method
//   - libc++ (Apple Clang/macOS): Kemp's modal method
//
// By linking against the host's C++ stdlib at build time, RustQC automatically
// uses the same algorithm as a preseq binary compiled on the same platform.

#include <cstdint>
#include <random>

extern "C" {

/// Opaque handle to a std::mt19937 instance.
struct Mt19937Handle {
    std::mt19937 rng;
    explicit Mt19937Handle(uint32_t seed) : rng(seed) {}
};

/// Create a new MT19937 RNG with the given seed.
/// Returns an opaque pointer. Caller must free with mt19937_free().
Mt19937Handle* mt19937_new(uint32_t seed) {
    return new Mt19937Handle(seed);
}

/// Free a MT19937 RNG instance.
void mt19937_free(Mt19937Handle* handle) {
    delete handle;
}

/// Draw a single sample from Binomial(trials, p) using the host stdlib's
/// std::binomial_distribution implementation.
///
/// This matches exactly how upstream preseq's multinomial() function works:
/// a temporary binomial_distribution is constructed per call, then sampled once.
uint64_t binomial_sample(Mt19937Handle* handle, uint64_t trials, double p) {
    // Match preseq's exact usage: binomial_distribution<uint32_t>
    // However, for large trial counts we need uint64_t to avoid overflow.
    // preseq uses uint32_t but its histograms rarely exceed 2^32 distinct.
    //
    // We use the same uint32_t type as preseq when trials fit, otherwise
    // fall back to uint64_t for safety.
    if (trials <= UINT32_MAX) {
        std::binomial_distribution<uint32_t> dist(
            static_cast<uint32_t>(trials), p);
        return static_cast<uint64_t>(dist(handle->rng));
    } else {
        std::binomial_distribution<uint64_t> dist(trials, p);
        return dist(handle->rng);
    }
}

} // extern "C"
