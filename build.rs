//! Build script — compiles the C++ RNG shim for preseq compatibility.
//!
//! The shim wraps `std::mt19937` and `std::binomial_distribution` so that
//! RustQC's preseq bootstrap uses the exact same random number generation
//! as upstream preseq compiled on the same platform.

fn main() {
    cc::Build::new()
        .cpp(true)
        .file("cpp/rng_shim.cpp")
        .std("c++17")
        .warnings(true)
        .compile("rng_shim");

    println!("cargo:rerun-if-changed=cpp/rng_shim.cpp");
}
