//! Build script — compiles the C++ RNG shim for preseq compatibility and
//! embeds git/build metadata as compile-time environment variables.
//!
//! The shim wraps `std::mt19937` and `std::binomial_distribution` so that
//! RustQC's preseq bootstrap uses the exact same random number generation
//! as upstream preseq compiled on the same platform.
//!
//! Embedded variables:
//! - `GIT_SHORT_HASH` — short commit hash (e.g. `84ec57f`), or `unknown`
//! - `BUILD_TIMESTAMP` — UTC timestamp of the build (e.g. `2026-03-07T12:34:56Z`)

use std::process::Command;

fn main() {
    // --- C++ RNG shim ---
    cc::Build::new()
        .cpp(true)
        .file("cpp/rng_shim.cpp")
        .std("c++17")
        .warnings(true)
        .compile("rng_shim");

    println!("cargo:rerun-if-changed=cpp/rng_shim.cpp");

    // --- Git short hash ---
    // First check the GIT_SHORT_HASH env var (set via Docker build arg),
    // then fall back to running `git rev-parse` (works in local/CI builds).
    let git_hash = std::env::var("GIT_SHORT_HASH")
        .ok()
        .filter(|s| !s.is_empty() && s != "unknown")
        .unwrap_or_else(|| {
            Command::new("git")
                .args(["rev-parse", "--short", "HEAD"])
                .output()
                .ok()
                .filter(|o| o.status.success())
                .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
                .unwrap_or_else(|| "unknown".to_string())
        });
    println!("cargo:rustc-env=GIT_SHORT_HASH={}", git_hash);

    // --- Build timestamp (UTC, ISO-8601) ---
    let build_ts = Command::new("date")
        .args(["-u", "+%Y-%m-%dT%H:%M:%SZ"])
        .output()
        .ok()
        .filter(|o| o.status.success())
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
        .unwrap_or_else(|| "unknown".to_string());
    println!("cargo:rustc-env=BUILD_TIMESTAMP={}", build_ts);

    // Rebuild when HEAD changes (new commits)
    println!("cargo:rerun-if-changed=.git/HEAD");
}
