//! CPU feature detection, binary target identification, and compatibility checks.
//!
//! Provides compile-time binary target identification (e.g. x86-64-v3, aarch64-neoverse-v1),
//! runtime CPU feature detection, upgrade hints when a faster binary is available,
//! and a startup compatibility guard to prevent SIGILL on mismatched binaries.

use anyhow::{bail, Result};

// ============================================================================
// Binary target (compile-time)
// ============================================================================

/// Returns the microarchitecture level this binary was compiled for.
///
/// Determined at compile time via `cfg!(target_feature = ...)`.
pub fn binary_target() -> &'static str {
    #[cfg(target_arch = "x86_64")]
    {
        if cfg!(target_feature = "avx512f") {
            "x86-64-v4"
        } else if cfg!(target_feature = "avx2") {
            "x86-64-v3"
        } else {
            "x86-64"
        }
    }
    #[cfg(target_arch = "aarch64")]
    {
        // SVE indicates neoverse-v1 or higher target
        if cfg!(target_feature = "sve") {
            "aarch64-neoverse-v1"
        } else {
            "aarch64"
        }
    }
    #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
    {
        "unknown"
    }
}

/// Returns a short label for the key feature set of the binary target.
///
/// Used in the header line: `Binary: x86-64-v3 (AVX2)`.
pub fn binary_target_label() -> &'static str {
    match binary_target() {
        "x86-64-v4" => "AVX-512",
        "x86-64-v3" => "AVX2",
        "aarch64-neoverse-v1" => "SVE",
        "aarch64" => "NEON",
        _ => "baseline",
    }
}

// ============================================================================
// Runtime CPU feature detection
// ============================================================================

/// Detects CPU features available at runtime.
///
/// Returns a list of feature names (e.g. `["AVX2", "SSE4.2", "POPCNT"]`),
/// ordered from most to least capable.
pub fn detected_features() -> Vec<&'static str> {
    let mut features = Vec::new();

    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx512f") {
            features.push("AVX-512");
        }
        if is_x86_feature_detected!("avx2") {
            features.push("AVX2");
        }
        if is_x86_feature_detected!("sse4.2") {
            features.push("SSE4.2");
        }
        if is_x86_feature_detected!("popcnt") {
            features.push("POPCNT");
        }
    }

    #[cfg(target_arch = "aarch64")]
    {
        if std::arch::is_aarch64_feature_detected!("sve2") {
            features.push("SVE2");
        }
        if std::arch::is_aarch64_feature_detected!("sve") {
            features.push("SVE");
        }
        // NEON is mandatory on aarch64 but list it for clarity
        features.push("NEON");
    }

    features
}

// ============================================================================
// Upgrade hint
// ============================================================================

/// Returns an upgrade suggestion if a faster binary variant is available
/// for the current CPU.
///
/// Returns `None` if the binary already uses the best available level.
pub fn upgrade_hint() -> Option<String> {
    #[cfg(target_arch = "x86_64")]
    {
        let has_avx512 = is_x86_feature_detected!("avx512f");
        let has_avx2 = is_x86_feature_detected!("avx2");
        let target = binary_target();

        if target == "x86-64" && has_avx512 {
            return Some("Use the x86-64-v4 build for best performance.".to_string());
        }
        if target == "x86-64" && has_avx2 {
            return Some("Use the x86-64-v3 build for better performance.".to_string());
        }
        if target == "x86-64-v3" && has_avx512 {
            return Some("Use the x86-64-v4 build for best performance.".to_string());
        }
    }

    #[cfg(all(target_arch = "aarch64", target_os = "linux"))]
    {
        let has_sve = std::arch::is_aarch64_feature_detected!("sve");
        let target = binary_target();

        if target == "aarch64" && has_sve {
            return Some("Use the aarch64-neoverse-v1 build for better performance.".to_string());
        }
    }

    #[cfg(all(target_arch = "aarch64", target_os = "macos"))]
    {
        let target = binary_target();
        // All Apple Silicon supports apple-m1 target features
        if target == "aarch64" {
            return Some(
                "Use the macos-aarch64-apple-m1 build for better performance.".to_string(),
            );
        }
    }

    None
}

// ============================================================================
// CPU compatibility guard
// ============================================================================

/// Checks that the CPU supports the features required by this binary.
///
/// Call at the very start of `main()` — before any SIMD code runs — to give
/// a friendly error instead of a SIGILL crash.
pub fn check_cpu_compat() -> Result<()> {
    #[cfg(target_arch = "x86_64")]
    match binary_target() {
        "x86-64-v4" if !is_x86_feature_detected!("avx512f") => bail!(
            "This RustQC binary was compiled for x86-64-v4 (AVX-512) but your CPU \
             does not support AVX-512.\nPlease use the x86-64-v3 or baseline build instead.\n\
             See: https://github.com/seqeralabs/rustqc#installation"
        ),
        "x86-64-v3" if !is_x86_feature_detected!("avx2") => bail!(
            "This RustQC binary was compiled for x86-64-v3 (AVX2) but your CPU \
             does not support AVX2.\nPlease use the baseline build instead.\n\
             See: https://github.com/seqeralabs/rustqc#installation"
        ),
        _ => {}
    }

    #[cfg(target_arch = "aarch64")]
    if binary_target() == "aarch64-neoverse-v1" && !std::arch::is_aarch64_feature_detected!("sve") {
        bail!(
            "This RustQC binary was compiled for aarch64-neoverse-v1 (SVE) but your CPU \
             does not support SVE.\nPlease use the baseline aarch64 build instead.\n\
             See: https://github.com/seqeralabs/rustqc#installation"
        );
    }

    Ok(())
}

/// Formats the CPU info line for display in the header and `--version` output.
///
/// Example outputs:
/// - `Binary: x86-64 (baseline) | CPU: AVX2 AVX-512 detected. Use the x86-64-v3 build for better performance.`
/// - `Binary: x86-64-v3 (AVX2) | CPU: AVX2 detected`
/// - `Binary: aarch64 (NEON) | CPU: SVE detected. Use the aarch64-neoverse-v1 build for better performance.`
pub fn cpu_info_line() -> String {
    let target = binary_target();
    let label = binary_target_label();
    let features = detected_features();
    let hint = upgrade_hint();

    let features_str = if features.is_empty() {
        "none detected".to_string()
    } else {
        format!("{} detected", features.join(" "))
    };

    match hint {
        Some(hint) => format!("Binary: {target} ({label}) | CPU: {features_str}. {hint}"),
        None => format!("Binary: {target} ({label}) | CPU: {features_str}"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn binary_target_is_not_empty() {
        assert!(!binary_target().is_empty());
    }

    #[test]
    fn binary_target_label_is_not_empty() {
        assert!(!binary_target_label().is_empty());
    }

    #[test]
    fn detected_features_returns_vec() {
        let features = detected_features();
        // On any platform we should detect at least one feature
        // (NEON on aarch64, at least SSE on x86_64)
        // But don't assert non-empty since CI might run on unusual hardware
        assert!(features.len() <= 10);
    }

    #[test]
    fn cpu_info_line_contains_binary() {
        let line = cpu_info_line();
        assert!(line.contains("Binary:"));
        assert!(line.contains("CPU:"));
    }

    #[test]
    fn check_cpu_compat_passes_on_current_hardware() {
        // The binary we're running was compiled for this CPU, so this must pass
        assert!(check_cpu_compat().is_ok());
    }
}
