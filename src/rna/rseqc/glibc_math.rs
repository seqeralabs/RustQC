//! Bit-exact ports of glibc's double-precision `log` and `exp`, FMA variant.
//!
//! `tin.py` takes its logarithms and exponentials from Python's `math`
//! module, which calls the C library. Rust's `f64::ln` and `f64::exp` call
//! the C library of the build target as well, so a binary built for macOS,
//! musl or another glibc configuration prints TINs that differ from `tin.py`'s
//! in the last digits for a small share of the transcripts, and the summary
//! statistics inherit those differences.
//!
//! Both functions are glibc 2.28+'s `log` and `exp`
//! (`sysdeps/ieee754/dbl-64/e_log.c` and `e_exp.c`, tables in `e_log_data.c`
//! and `e_exp_data.c`, taken from glibc 2.39) as compiled with FMA enabled:
//! `__ieee754_log_fma` and `__ieee754_exp_fma`. glibc on x86-64 selects those
//! through an ifunc on every CPU with FMA and AVX2 (anything since Intel
//! Haswell and AMD Zen), and the aarch64 build has FMA in its baseline and
//! returns the same bits (checked on 300 000 inputs). The generic build used
//! on older x86-64 CPUs rounds differently in the last bit for about 1 in
//! 10 000 inputs, so `tin.py` on such a CPU prints slightly different TINs;
//! its output is the one thing RustQC does not reproduce.
//!
//! The code below follows the C code operation by operation, with
//! [`f64::mul_add`] (a correctly rounded fused multiply-add on every target)
//! at exactly the places where GCC contracted `a * b + c` when it built
//! glibc's FMA variants. The contraction points were read from the machine
//! code of Ubuntu 24.04's `libm.so.6` (glibc 2.39), and each is noted next to
//! the code. Plain `*` and `+` are never contracted by Rust, so every other
//! operation rounds exactly once, as in the machine code.
//!
//! Only the arguments `compute_tin` can produce are handled here: `log` of a
//! positive finite number and `exp` of a finite number below 512 in magnitude
//! (an entropy of 512 nats would need more than 10^222 sampled positions).
//! Anything else is passed on to the platform's `ln` and `exp`.

use super::glibc_math_tables::{EXP_TAB, LOG_TAB};

const LN2_HI: f64 = f64::from_bits(0x3fe62e42fefa3800);
const LN2_LO: f64 = f64::from_bits(0x3d2ef35793c76730);

/// `__log_data.poly1`: polynomial for inputs close to 1.
const LOG_POLY1: [f64; 11] = [
    f64::from_bits(0xbfe0000000000000),
    f64::from_bits(0x3fd5555555555577),
    f64::from_bits(0xbfcffffffffffdcb),
    f64::from_bits(0x3fc999999995dd0c),
    f64::from_bits(0xbfc55555556745a7),
    f64::from_bits(0x3fc24924a344de30),
    f64::from_bits(0xbfbfffffa4423d65),
    f64::from_bits(0x3fbc7184282ad6ca),
    f64::from_bits(0xbfb999eb43b068ff),
    f64::from_bits(0x3fb78182f7afd085),
    f64::from_bits(0xbfb5521375d145cd),
];

/// `__log_data.poly`: polynomial for the table-driven path.
const LOG_POLY: [f64; 5] = [
    f64::from_bits(0xbfe0000000000001),
    f64::from_bits(0x3fd555555551305b),
    f64::from_bits(0xbfcfffffffeb4590),
    f64::from_bits(0x3fc999b324f10111),
    f64::from_bits(0xbfc55575e506c89f),
];

const EXP_INV_LN2_N: f64 = f64::from_bits(0x40671547652b82fe);
const EXP_NEG_LN2_HI_N: f64 = f64::from_bits(0xbf762e42fefa0000);
const EXP_NEG_LN2_LO_N: f64 = f64::from_bits(0xbd0cf79abc9e3b3a);
const EXP_SHIFT: f64 = f64::from_bits(0x4338000000000000);

/// `__exp_data.poly` (`C2` to `C5` in `e_exp.c`).
const EXP_POLY: [f64; 4] = [
    f64::from_bits(0x3fdffffffffffdbd),
    f64::from_bits(0x3fc555555555543c),
    f64::from_bits(0x3fa55555cf172b91),
    f64::from_bits(0x3f81111167a4d017),
];

/// glibc's `log(x)` (FMA build), for `x > 0`. Zero, negative, infinite and NaN
/// inputs are passed to the platform `ln`.
pub(crate) fn log(x: f64) -> f64 {
    let mut ix = x.to_bits();
    let top = (ix >> 48) as u32;

    const LO: u64 = 0x3fee000000000000; // 1.0 - 0x1p-4
    const HI: u64 = 0x3ff1090000000000; // 1.0 + 0x1.09p-4
    if ix.wrapping_sub(LO) < HI - LO {
        // Close to 1.0.
        if ix == 1.0f64.to_bits() {
            return 0.0;
        }
        let b = &LOG_POLY1;
        let r = x - 1.0;
        let r2 = r * r;
        let r3 = r * r2;
        // Horner steps of the degree-10 polynomial, each `a + r^j * b` fused.
        let p1 = r.mul_add(b[2], b[1]);
        let p2 = r.mul_add(b[5], b[4]);
        let p3 = r.mul_add(b[8], b[7]);
        let p1 = r2.mul_add(b[3], p1);
        let p2 = r2.mul_add(b[6], p2);
        let q = r2.mul_add(b[9], p3);
        let q = r3.mul_add(b[10], q);
        let q = q.mul_add(r3, p2);
        let q = q.mul_add(r3, p1);
        // r + r*2^27 and (r + w) - w are exact in either rounding mode.
        let w = r * f64::from_bits(0x41a0000000000000); // 0x1p27
        let rhi = r + w - w;
        let rlo = r - rhi;
        let rhi2 = rhi * rhi;
        let hi = rhi2.mul_add(b[0], r); // r + rhi^2 * b0
        let lo = rhi2.mul_add(b[0], r - hi); // (r - hi) + rhi^2 * b0
        let lo = (b[0] * rlo).mul_add(rhi + r, lo); // lo + b0 * rlo * (rhi + r)
        let y = q.mul_add(r3, lo); // lo + r^3 * poly
        return y + hi;
    }
    if top.wrapping_sub(0x0010) >= 0x7ff0 - 0x0010 {
        if ix == 0 || ix >> 52 != 0 {
            // Zero, negative, infinite or NaN.
            return x.ln();
        }
        // Subnormal: scale into the normal range, as glibc does.
        ix = (x * f64::from_bits(0x4330000000000000)).to_bits(); // 0x1p52
        ix = ix.wrapping_sub(52u64 << 52);
    }

    // x = 2^k z, with z in [OFF, 2 OFF), split into 128 subintervals.
    const OFF: u64 = 0x3fe6000000000000;
    let tmp = ix.wrapping_sub(OFF);
    let i = ((tmp >> (52 - 7)) % 128) as usize;
    let k = (tmp as i64) >> 52;
    let iz = ix.wrapping_sub(tmp & (0xfffu64 << 52));
    let invc = f64::from_bits(LOG_TAB[i].0);
    let logc = f64::from_bits(LOG_TAB[i].1);
    let z = f64::from_bits(iz);

    // r ~= z/c - 1, rounded once.
    let r = z.mul_add(invc, -1.0);
    let kd = k as f64;

    // hi + lo = r + log(c) + k*Ln2.
    let w = kd.mul_add(LN2_HI, logc);
    let hi = r + w;
    let lo = kd.mul_add(LN2_LO, w - hi + r);

    // log(x) = lo + (log1p(r) - r) + hi.
    let a = &LOG_POLY;
    let r2 = r * r;
    let r3 = r * r2;
    let t1 = r.mul_add(a[2], a[1]);
    let t2 = r.mul_add(a[4], a[3]);
    let lo = r2.mul_add(a[0], lo);
    let poly = t2.mul_add(r2, t1);
    r3.mul_add(poly, lo) + hi
}

/// glibc's `exp(x)` (FMA build), for |x| < 512. Larger magnitudes and
/// non-finite inputs are passed to the platform `exp`.
pub(crate) fn exp(x: f64) -> f64 {
    let abstop = ((x.to_bits() >> 52) & 0x7ff) as u32;
    const TOP_TINY: u32 = 0x3c9; // top 12 bits of 0x1p-54
    const TOP_512: u32 = 0x408; // top 12 bits of 512.0
    if abstop.wrapping_sub(TOP_TINY) >= TOP_512 - TOP_TINY {
        if abstop.wrapping_sub(TOP_TINY) >= 0x8000_0000 {
            // Tiny |x| (including 0): 1 + x.
            return 1.0 + x;
        }
        return x.exp();
    }

    // exp(x) = 2^(k/N) * exp(r), with N = 128 and r in [-ln2/2N, ln2/2N].
    let mut kd = x.mul_add(EXP_INV_LN2_N, EXP_SHIFT); // InvLn2N * x + Shift
    let ki = kd.to_bits();
    kd -= EXP_SHIFT;
    let r = kd.mul_add(EXP_NEG_LN2_HI_N, x);
    let r = kd.mul_add(EXP_NEG_LN2_LO_N, r);
    // 2^(k/N) ~= scale * (1 + tail).
    let idx = 2 * (ki % 128) as usize;
    let top = ki << (52 - 7);
    let tail = f64::from_bits(EXP_TAB[idx]);
    let sbits = EXP_TAB[idx + 1].wrapping_add(top);
    let r2 = r * r;
    let c = &EXP_POLY;
    let ta = r.mul_add(c[1], c[0]);
    let tb = r.mul_add(c[3], c[2]);
    let tmp = ta.mul_add(r2, tail + r);
    let tmp = (r2 * r2).mul_add(tb, tmp);
    let scale = f64::from_bits(sbits);
    scale.mul_add(tmp, scale)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Inputs and the results of glibc 2.39's `log` and `exp` on a CPU with FMA
    // (what Python's `math.log` and `math.exp` return there), as IEEE-754 bit
    // patterns. Entries marked "generic build differs" are inputs for which
    // glibc's generic (non-FMA) build returns a neighbouring double. The
    // subnormal and tiny inputs exercise the special paths.
    const LOG_VECTORS: &[(u64, u64)] = &[
        (0x3fecec2bb38afa34, 0xbfb9e2a5c7d0b764), // generic build differs
        (0x3fe65be5b00e8522, 0xbfd6f1cc39044bf3), // generic build differs
        (0x3ff14ee7e9f500a4, 0x3fb41eae99645126), // generic build differs
        (0x3fe547332e20495c, 0xbfda1d9874251361), // generic build differs
        (0x3fed9d33141fad8d, 0xbfb3d629f7ed4440), // generic build differs
        (0x3fe9ccc8547ab8a2, 0xbfcb910ca183e7c0), // generic build differs
        (0x3fed21b5116f5de0, 0xbfb80a7ced5115d4), // generic build differs
        (0x3fecdc66cf1f9d9b, 0xbfba6e5f562538b3), // generic build differs
        (0x3fee3c402fc0214a, 0xbfad0aed2d115496), // generic build differs
        (0x3fee5c36d76842ef, 0xbfaaeec6e62e1b9d), // generic build differs
        (0x3fee1ef3e671abc4, 0xbfaefbfea4e6a3a8), // generic build differs
        (0x3ff0c7e4dfeaefc8, 0x3fa8656c944572f0), // generic build differs (near 1)
        (0x3fee0b57736386d2, 0xbfb024e2df71533f), // generic build differs (near 1)
        (0x3fee8269f4f4327e, 0xbfa86c261cee6484), // generic build differs (near 1)
        (0x3fd5935f686ee044, 0xbff165a612420b4c),
        (0x3fe5eede501e021a, 0xbfd82ce40a36cc16),
        (0x369c72d601975bb2, 0xc059d9615d3293bd),
        (0x7604fd61992d17f0, 0x4082bec004044773),
        (0x3e6dab7233c0f7cc, 0xc030b60d5f12f4f7),
        (0x3fefffffffff562c, 0xbd953a8000003855),
        (0x3feff1ac5d6874da, 0xbf5cadb122b0b270),
        (0x3ff0e6f463521c65, 0x3fac15b44c0b260a),
        (0x000ee1780f9c0e82, 0xc08623c06369870e),
        (0x0007e6de644a3e99, 0xc08628d0b7a2ee4b),
        (0x000fc3341b23eb93, 0xc086234a7dc1f43c),
    ];
    const EXP_VECTORS: &[(u64, u64)] = &[
        (0xc0485bfef4264fc0, 0x3b8a3dadb406b785), // generic build differs
        (0x3fee107bbde0d838, 0x40047848e4d06db1), // generic build differs
        (0x40392d7e01f530c4, 0x4234067712417e53), // generic build differs
        (0xc07c6f6b1d36e92e, 0x16e8b223988e78a1), // generic build differs
        (0xbfef168fbcec9eb0, 0x3fd8398fd1a4a781), // generic build differs
        (0xc070e8d96829ff3e, 0x27898903d421ce21), // generic build differs
        (0x402ffa60ee3abc5a, 0x4160c38b3056e4b9), // generic build differs
        (0x40382d9f21bef4d0, 0x421d7b7f90b1f779), // generic build differs
        (0xbfeb18418c2966e0, 0x3fdb71d62f7503a8),
        (0x4068d4fa8ae2a6bc, 0x51d83de75a00a782),
        (0xbf868915c87e7f00, 0x3fefa65a284488a7),
        (0x3ffe3a4bc1c4c1f0, 0x401a750857119ada),
        (0x40223e02dd073db4, 0x40c1dd3326870d7c),
        (0xbbf89cf39c0eeeb0, 0x3ff0000000000000),
        (0xa7cc3ef0cbb3af1c, 0x3ff0000000000000),
    ];

    #[test]
    fn test_log_matches_glibc() {
        for &(x, want) in LOG_VECTORS {
            let x = f64::from_bits(x);
            assert_eq!(log(x).to_bits(), want, "log({x:e})");
        }
    }

    #[test]
    fn test_exp_matches_glibc() {
        for &(x, want) in EXP_VECTORS {
            let x = f64::from_bits(x);
            assert_eq!(exp(x).to_bits(), want, "exp({x:e})");
        }
    }

    #[test]
    fn test_special_values() {
        assert_eq!(log(1.0).to_bits(), 0.0f64.to_bits());
        assert_eq!(exp(0.0), 1.0);
        assert_eq!(exp(-0.0), 1.0);
        assert_eq!(exp(1e-300), 1.0);
        assert!(log(0.0).is_infinite());
        assert!(log(-1.0).is_nan());
        assert!(exp(1000.0).is_infinite());
        assert_eq!(exp(-1000.0), 0.0);
    }

    #[test]
    fn test_close_to_platform_libm() {
        // Both are correct to within an ulp or so; this only guards
        // against a transcription error in a table or coefficient.
        for i in 1..20_000u32 {
            let x = f64::from(i) * 0.000_731;
            let l = log(x);
            assert!(
                (l - x.ln()).abs() <= 4.0 * f64::EPSILON * l.abs().max(1e-3),
                "log({x})"
            );
            let e = exp(x - 7.0);
            assert!(
                (e / (x - 7.0).exp() - 1.0).abs() <= 4.0 * f64::EPSILON,
                "exp({})",
                x - 7.0
            );
        }
    }
}
