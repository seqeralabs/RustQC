//! Pure-Rust emulation of CPython's `random` module, for RSeQC reproducibility.
//!
//! RSeQC's `junction_saturation.py` shuffles its splice-junction observations
//! with Python's module-level `random.shuffle()`. To reproduce its output
//! exactly (including the intermediate sampling points), RustQC must draw the
//! same random stream as CPython after `random.seed(n)`. This module
//! reimplements the parts of CPython's `random.Random` that path uses:
//!
//! - `seed(int)`: MT19937 seeded with `init_by_array()` over the 32-bit words
//!   of `abs(n)`, least significant word first (`Modules/_randommodule.c`,
//!   `random_seed`).
//! - `getrandbits(k)`: the top `k` bits of one MT19937 output for `k <= 32`;
//!   for larger `k`, 32-bit words filled from least to most significant, with
//!   the last word's low bits dropped.
//! - `_randbelow(n)`: `getrandbits(n.bit_length())` with rejection sampling
//!   (`Lib/random.py`, `_randbelow_with_getrandbits`).
//! - `shuffle(x)`: Fisher-Yates from the end of the list, swapping `x[i]` with
//!   `x[_randbelow(i + 1)]`.
//!
//! These algorithms have been stable across CPython 3.2 onwards; the unit
//! tests check them against values printed by CPython 3.12.

/// MT19937 state size.
const N: usize = 624;
/// MT19937 shift size.
const M: usize = 397;
const MATRIX_A: u32 = 0x9908_b0df;
const UPPER_MASK: u32 = 0x8000_0000;
const LOWER_MASK: u32 = 0x7fff_ffff;

/// A `random.Random` instance as CPython implements it.
#[derive(Clone)]
pub struct PyRandom {
    mt: Box<[u32; N]>,
    mti: usize,
}

impl std::fmt::Debug for PyRandom {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PyRandom").field("mti", &self.mti).finish()
    }
}

impl PyRandom {
    /// Equivalent to `random.Random(seed)` / `random.seed(seed)` for a
    /// non-negative integer seed.
    pub fn new(seed: u64) -> Self {
        // CPython splits abs(seed) into 32-bit words, least significant first,
        // using at least one word (so seed 0 is the key [0]).
        let key: Vec<u32> = if seed >> 32 == 0 {
            vec![seed as u32]
        } else {
            vec![seed as u32, (seed >> 32) as u32]
        };
        let mut rng = PyRandom {
            mt: Box::new([0u32; N]),
            mti: N + 1,
        };
        rng.init_by_array(&key);
        rng
    }

    fn init_genrand(&mut self, s: u32) {
        let mt = &mut self.mt;
        mt[0] = s;
        for i in 1..N {
            mt[i] = 1_812_433_253u32
                .wrapping_mul(mt[i - 1] ^ (mt[i - 1] >> 30))
                .wrapping_add(i as u32);
        }
        self.mti = N;
    }

    fn init_by_array(&mut self, key: &[u32]) {
        self.init_genrand(19_650_218);
        let mt = &mut self.mt;
        let mut i = 1usize;
        let mut j = 0usize;
        for _ in 0..N.max(key.len()) {
            mt[i] = (mt[i] ^ (mt[i - 1] ^ (mt[i - 1] >> 30)).wrapping_mul(1_664_525))
                .wrapping_add(key[j])
                .wrapping_add(j as u32);
            i += 1;
            j += 1;
            if i >= N {
                mt[0] = mt[N - 1];
                i = 1;
            }
            if j >= key.len() {
                j = 0;
            }
        }
        for _ in 0..N - 1 {
            mt[i] = (mt[i] ^ (mt[i - 1] ^ (mt[i - 1] >> 30)).wrapping_mul(1_566_083_941))
                .wrapping_sub(i as u32);
            i += 1;
            if i >= N {
                mt[0] = mt[N - 1];
                i = 1;
            }
        }
        mt[0] = 0x8000_0000;
    }

    /// One raw 32-bit MT19937 output (`genrand_uint32`).
    pub fn genrand_u32(&mut self) -> u32 {
        if self.mti >= N {
            let mt = &mut self.mt;
            let mag01 = |y: u32| if y & 1 == 0 { 0 } else { MATRIX_A };
            for kk in 0..N - M {
                let y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
                mt[kk] = mt[kk + M] ^ (y >> 1) ^ mag01(y);
            }
            for kk in N - M..N - 1 {
                let y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
                mt[kk] = mt[kk + M - N] ^ (y >> 1) ^ mag01(y);
            }
            let y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
            mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ mag01(y);
            self.mti = 0;
        }
        let mut y = self.mt[self.mti];
        self.mti += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^= y >> 18;
        y
    }

    /// `random.getrandbits(k)` for `0 <= k <= 64`.
    pub fn getrandbits(&mut self, k: u32) -> u64 {
        assert!(k <= 64, "getrandbits supports at most 64 bits, got {k}");
        if k == 0 {
            0
        } else if k <= 32 {
            u64::from(self.genrand_u32() >> (32 - k))
        } else {
            // Least significant word first; the last word keeps its top bits.
            let lo = u64::from(self.genrand_u32());
            let hi = u64::from(self.genrand_u32() >> (64 - k));
            (hi << 32) | lo
        }
    }

    /// `random._randbelow(n)`: a uniform integer in `[0, n)`, for `n > 0`.
    pub fn randbelow(&mut self, n: u64) -> u64 {
        assert!(n > 0, "randbelow requires n > 0");
        let k = u64::BITS - n.leading_zeros();
        let mut r = self.getrandbits(k);
        while r >= n {
            r = self.getrandbits(k);
        }
        r
    }

    /// `random.shuffle(x)`: shuffle a slice in place.
    pub fn shuffle<T>(&mut self, x: &mut [T]) {
        for i in (1..x.len()).rev() {
            let j = self.randbelow(i as u64 + 1) as usize;
            x.swap(i, j);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Expected values were printed by CPython 3.12.3:
    //   r = random.Random(seed)
    //   [r.getrandbits(32) for _ in range(3)], r.getrandbits(7),
    //   r.getrandbits(33), r.getrandbits(64),
    //   [r._randbelow(n) for n in (1, 2, 3, 1000, 2**31, 2**32 + 1)]
    //   x = list(range(20)); r.shuffle(x)
    struct Fixture {
        seed: u64,
        u32s: [u64; 3],
        b7: u64,
        b33: u64,
        b64: u64,
        below: [u64; 6],
        shuffle20: [usize; 20],
    }

    const FIXTURES: &[Fixture] = &[
        Fixture {
            seed: 0,
            u32s: [3626764237, 1654615998, 3255389356],
            b7: 113,
            b33: 1806341205,
            b64: 17809683713383489082,
            below: [0, 1, 1, 991, 2046968324, 598176026],
            shuffle20: [
                0, 13, 16, 15, 6, 18, 2, 12, 17, 7, 5, 10, 14, 11, 1, 9, 19, 8, 3, 4,
            ],
        },
        Fixture {
            seed: 1,
            u32s: [577090037, 2444712010, 3639700191],
            b7: 102,
            b33: 3280387012,
            b64: 2175216119781798972,
            below: [0, 1, 1, 667, 1630434966, 3387541014],
            shuffle20: [
                10, 8, 6, 14, 5, 11, 2, 1, 19, 4, 7, 17, 18, 16, 9, 13, 12, 0, 15, 3,
            ],
        },
        Fixture {
            seed: 42,
            u32s: [2746317213, 478163327, 107420369],
            b7: 94,
            b33: 1181241943,
            b64: 2574020394472462046,
            below: [0, 0, 2, 432, 136505587, 127978094],
            shuffle20: [
                5, 2, 12, 9, 1, 17, 13, 4, 18, 3, 15, 8, 14, 10, 11, 19, 0, 16, 7, 6,
            ],
        },
        Fixture {
            seed: 4294967295,
            u32s: [2728839433, 2661025012, 872737089],
            b7: 27,
            b33: 6903586996,
            b64: 11124323167451281159,
            below: [0, 1, 2, 705, 1201118819, 2135565362],
            shuffle20: [
                8, 1, 3, 9, 14, 10, 17, 7, 18, 12, 5, 0, 13, 16, 19, 11, 15, 4, 6, 2,
            ],
        },
        Fixture {
            seed: 4294967296,
            u32s: [485306839, 1508871100, 1794561286],
            b7: 119,
            b33: 71624475,
            b64: 16161065642225157897,
            below: [0, 1, 0, 141, 325445696, 4109391884],
            shuffle20: [
                3, 6, 0, 5, 17, 13, 8, 9, 19, 4, 1, 11, 2, 7, 14, 15, 12, 18, 10, 16,
            ],
        },
        Fixture {
            seed: u64::MAX,
            u32s: [93740670, 1068495656, 1452108352],
            b7: 79,
            b33: 910393425,
            b64: 1752763915482752451,
            below: [0, 0, 2, 197, 2017374292, 4079183629],
            shuffle20: [
                3, 19, 0, 16, 2, 18, 6, 14, 7, 11, 12, 9, 13, 4, 5, 1, 17, 8, 15, 10,
            ],
        },
    ];

    #[test]
    fn test_matches_cpython_stream() {
        for fx in FIXTURES {
            let mut r = PyRandom::new(fx.seed);
            let u32s: Vec<u64> = (0..3).map(|_| r.getrandbits(32)).collect();
            assert_eq!(u32s, fx.u32s, "getrandbits(32), seed {}", fx.seed);
            assert_eq!(r.getrandbits(7), fx.b7, "getrandbits(7), seed {}", fx.seed);
            assert_eq!(
                r.getrandbits(33),
                fx.b33,
                "getrandbits(33), seed {}",
                fx.seed
            );
            assert_eq!(
                r.getrandbits(64),
                fx.b64,
                "getrandbits(64), seed {}",
                fx.seed
            );
            let below: Vec<u64> = [1, 2, 3, 1000, 1 << 31, (1 << 32) + 1]
                .iter()
                .map(|&n| r.randbelow(n))
                .collect();
            assert_eq!(below, fx.below, "_randbelow, seed {}", fx.seed);
            let mut x: Vec<usize> = (0..20).collect();
            r.shuffle(&mut x);
            assert_eq!(x, fx.shuffle20, "shuffle, seed {}", fx.seed);
        }
    }

    #[test]
    fn test_large_shuffle_matches_cpython() {
        // y = list(range(100000)); random.Random(seed).shuffle(y)
        // y[:8], y[-4:], sum(i * v for i, v in enumerate(y)) % (2**61 - 1)
        let cases: [(u64, [usize; 8], [usize; 4], u128); 3] = [
            (
                0,
                [74921, 6849, 99252, 14859, 14403, 24051, 23134, 28153],
                [5306, 55125, 99346, 50494],
                249610304483173,
            ),
            (
                1,
                [63229, 7871, 74587, 59445, 71381, 5404, 56721, 41863],
                [33432, 8271, 74606, 17611],
                250156808791575,
            ),
            (
                42,
                [3813, 99183, 80737, 65674, 53670, 97104, 38743, 36888],
                [97196, 3278, 14592, 83810],
                250177806368157,
            ),
        ];
        for (seed, head, tail, wsum) in cases {
            let mut y: Vec<usize> = (0..100_000).collect();
            PyRandom::new(seed).shuffle(&mut y);
            assert_eq!(&y[..8], &head, "head, seed {seed}");
            assert_eq!(&y[y.len() - 4..], &tail, "tail, seed {seed}");
            let sum: u128 = y
                .iter()
                .enumerate()
                .map(|(i, &v)| i as u128 * v as u128)
                .sum::<u128>()
                % ((1u128 << 61) - 1);
            assert_eq!(sum, wsum, "weighted sum, seed {seed}");
        }
    }

    #[test]
    fn test_shuffle_trivial_lengths() {
        let mut r = PyRandom::new(7);
        let mut empty: Vec<u8> = Vec::new();
        r.shuffle(&mut empty);
        let mut one = vec![9u8];
        r.shuffle(&mut one);
        assert_eq!(one, vec![9]);
        assert_eq!(r.getrandbits(0), 0);
    }
}
