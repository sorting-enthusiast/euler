use std::time::Instant;

use super::prime_sieves::{BIT64TOVAL240, MOD30_TO_MASK, WHEEL_2_3_5, WHEEL_2_3_5_7, sift};
use crate::utils::{
    FIArray::{FIArray, FIArrayU64},
    bit_array::BitArray,
    fenwick::{FenwickTree, FenwickTreeU32, FenwickTreeUsize},
    primes::{
        log_zeta::{log_zeta, log_zeta_reordered},
        primepi_approx::{Li, R},
    },
};
use fastdivide::DividerU64;
use itertools::Itertools;
const N: usize = 5e8 as usize;
// todo:
// calculate pi(n) the following way:
// using fenwick tree fiarrays, multiple zeta(s) by 1-p^-s for each p < n^1/3
// flatten fenwick tree into normal fiarray
// this computes phi(n, pi(n^1/3)) for each floor(n/d), and we know pi(n^1/3).
// all that is left is to evaluate P2(n, n^1/3), can be computed in O(sqrt(n))
// first step takes \Theta(n^2/3), so the algorithm in total runs in that timee complexity

// alternatively, try using ecnerwala's approach: sieve up to n^1/4, flatten, and compute P2 and P3

// repeated convolution of the prefix sum representation of u with mu_p for p below sqrt(n)
// I guess this is essentially legendre's formula for prime counting, implemented using bottom-up dp
// not efficient, lucy is essentially a smarter version of this, reducing the complexity from O(n/logn) to O(n^0.75/logn)
// can be optimised using fenwick trees and the sqrt trick
#[must_use]
pub fn legendre(x: usize) -> usize {
    let primes = sift(x.isqrt() as u64);
    let mut s = FIArray::unit(x);
    let keys = FIArray::keys(x).collect_vec();

    for &p in &primes {
        let p = p as usize;
        for (i, &v) in keys.iter().enumerate().rev() {
            if v < p {
                break;
            }
            s.arr[i] -= s[v / p];
        }
    }
    s[x] + primes.len() - 1
}

// O(n^(3/4)/log(n)) time, O(sqrt(n)) space prime counting function
// kinda simulates pritchard's wheel sieve
// 1e17: prime counting took 3071.1531213s: 2623557157654233
// 1e16: prime counting took 501.5590413s: 279238341033925
// 1e15: prime counting took 102.6765764s: 29844570422669
// 1e14: prime counting took 19.8688139s: 3204941750802
// 1e13: prime counting took 3.527102s: 346065536839
// 1e12: prime counting took 649.5083ms: 37607912018
// 1e11: prime counting took 117.2172ms: 4118054813
// 1e10: prime counting took 23.0519ms: 455052511
// 1e9: prime counting took 4.6073ms: 50847534
#[must_use]
pub fn lucy(x: usize) -> FIArray {
    const LUT: [usize; 30] = [
        0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 8,
    ];
    let mut s = FIArray::new(x);
    let keys = FIArray::keys(x).collect_vec();
    unsafe { core::hint::assert_unchecked(s.arr.len() == keys.len()) };
    for (i, v) in keys.iter().enumerate() {
        //s.arr[i] = (v + 1) >> 1;
        s.arr[i] = ((v / 30) << 3) + LUT[v % 30] - 1 + 3;
    }
    s.arr[0] = 0;
    s.arr[1] = 1;
    s.arr[2] = 2;
    s.arr[3] = 2;

    let mut pp = 25;
    unsafe { core::hint::assert_unchecked(s.arr.len() > x.isqrt()) };
    for p in (7..=x.isqrt()).step_by(2) {
        pp += (p << 2) - 4;
        let sp = s.arr[p - 2];
        if s.arr[p - 1] == sp {
            continue;
        }

        for (i, &v) in keys.iter().enumerate().rev() {
            if v < pp {
                break;
            }
            s.arr[i] -= s[v / p] - sp;
        }
    }
    s
}

#[must_use]
pub fn lucy_non_fiarray_alt(x: usize) -> usize {
    let isqrt = x.isqrt();
    let primes = sift(isqrt as u64);

    let mut small_s = vec![0; isqrt].into_boxed_slice();
    let mut large_s = vec![0; isqrt].into_boxed_slice();
    for v in 1..=isqrt {
        small_s[v - 1] = (v + 1) >> 1;
        large_s[v - 1] = (x / v + 1) >> 1;
    }
    small_s[0] = 0;
    for &p in &primes[1..] {
        let p = p as usize;
        let pp = p * p;
        let sp = small_s[p - 2];
        let mut dp = 0;
        let mut dpp = 0;
        for d in 1..=isqrt {
            dp += p;
            dpp += pp;

            if x < dpp {
                break;
            }
            large_s[d - 1] -= if x / dp <= isqrt {
                small_s[(x / dp) - 1]
            } else {
                large_s[dp - 1]
            } - sp;
        }
        for v in (1..=isqrt).rev() {
            if v < pp {
                break;
            }
            small_s[v - 1] -= small_s[(v / p) - 1] - sp;
        }
    }
    large_s[0]
}

#[must_use]
pub fn lucy_non_fiarray(x: usize) -> usize {
    const LUT: [usize; 30] = [
        0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 8,
    ];
    let isqrt = x.isqrt();

    let mut small_s = vec![0; isqrt].into_boxed_slice();
    let mut large_s = vec![0; isqrt - usize::from(isqrt == x / isqrt)].into_boxed_slice();

    for v in 1..isqrt {
        small_s[v - 1] = ((v / 30) << 3) + LUT[v % 30] - 1 + 3; //(v + 1) >> 1;
        large_s[v - 1] = (((x / v) / 30) << 3) + LUT[(x / v) % 30] - 1 + 3; //(x / v + 1) >> 1;
    }
    small_s[isqrt - 1] = ((isqrt / 30) << 3) + LUT[isqrt % 30] - 1 + 3; //(isqrt + 1) >> 1;
    if isqrt != x / isqrt {
        large_s[isqrt - 1] = (((x / isqrt) / 30) << 3) + LUT[(x / isqrt) % 30] - 1 + 3; //(x / isqrt + 1) >> 1;
    }

    small_s[0] = 0;
    small_s[1] = 1;
    small_s[2] = 2;
    small_s[3] = 2;

    let mut pp = 25;

    for p in (7..=x.isqrt()).step_by(2) {
        pp += (p << 2) - 4;
        let sp = small_s[p - 2];
        if small_s[p - 1] == sp {
            continue;
        }
        let mut dp = 0;
        let xp = x / p;
        for d in (1..isqrt).chain((isqrt != x / isqrt).then_some(isqrt)) {
            dp += p;

            if xp < dp {
                break;
            }
            large_s[d - 1] -= if xp / d <= isqrt {
                small_s[(xp / d) - 1]
            } else {
                large_s[dp - 1]
            } - sp;
        }
        for v in (pp..=isqrt).rev() {
            small_s[v - 1] -= small_s[(v / p) - 1] - sp;
        }
    }
    large_s[0]
}

#[must_use]
pub fn lucy_alt(x: usize) -> FIArray {
    let primes = sift(x.isqrt() as u64);
    let mut s = FIArray::new(x);
    let keys = FIArray::keys(x).collect_vec();

    unsafe { core::hint::assert_unchecked(s.arr.len() == keys.len()) };
    for (i, v) in keys.iter().enumerate() {
        s.arr[i] = (v + 1) >> 1;
    }
    s.arr[0] = 0;

    for &p in &primes[1..] {
        let p = p as usize;
        let sp = s.arr[p - 2];
        for (i, &v) in keys.iter().enumerate().rev() {
            if v < p * p {
                break;
            }
            s.arr[i] -= s[v / p] - sp;
        }
    }
    s
}

#[must_use]
pub fn lucy_alt_single(x: usize) -> usize {
    const LUT: [usize; 30] = [
        0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 8,
    ];
    let primes = sift(x.isqrt() as u64);
    let mut s = FIArray::new(x);
    let keys = FIArray::keys(x).collect_vec();

    unsafe { core::hint::assert_unchecked(s.arr.len() == keys.len()) };
    for (i, v) in keys.iter().enumerate() {
        //s.arr[i] = (v + 1) >> 1;
        s.arr[i] = ((v / 30) << 3) + LUT[v % 30] - 1 + 3;
    }
    s[1] = 0;
    s[2] = 1;
    s[3] = 2;
    s[4] = 2;
    s[5] = 3;

    let lim = primes.partition_point(|p| p.pow(3) <= x as u64);
    for &p in &primes[3..lim] {
        let p = p as usize;
        let sp = s.arr[p - 2];
        for (i, &v) in keys.iter().enumerate().rev() {
            if v < p * p {
                break;
            }
            s.arr[i] -= s[v / p] - sp;
        }
    }
    let mut res = s[x];
    for (i, &p) in primes[lim..].iter().enumerate() {
        let p = p as usize;
        res -= s[x / p] - lim - i;
    }
    res
}

// 1e17: 1506.2394159s
// 1e16: 300s
// 1e15: 58.8811143s
// 1e12: 345.7069ms
#[must_use]
pub fn prime_pi(x: usize) -> usize {
    const LUT: [usize; 30] = [
        0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 8,
    ];
    let isqrt = x.isqrt();
    let primes = sift(isqrt as u64); // technically unnecessary, algorithm already sieves

    let mut small_s = vec![0; isqrt].into_boxed_slice();
    let mut large_s = vec![0; isqrt].into_boxed_slice();
    for v in 1..=isqrt {
        //((v / 30) << 3) + LUT[v % 30] - 1 + 3
        small_s[v - 1] = ((v / 30) << 3) + LUT[v % 30] - 1 + 3; //(v + 1) >> 1;
        large_s[v - 1] = (((x / v) / 30) << 3) + LUT[(x / v) % 30] - 1 + 3; //(x / v + 1) >> 1;
    }
    small_s[0] = 0;
    small_s[1] = 1;
    small_s[2] = 2;
    small_s[3] = 2;
    let pi_4th_root = primes.partition_point(|p| p.pow(2) <= isqrt as u64);
    let pi_cbrt = primes.partition_point(|p| p.pow(3) <= x as u64);
    for &p in &primes[3..pi_4th_root] {
        let p = p as usize;
        let xp = x / p;
        let pp = p * p;
        let sp = small_s[p - 2];
        let mut d = p;
        let mut dp = pp;
        large_s[0] -= large_s[p - 1] - sp;
        if p == 7 {
            let mut incrs = WHEEL_2_3_5.into_iter().cycle();
            incrs.next();
            let incr = usize::from(unsafe { incrs.next().unwrap_unchecked() });
            d += incr;
            dp += incr * p;
            while d <= isqrt {
                //if small_s[d - 1] != small_s[d - 2] {
                large_s[d - 1] -= if xp / d <= isqrt {
                    small_s[(xp / d) - 1]
                } else {
                    large_s[dp - 1]
                } - sp;
                //}
                let incr = usize::from(unsafe { incrs.next().unwrap_unchecked() });
                d += incr;
                dp += incr * p;
            }
        } else {
            let mut incrs = WHEEL_2_3_5_7.into_iter().cycle();
            let mut m = (p % 210) - 1;
            while m != 0 {
                m -= usize::from(unsafe { incrs.next().unwrap_unchecked() });
            }
            let incr = usize::from(unsafe { incrs.next().unwrap_unchecked() });
            d += incr;
            dp += incr * p;
            //assert!(isqrt < xpp);
            while d <= isqrt {
                //if small_s[d - 1] != small_s[d - 2] {
                large_s[d - 1] -= if xp / d <= isqrt {
                    small_s[(xp / d) - 1]
                } else {
                    large_s[dp - 1]
                } - sp;
                //}
                let incr = usize::from(unsafe { incrs.next().unwrap_unchecked() });
                d += incr;
                dp += incr * p;
            }
        }
        for v in (pp..=isqrt).rev() {
            small_s[v - 1] -= small_s[(v / p) - 1] - sp;
        }
    }
    for (i, &p) in primes[pi_4th_root..pi_cbrt].iter().enumerate() {
        let p = p as usize;
        let xp = x / p;
        let xpp = xp / p;
        let sp = small_s[p - 2];
        large_s[0] -= large_s[p - 1] - sp;
        // each iteration does pi(x/(p*p)) - pi(p) work, x^1/4<=p<x^1/3
        for &d in &primes[pi_4th_root..][i + 1..] {
            let d = d as usize;
            if d > xpp {
                break;
            }
            large_s[d - 1] -= small_s[(xp / d) - 1] - sp;
        }
    }
    let mut res = large_s[0];

    // compute P2
    for &p in &primes[pi_cbrt..] {
        let p = p as usize;
        res -= large_s[p - 1] - small_s[p - 1] + 1;
    }
    res
}

// fucks up for small inputs
// starts being faster at around 10^11
// 1e17: 1231.2415746s
// 1e16: 230.4943082s
// todo: mix in wheel sieve, should noticeably improve first stage of the algorithm
// todo: store prime gaps instead of primes, only ever access them sequentially
// compress large_s,
#[must_use]
pub fn prime_pi_fenwick(x: usize) -> usize {
    const WHEEL: u32 =
        (1 << 1) | (1 << 7) | (1 << 11) | (1 << 13) | (1 << 17) | (1 << 19) | (1 << 23) | (1 << 29);
    const LUT: [usize; 30] = [
        0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 8,
    ];

    let isqrt = x.isqrt();
    let bitmap_size = (((isqrt / 30) >> 3) + 1) << 3;
    let mut sieve_raw = vec![0u8; bitmap_size].into_boxed_slice();
    sieve_raw[0] = 1;

    let bitmap = sieve_raw.as_mut_ptr();

    let mut large_s = vec![0; isqrt].into_boxed_slice();
    for v in 1..=isqrt {
        large_s[v - 1] = (((x / v) / 30) << 3) + LUT[(x / v) % 30] - 1 + 3;
    }
    let mut count = 0usize;
    let mut sieve =
        FenwickTreeUsize::new_with(isqrt, |i| usize::from((WHEEL >> ((i + 1) % 30)) & 1 == 1));
    sieve.dec(0);
    sieve.inc(1);
    sieve.inc(2);
    sieve.inc(4);

    let p = 7;
    let xp = x / p;
    let pp = 49;
    let sp = 3;
    let mut d = p;
    let mut dp = pp;
    large_s[0] -= large_s[p - 1] - sp;
    let mut incrs = WHEEL_2_3_5.into_iter().cycle();
    incrs.next();
    let unmark_incrs = incrs.clone();
    let incr = usize::from(unsafe { incrs.next().unwrap_unchecked() });
    d += incr;
    dp += incr * p;
    while d <= isqrt {
        large_s[d - 1] -= if xp / d <= isqrt {
            count += 1;
            sieve.sum((xp / d) - 1)
        } else {
            large_s[dp - 1]
        } - sp;
        let incr = usize::from(unsafe { incrs.next().unwrap_unchecked() });
        d += incr;
        dp += incr * p;
    }
    let precomp = [0, 0, 14, 0, 28, 0, 42];
    let mut multiple = pp;
    for incr in unmark_incrs {
        if multiple > isqrt {
            break;
        }
        if unsafe { *bitmap.add(multiple / 30) } & MOD30_TO_MASK[multiple % 30] == 0 {
            unsafe { *bitmap.add(multiple / 30) |= MOD30_TO_MASK[multiple % 30] };
            sieve.dec(multiple - 1);
        }
        multiple += precomp[incr as usize];
    }

    let mut pi = 4;
    let mut p = 11;
    let mut wheel_incr = WHEEL_2_3_5_7.into_iter().cycle();
    wheel_incr.next();
    while p * p <= isqrt {
        if unsafe { *bitmap.add(p / 30) } & MOD30_TO_MASK[p % 30] == 0 {
            let xp = x / p;
            let pp = p * p;
            let sp = pi;
            pi += 1;
            let mut d = p;
            let mut dp = pp;
            large_s[0] -= large_s[p - 1] - sp;
            let precomp = [
                0,
                0,
                p << 1,
                0,
                p << 2,
                0,
                (p << 2) + (p << 1),
                0,
                p << 3,
                0,
                (p << 1) + (p << 3),
            ];
            let mut incrs = wheel_incr.clone();

            let incr = usize::from(unsafe { incrs.next().unwrap_unchecked() });
            d += incr;
            dp += incr * p;
            //assert!(isqrt < xpp);
            while d <= isqrt {
                if unsafe { *bitmap.add(d / 30) } & MOD30_TO_MASK[d % 30] == 0 {
                    large_s[d - 1] -= if xp / d <= isqrt {
                        count += 1;
                        sieve.sum((xp / d) - 1)
                    } else {
                        large_s[dp - 1]
                    } - sp;
                }
                let incr = usize::from(unsafe { incrs.next().unwrap_unchecked() });
                d += incr;
                dp += precomp[incr];
            }
            let mut multiple = pp;
            for incr in wheel_incr.clone() {
                if multiple > isqrt {
                    break;
                }
                if unsafe { *bitmap.add(multiple / 30) } & MOD30_TO_MASK[multiple % 30] == 0 {
                    unsafe { *bitmap.add(multiple / 30) |= MOD30_TO_MASK[multiple % 30] };
                    sieve.dec(multiple - 1);
                }
                multiple += precomp[incr as usize];
            }
        }
        p += usize::from(unsafe { wheel_incr.next().unwrap_unchecked() });
    }
    // flatten tree here, no longer update small_s
    let small_s = sieve.flatten();
    dbg!(count, pi);
    let mut primes = Vec::with_capacity(small_s[isqrt - 1]);
    //primes.extend([2, 3, 5]);
    let bitmap64: *const u64 = bitmap.cast();
    let mut base = 0;
    for k in 0..bitmap_size / 8 - 1 {
        let mut bitset = unsafe { *bitmap64.add(k) };
        bitset = !bitset;
        while bitset != 0 {
            let r = bitset.trailing_zeros() as usize;
            primes.push(base + BIT64TOVAL240[r] as usize);
            bitset &= bitset - 1;
        }

        base += 240;
    }
    let mut bitset = unsafe { *bitmap64.add(bitmap_size / 8 - 1) };
    bitset = !bitset;
    while bitset != 0 {
        let r = bitset.trailing_zeros() as usize;
        let prime_cand = base + BIT64TOVAL240[r] as usize;
        if prime_cand > isqrt {
            break;
        }
        primes.push(prime_cand);
        bitset &= bitset - 1;
    }
    let pi_4th_root = pi - 3;
    let pi_cbrt = primes.partition_point(|p| p.pow(3) <= x);
    dbg!(pi_4th_root, pi_cbrt);
    for (i, &p) in primes[pi_4th_root..pi_cbrt].iter().enumerate() {
        let xp = x / p;
        let xpp = xp / p;
        let sp = pi;
        pi += 1;
        large_s[0] -= large_s[p - 1] - sp;
        // each iteration does pi(x/(p*p)) - pi(p) work, x^1/4<=p<x^1/3
        for &d in &primes[pi_4th_root..][i + 1..] {
            if d > xpp {
                break;
            }
            large_s[d - 1] -= small_s[(xp / d) - 1] - sp;
        }
    }
    let mut res = large_s[0];
    dbg!(res);
    // compute P2
    for &p in &primes[pi_cbrt..] {
        res -= large_s[p - 1] - pi;
        pi += 1;
    }
    res
}
const REMOVE_LESS: [u64; 240] = [
    0x0,
    0x1,
    0x1,
    0x1,
    0x1,
    0x1,
    0x1,
    0x3,
    0x3,
    0x3,
    0x3,
    0x7,
    0x7,
    0xf,
    0xf,
    0xf,
    0xf,
    0x1f,
    0x1f,
    0x3f,
    0x3f,
    0x3f,
    0x3f,
    0x7f,
    0x7f,
    0x7f,
    0x7f,
    0x7f,
    0x7f,
    0xff,
    0xff,
    0x1ff,
    0x1ff,
    0x1ff,
    0x1ff,
    0x1ff,
    0x1ff,
    0x3ff,
    0x3ff,
    0x3ff,
    0x3ff,
    0x7ff,
    0x7ff,
    0xfff,
    0xfff,
    0xfff,
    0xfff,
    0x1fff,
    0x1fff,
    0x3fff,
    0x3fff,
    0x3fff,
    0x3fff,
    0x7fff,
    0x7fff,
    0x7fff,
    0x7fff,
    0x7fff,
    0x7fff,
    0xffff,
    0xffff,
    0x1ffff,
    0x1ffff,
    0x1ffff,
    0x1ffff,
    0x1ffff,
    0x1ffff,
    0x3ffff,
    0x3ffff,
    0x3ffff,
    0x3ffff,
    0x7ffff,
    0x7ffff,
    0xfffff,
    0xfffff,
    0xfffff,
    0xfffff,
    0x1fffff,
    0x1fffff,
    0x3fffff,
    0x3fffff,
    0x3fffff,
    0x3fffff,
    0x7fffff,
    0x7fffff,
    0x7fffff,
    0x7fffff,
    0x7fffff,
    0x7fffff,
    0xffffff,
    0xffffff,
    0x1ffffff,
    0x1ffffff,
    0x1ffffff,
    0x1ffffff,
    0x1ffffff,
    0x1ffffff,
    0x3ffffff,
    0x3ffffff,
    0x3ffffff,
    0x3ffffff,
    0x7ffffff,
    0x7ffffff,
    0xfffffff,
    0xfffffff,
    0xfffffff,
    0xfffffff,
    0x1fffffff,
    0x1fffffff,
    0x3fffffff,
    0x3fffffff,
    0x3fffffff,
    0x3fffffff,
    0x7fffffff,
    0x7fffffff,
    0x7fffffff,
    0x7fffffff,
    0x7fffffff,
    0x7fffffff,
    0xffffffff,
    0xffffffff,
    0x1ffffffff,
    0x1ffffffff,
    0x1ffffffff,
    0x1ffffffff,
    0x1ffffffff,
    0x1ffffffff,
    0x3ffffffff,
    0x3ffffffff,
    0x3ffffffff,
    0x3ffffffff,
    0x7ffffffff,
    0x7ffffffff,
    0xfffffffff,
    0xfffffffff,
    0xfffffffff,
    0xfffffffff,
    0x1fffffffff,
    0x1fffffffff,
    0x3fffffffff,
    0x3fffffffff,
    0x3fffffffff,
    0x3fffffffff,
    0x7fffffffff,
    0x7fffffffff,
    0x7fffffffff,
    0x7fffffffff,
    0x7fffffffff,
    0x7fffffffff,
    0xffffffffff,
    0xffffffffff,
    0x1ffffffffff,
    0x1ffffffffff,
    0x1ffffffffff,
    0x1ffffffffff,
    0x1ffffffffff,
    0x1ffffffffff,
    0x3ffffffffff,
    0x3ffffffffff,
    0x3ffffffffff,
    0x3ffffffffff,
    0x7ffffffffff,
    0x7ffffffffff,
    0xfffffffffff,
    0xfffffffffff,
    0xfffffffffff,
    0xfffffffffff,
    0x1fffffffffff,
    0x1fffffffffff,
    0x3fffffffffff,
    0x3fffffffffff,
    0x3fffffffffff,
    0x3fffffffffff,
    0x7fffffffffff,
    0x7fffffffffff,
    0x7fffffffffff,
    0x7fffffffffff,
    0x7fffffffffff,
    0x7fffffffffff,
    0xffffffffffff,
    0xffffffffffff,
    0x1ffffffffffff,
    0x1ffffffffffff,
    0x1ffffffffffff,
    0x1ffffffffffff,
    0x1ffffffffffff,
    0x1ffffffffffff,
    0x3ffffffffffff,
    0x3ffffffffffff,
    0x3ffffffffffff,
    0x3ffffffffffff,
    0x7ffffffffffff,
    0x7ffffffffffff,
    0xfffffffffffff,
    0xfffffffffffff,
    0xfffffffffffff,
    0xfffffffffffff,
    0x1fffffffffffff,
    0x1fffffffffffff,
    0x3fffffffffffff,
    0x3fffffffffffff,
    0x3fffffffffffff,
    0x3fffffffffffff,
    0x7fffffffffffff,
    0x7fffffffffffff,
    0x7fffffffffffff,
    0x7fffffffffffff,
    0x7fffffffffffff,
    0x7fffffffffffff,
    0xffffffffffffff,
    0xffffffffffffff,
    0x1ffffffffffffff,
    0x1ffffffffffffff,
    0x1ffffffffffffff,
    0x1ffffffffffffff,
    0x1ffffffffffffff,
    0x1ffffffffffffff,
    0x3ffffffffffffff,
    0x3ffffffffffffff,
    0x3ffffffffffffff,
    0x3ffffffffffffff,
    0x7ffffffffffffff,
    0x7ffffffffffffff,
    0xfffffffffffffff,
    0xfffffffffffffff,
    0xfffffffffffffff,
    0xfffffffffffffff,
    0x1fffffffffffffff,
    0x1fffffffffffffff,
    0x3fffffffffffffff,
    0x3fffffffffffffff,
    0x3fffffffffffffff,
    0x3fffffffffffffff,
    0x7fffffffffffffff,
    0x7fffffffffffffff,
    0x7fffffffffffffff,
    0x7fffffffffffffff,
    0x7fffffffffffffff,
    0x7fffffffffffffff,
    0xffffffffffffffff,
];
// 1e18: res = 24739954287740860, took 7059.9063396s
// 1e17: 1185.792844s
// 1e16: res = 279238341033925, took 214.8473257s
// TODO: utilise sqrt trick, few possible values for xp/d
#[must_use]
pub fn prime_pi_fenwick_2(x: usize) -> usize {
    const LUT: [usize; 30] = [
        0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 8,
    ];
    let isqrt = x.isqrt();
    let bitmap_size = ((isqrt / 30) >> 3) + 1;
    let mut sieve_raw = vec![0u64; bitmap_size].into_boxed_slice();
    sieve_raw[0] = 1;

    let bitmap: *mut u8 = sieve_raw.as_mut_ptr().cast();

    let mut large_s = vec![0; isqrt].into_boxed_slice();
    for v in 1..=isqrt {
        large_s[v - 1] = (((x / v) / 30) << 3) + LUT[(x / v) % 30] - 1 + 3;
    }
    let mut count = 0usize;
    let mut sieve = FenwickTreeU32::new(bitmap_size, 64);
    sieve.dec(0); // remove 1
    // add 2,3,5 as necessary later
    let p = 7;
    let xp = x / p;
    let pp = 49;
    let sp = 3;
    let mut d = p;
    let mut dp = pp;
    large_s[0] -= large_s[p - 1] - sp;
    let mut incrs = WHEEL_2_3_5.into_iter().cycle();
    incrs.next();
    let unmark_incrs = incrs.clone();
    let incr = usize::from(unsafe { incrs.next().unwrap_unchecked() });
    d += incr;
    dp += incr * p;
    while d <= isqrt {
        let xpd = xp / d;
        large_s[d - 1] -= if xpd <= isqrt {
            count += 1;
            let mut ret = sieve.sum(xpd / 240) as usize
                - (sieve_raw[xpd / 240] | REMOVE_LESS[xpd % 240]).count_zeros() as usize;
            ret += usize::from(xpd > 1) + usize::from(xpd > 2) + usize::from(xpd > 4);
            ret
        } else {
            large_s[dp - 1]
        } - sp;
        let incr = usize::from(unsafe { incrs.next().unwrap_unchecked() });
        d += incr;
        dp += incr * p;
    }
    let precomp = [0, 0, 14, 0, 28, 0, 42];
    let mut multiple = pp;
    for incr in unmark_incrs {
        if multiple > isqrt {
            break;
        }
        if unsafe { *bitmap.add(multiple / 30) } & MOD30_TO_MASK[multiple % 30] == 0 {
            unsafe { *bitmap.add(multiple / 30) |= MOD30_TO_MASK[multiple % 30] };
            sieve.dec(multiple / 240);
        }
        multiple += precomp[incr as usize];
    }

    let mut pi = 4;
    let mut p = 11;
    let mut wheel_incr = WHEEL_2_3_5_7.into_iter().cycle();
    wheel_incr.next();
    while p * p <= isqrt {
        if unsafe { *bitmap.add(p / 30) } & MOD30_TO_MASK[p % 30] == 0 {
            //dbg!(p,pi);
            let xp = x / p;
            let pp = p * p;
            let sp = pi;
            pi += 1;
            let mut d = p;
            let mut dp = pp;
            large_s[0] -= large_s[p - 1] - sp;
            let precomp = [
                0,
                0,
                p << 1,
                0,
                p << 2,
                0,
                (p << 2) + (p << 1),
                0,
                p << 3,
                0,
                (p << 1) + (p << 3),
            ];
            let mut incrs = wheel_incr.clone();

            let incr = usize::from(unsafe { incrs.next().unwrap_unchecked() });
            d += incr;
            dp += incr * p;
            //assert!(isqrt < xpp);
            while d <= isqrt {
                if unsafe { *bitmap.add(d / 30) } & MOD30_TO_MASK[d % 30] == 0 {
                    let xpd = xp / d;
                    large_s[d - 1] -= if xpd <= isqrt {
                        count += 1;
                        let mut ret = sieve.sum(xpd / 240) as usize
                            - (sieve_raw[xpd / 240] | REMOVE_LESS[xpd % 240]).count_zeros()
                                as usize;
                        ret += usize::from(xpd > 1) + usize::from(xpd > 2) + usize::from(xpd > 4);
                        ret
                    } else {
                        large_s[dp - 1]
                    } - sp;
                }
                let incr = usize::from(unsafe { incrs.next().unwrap_unchecked() });
                d += incr;
                dp += precomp[incr];
            }
            let mut multiple = pp;
            for incr in wheel_incr.clone() {
                if multiple > isqrt {
                    break;
                }
                if unsafe { *bitmap.add(multiple / 30) } & MOD30_TO_MASK[multiple % 30] == 0 {
                    unsafe { *bitmap.add(multiple / 30) |= MOD30_TO_MASK[multiple % 30] };
                    sieve.dec(multiple / 240);
                }
                multiple += precomp[incr as usize];
            }
        }
        p += usize::from(unsafe { wheel_incr.next().unwrap_unchecked() });
    }
    // flatten tree here, no longer update small_s
    let small_s = sieve.flatten();
    dbg!(count, pi);
    let mut primes = Vec::with_capacity(small_s[isqrt / 240] as usize);
    //primes.extend([2, 3, 5]);
    let bitmap64: *const u64 = bitmap.cast();
    let mut base = 0;
    for k in 0..bitmap_size - 1 {
        let mut bitset = unsafe { *bitmap64.add(k) };
        bitset = !bitset;
        while bitset != 0 {
            let r = bitset.trailing_zeros() as usize;
            primes.push((base + BIT64TOVAL240[r] as usize) as u32);
            bitset &= bitset - 1;
        }

        base += 240;
    }
    let mut bitset = unsafe { *bitmap64.add(bitmap_size - 1) };
    bitset = !bitset;
    while bitset != 0 {
        let r = bitset.trailing_zeros() as usize;
        let prime_cand = base + BIT64TOVAL240[r] as usize;
        if prime_cand > isqrt {
            break;
        }
        primes.push(prime_cand as u32);
        bitset &= bitset - 1;
    }
    let pi_4th_root = pi - 3;
    let pi_cbrt = primes.partition_point(|&p| (p as usize).pow(3) <= x);
    dbg!(pi_4th_root, pi_cbrt);
    for (i, &p) in primes[pi_4th_root..pi_cbrt].iter().enumerate() {
        let p = p as usize;
        let xp = x / p;
        let xpp = xp / p;
        let sp = pi;
        pi += 1;
        large_s[0] -= large_s[p - 1] - sp;
        // each iteration does pi(x/(p*p)) - pi(p) work, x^1/4<=p<x^1/3
        for &d in &primes[pi_4th_root..][i + 1..] {
            let d = d as usize;
            if d > xpp {
                break;
            }
            let xpd = xp / d;
            let ret = 3 + small_s[xpd / 240] as usize
                - (sieve_raw[xpd / 240] | REMOVE_LESS[xpd % 240]).count_zeros() as usize;
            large_s[d - 1] -= ret - sp;
        }
    }
    let mut res = large_s[0];
    dbg!(res);
    // compute P2
    for &p in &primes[pi_cbrt..] {
        let p = p as usize;
        res -= large_s[p - 1] - pi;
        pi += 1;
    }
    res
}

#[must_use]
pub fn lucy_wheel(x: usize) -> FIArray {
    const LUT: [usize; 30] = [
        0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 8,
    ];
    let mut s = FIArray::new(x);
    let keys = FIArray::keys(x).collect_vec();

    unsafe { core::hint::assert_unchecked(s.arr.len() == keys.len()) };
    for (i, v) in keys.iter().enumerate() {
        //s.arr[i] = (v + 1) >> 1;
        s.arr[i] = ((v / 30) << 3) + LUT[v % 30] - 1 + 3;
    }
    s[1] = 0;
    s[2] = 1;
    s[3] = 2;
    s[4] = 2;
    s[5] = 3;
    unsafe { core::hint::assert_unchecked(s.arr.len() > x.isqrt()) };
    let lim = x.isqrt();
    assert!(lim >= 5);
    let mut p = 1;
    let mut incrs = WHEEL_2_3_5.into_iter().cycle();
    loop {
        let incr = usize::from(unsafe { incrs.next().unwrap_unchecked() });
        p += incr;
        if p > lim {
            break;
        }

        if s.arr[p - 1] == s.arr[p - 2] {
            continue;
        }

        let sp = s.arr[p - 2];

        for (i, &v) in keys.iter().enumerate().rev() {
            if v < p * p {
                break;
            }
            s.arr[i] -= s[v / p] - sp;
        }
    }
    s
}
#[must_use]
pub fn lucy_wheel210(x: usize) -> FIArray {
    let mut s = FIArray::new(x);
    let keys = FIArray::keys(x).collect_vec();

    unsafe { core::hint::assert_unchecked(s.arr.len() == keys.len()) };
    for (i, v) in keys.iter().enumerate() {
        s.arr[i] = v - 1;
    }
    unsafe { core::hint::assert_unchecked(s.arr.len() > x.isqrt()) };
    let lim = x.isqrt();
    assert!(lim >= 7);
    for p in [2, 3, 5, 7] {
        let sp = s.arr[p - 2];

        for (i, &v) in keys.iter().enumerate().rev() {
            if v < p * p {
                break;
            }
            s.arr[i] -= s[v / p] - sp;
        }
    }
    let mut p = 1;
    let mut incrs = WHEEL_2_3_5_7.into_iter().cycle();
    loop {
        let incr = usize::from(unsafe { incrs.next().unwrap_unchecked() });
        p += incr;
        if p > lim {
            break;
        }

        if s.arr[p - 1] == s.arr[p - 2] {
            continue;
        }

        let sp = s.arr[p - 2];

        for (i, &v) in keys.iter().enumerate().rev() {
            if v < p * p {
                break;
            }
            s.arr[i] -= s[v / p] - sp;
        }
    }
    s
}

#[must_use]
pub fn lucy_fastdivide(x: u64) -> FIArrayU64 {
    let mut s = FIArrayU64::new(x);
    let keys = FIArrayU64::keys(x).collect_vec();

    unsafe { core::hint::assert_unchecked(s.arr.len() == keys.len()) };
    for (i, v) in keys.iter().enumerate() {
        s.arr[i] = (v + 1) >> 1;
    }
    s.arr[0] = 0;
    s.arr[2] = 2;
    let mut pp = 1;
    unsafe { core::hint::assert_unchecked(s.arr.len() as u64 > x.isqrt()) };
    for p in (3..=x.isqrt()).step_by(2) {
        pp += (p << 2) - 4;
        let sp = s.arr[p as usize - 2];
        if s.arr[p as usize - 1] == sp {
            continue;
        }

        let pdiv = DividerU64::divide_by(p);
        for (i, &v) in keys.iter().enumerate().rev() {
            if v < pp {
                break;
            }
            s.arr[i] -= s[v / pdiv] - sp;
        }
    }
    s
}
#[must_use]
pub fn lucy_fastdivide_alt(x: u64) -> FIArrayU64 {
    let primes = sift(x.isqrt());

    let mut s = FIArrayU64::new(x);
    let keys = FIArrayU64::keys(x).collect_vec();

    unsafe { core::hint::assert_unchecked(s.arr.len() == keys.len()) };
    for (i, v) in keys.iter().enumerate() {
        s.arr[i] = (v + 1) >> 1;
    }
    s.arr[0] = 0;
    s.arr[2] = 2;
    unsafe { core::hint::assert_unchecked(s.arr.len() as u64 > x.isqrt()) };
    for &p in &primes[1..] {
        let sp = s.arr[p as usize - 2];

        let pdiv = DividerU64::divide_by(p);
        for (i, &v) in keys.iter().enumerate().rev() {
            if v < p * p {
                break;
            }
            s.arr[i] -= s[v / pdiv] - sp;
        }
    }
    s
}
#[must_use]
pub fn lucy_fastdivide_wheel(x: u64) -> FIArrayU64 {
    let mut s = FIArrayU64::new(x);
    let keys = FIArrayU64::keys(x).collect_vec();

    unsafe { core::hint::assert_unchecked(s.arr.len() == keys.len()) };
    for (i, v) in keys.iter().enumerate() {
        s.arr[i] = (v + 1) >> 1;
    }
    s.arr[0] = 0;
    s.arr[2] = 2;
    unsafe { core::hint::assert_unchecked(s.arr.len() as u64 > x.isqrt()) };
    let lim = x.isqrt();
    assert!(lim >= 5);
    for p in [3, 5] {
        let sp = s.arr[p as usize - 2];

        let pdiv = DividerU64::divide_by(p);
        for (i, &v) in keys.iter().enumerate().rev() {
            if v < p * p {
                break;
            }
            s.arr[i] -= s[v / pdiv] - sp;
        }
    }
    let mut p = 1;
    let mut incrs = WHEEL_2_3_5.into_iter().cycle();
    loop {
        p += u64::from(unsafe { incrs.next().unwrap_unchecked() });
        if p > lim {
            break;
        }

        if s.arr[p as usize - 1] == s.arr[p as usize - 2] {
            continue;
        }

        let sp = s.arr[p as usize - 2];

        let pdiv = DividerU64::divide_by(p);
        for (i, &v) in keys.iter().enumerate().rev() {
            if v < p * p {
                break;
            }
            s.arr[i] -= s[v / pdiv] - sp;
        }
    }
    s
}
#[must_use]
pub fn lucy_fastdivide_wheel210(x: u64) -> FIArrayU64 {
    let mut s = FIArrayU64::new(x);
    let keys = FIArrayU64::keys(x).collect_vec();

    unsafe { core::hint::assert_unchecked(s.arr.len() == keys.len()) };
    for (i, v) in keys.iter().enumerate() {
        s.arr[i] = (v + 1) >> 1;
    }
    s.arr[0] = 0;
    s.arr[2] = 2;
    unsafe { core::hint::assert_unchecked(s.arr.len() as u64 > x.isqrt()) };
    let lim = x.isqrt();
    assert!(lim >= 7);
    for p in [3, 5, 7] {
        let sp = s.arr[p as usize - 2];

        let pdiv = DividerU64::divide_by(p);
        for (i, &v) in keys.iter().enumerate().rev() {
            if v < p * p {
                break;
            }
            s.arr[i] -= s[v / pdiv] - sp;
        }
    }
    let mut p = 1;
    let mut incrs = WHEEL_2_3_5_7.into_iter().cycle();
    loop {
        p += u64::from(unsafe { incrs.next().unwrap_unchecked() });
        if p > lim {
            break;
        }
        if s.arr[p as usize - 1] == s.arr[p as usize - 2] {
            continue;
        }
        let sp = s.arr[p as usize - 2];

        let pdiv = DividerU64::divide_by(p);
        for (i, &v) in keys.iter().enumerate().rev() {
            if v < p * p {
                break;
            }
            s.arr[i] -= s[v / pdiv] - sp;
        }
    }
    s
}

#[must_use]
pub fn lucy_strengthreduce(x: usize) -> FIArray {
    let mut s = FIArray::new(x);
    let keys = FIArray::keys(x).collect_vec();

    unsafe { core::hint::assert_unchecked(s.arr.len() == keys.len()) };
    for (i, v) in keys.iter().enumerate() {
        s.arr[i] = (v + 1) >> 1;
    }

    s.arr[0] = 0;

    let mut pp = 1;
    unsafe { core::hint::assert_unchecked(s.arr.len() > x.isqrt()) };
    for p in (3..=x.isqrt()).step_by(2) {
        pp += (p << 2) - 4;
        let sp = s.arr[p - 2];
        if s.arr[p - 1] == sp {
            continue;
        }

        let pdiv = strength_reduce::StrengthReducedUsize::new(p);
        for (i, &v) in keys.iter().enumerate().rev() {
            if v < pp {
                break;
            }
            s.arr[i] -= s[v / pdiv] - sp;
        }
    }
    s
}
#[must_use]
pub fn lucy_strengthreduce_alt(x: usize) -> FIArray {
    let primes = sift(x.isqrt() as u64);
    let mut s = FIArray::new(x);
    let keys = FIArray::keys(x).collect_vec();

    unsafe { core::hint::assert_unchecked(s.arr.len() == keys.len()) };
    for (i, v) in keys.iter().enumerate() {
        s.arr[i] = (v + 1) >> 1;
    }
    s.arr[0] = 0;
    s.arr[2] = 2;
    unsafe { core::hint::assert_unchecked(s.arr.len() > x.isqrt()) };
    for &p in &primes[1..] {
        let p = p as usize;
        let sp = s.arr[p - 2];

        let pdiv = strength_reduce::StrengthReducedUsize::new(p);
        for (i, &v) in keys.iter().enumerate().rev() {
            if v < p * p {
                break;
            }
            s.arr[i] -= s[v / pdiv] - sp;
        }
    }
    s
}

// easier to understand, but completely inferior due to many more integer divisions: never use
#[must_use]
pub fn lucy_dumber(x: usize) -> FIArray {
    let mut s = FIArray::new(x);
    for v in FIArray::keys(x) {
        s[v] = v - 1;
    }
    for p in 2..=x.isqrt() {
        if s[p] == s[p - 1] {
            continue;
        }
        for v in FIArray::keys(x).rev() {
            if v < p * p {
                break;
            }
            s[v] -= s[v / p] - s[p - 1];
        }
    }
    s
}

pub fn main() {
    println!("{N}");
    println!("logarithmic integral:");
    let start = Instant::now();
    let count = Li(N as _);
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    println!("Riemann R:");
    let start = Instant::now();
    let count = R(N as _);
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    /*  println!("lucy fenwick:");
    let start = Instant::now();
    let count = lucy_fenwick_trick(N as _);
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    let start = Instant::now();
    let count = lucy_fenwick(N as _)[N as _];
    let end = start.elapsed();
    println!("res = {count}, took {end:?}"); */

    /* println!("legendre:");
    let start = Instant::now();
    let count = legendre(N as _);
    let end = start.elapsed();
    println!("res = {count}, took {end:?}"); */

    println!("standard-ish lucy");
    let start = Instant::now();
    let count = prime_pi_fenwick_2(N as _);
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    let start = Instant::now();
    let count = prime_pi_fenwick(N as _);
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    let start = Instant::now();
    let count = prime_pi(N as _);
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    let start = Instant::now();
    let count = lucy_alt_single(N as _);
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    let start = Instant::now();
    let count = lucy(N as _)[N as _];
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    let start = Instant::now();
    let count = lucy_non_fiarray(N as _);
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    println!("prime counting using the logarithm of the zeta function:");
    let start = Instant::now();
    let count = log_zeta_reordered(N as _)[N as _]; // n^(2/3)
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    let start = Instant::now();
    let count = log_zeta(N as _)[N as _]; // n^(2/3)
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");
}

//never faster for me, though asymptotically better: O(x^(2/3) (logx)(loglogx)^(1/3)) vs O(x^(3/4) / log(x))
#[must_use]
pub fn lucy_fenwick(x: usize) -> FIArray {
    let mut s = FIArray::new(x);
    let sqrtx = x.isqrt();
    let len = s.arr.len();
    let y = (2e9 as usize)
        .min((((x as f64) / (x as f64).ln().ln()).powf(2. / 3.)).round() as usize >> 3)
        .max(sqrtx + 1);
    let mut sieve_raw = BitArray::zeroed(y + 1);
    let mut sieve = FenwickTree::new(y + 1, 1);
    sieve.add(1, -1);
    sieve.add(0, -1);

    for (i, v) in FIArray::keys(x).enumerate() {
        s.arr[i] = v - 1;
    }

    for p in 2..=sqrtx {
        if !sieve_raw.get(p) {
            let sp = sieve.sum(p - 1) as usize;
            let lim = (x / y).min(x / (p * p));
            for i in 1..=lim {
                let xip = x / (i * p);
                s.arr[len - i] -= if xip <= y {
                    sieve.sum(xip) as usize
                } else {
                    s[xip]
                } - sp;
            }
            for j in (p * p..=y).step_by(p) {
                if !sieve_raw.get(j) {
                    sieve_raw.set(j);
                    sieve.add(j, -1);
                }
            }
        }
    }
    for (i, v) in FIArray::keys(x).take_while(|&v| v <= y).enumerate() {
        s.arr[i] = sieve.sum(v) as usize;
    }
    /* dbg!((count_sum, count_add));
    dbg!(y * (y as f64).ln().ln() as usize); */
    s
}

#[must_use]
pub fn lucy_fenwick_trick(x: usize) -> usize {
    let mut s = FIArray::new(x);
    let sqrtx = x.isqrt();
    let len = s.arr.len();
    let y = (2e9 as usize)
        .min((((x as f64) / (x as f64).ln().ln()).powf(2. / 3.)).round() as usize >> 4)
        .max(sqrtx + 1);
    let mut sieve_raw = BitArray::zeroed(y + 1);
    let mut sieve = FenwickTree::new(y + 1, 1);
    sieve.add(1, -1);
    sieve.add(0, -1);

    for (i, v) in FIArray::keys(x).enumerate() {
        s.arr[i] = v - 1;
    }

    for p in 2..=sqrtx {
        if sieve_raw.get(p) {
            continue;
        }
        let sp = sieve.sum(p - 1) as usize;

        let lim = (x / y).min(x / (p * p));
        let xp = x / p;
        s.arr[len - 1] -= if xp <= y {
            sieve.sum(xp) as usize
        } else {
            s[xp]
        } - sp;
        for i in p..=lim {
            if sieve_raw.get(i) {
                continue;
            }
            let xip = xp / i;
            s.arr[len - i] -= if xip <= y {
                sieve.sum(xip) as usize
            } else {
                s[xip]
            } - sp;
        }
        for j in (p * p..=y).step_by(p) {
            if !sieve_raw.get(j) {
                sieve_raw.set(j);
                sieve.add(j, -1);
            }
        }
    }
    /* dbg!((count_sum, count_add));
       dbg!(y * (y as f64).ln().ln() as usize);
    */
    s[x]
}

//TODO: finish
/*#[must_use]
pub fn prime_pi_fenwick_3(x: usize) -> usize {
    const LUT: [usize; 30] = [
        0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 8,
    ];
    let isqrt = x.isqrt();
    let y = (1e9 as usize)
        .min((((x as f64) / (x as f64).ln().ln()).powf(2. / 3.)).round() as usize)
        .max(isqrt);
    let bitmap_size = ((y / 30) >> 3) + 1;
    let mut sieve_raw = vec![0u64; bitmap_size].into_boxed_slice();
    sieve_raw[0] = 1;

    let bitmap: *mut u8 = sieve_raw.as_mut_ptr().cast();

    let mut large_s = vec![0; isqrt].into_boxed_slice();
    for v in 1..=isqrt {
        large_s[v - 1] = (((x / v) / 30) << 3) + LUT[(x / v) % 30] - 1 + 3;
    }
    let mut count = 0usize;
    let mut sieve = FenwickTreeU32::new(bitmap_size, 64);
    sieve.dec(0); // remove 1
    // add 2,3,5 as necessary later
    let p = 7;
    let xp = x / p;
    let pp = 49;
    let sp = 3;
    let mut d = p;
    let mut dp = pp;
    large_s[0] -= large_s[p - 1] - sp;
    let mut incrs = WHEEL_2_3_5.into_iter().cycle();
    incrs.next();
    let unmark_incrs = incrs.clone();
    let incr = usize::from(unsafe { incrs.next().unwrap_unchecked() });
    d += incr;
    dp += incr * p;
    while d <= isqrt {
        let xpd = xp / d;
        large_s[d - 1] -= if xpd <= isqrt {
            count += 1;
            let mut ret = sieve.sum(xpd / 240) as usize
                - (sieve_raw[xpd / 240] | REMOVE_LESS[xpd % 240]).count_zeros() as usize;
            ret += usize::from(xpd > 1) + usize::from(xpd > 2) + usize::from(xpd > 4);
            ret
        } else {
            large_s[dp - 1]
        } - sp;
        let incr = usize::from(unsafe { incrs.next().unwrap_unchecked() });
        d += incr;
        dp += incr * p;
    }
    let precomp = [0, 0, 14, 0, 28, 0, 42];
    let mut multiple = pp;
    for incr in unmark_incrs {
        if multiple > y {
            break;
        }
        if unsafe { *bitmap.add(multiple / 30) } & MOD30_TO_MASK[multiple % 30] == 0 {
            unsafe { *bitmap.add(multiple / 30) |= MOD30_TO_MASK[multiple % 30] };
            sieve.dec(multiple / 240);
        }
        multiple += precomp[incr as usize];
    }

    let mut pi = 4;
    let mut p = 11;
    let mut wheel_incr = WHEEL_2_3_5_7.into_iter().cycle();
    wheel_incr.next();
    while p * p <= isqrt {
        if unsafe { *bitmap.add(p / 30) } & MOD30_TO_MASK[p % 30] == 0 {
            //dbg!(p,pi);
            let xp = x / p;
            let pp = p * p;
            let sp = pi;
            pi += 1;
            let mut d = p;
            let mut dp = pp;
            large_s[0] -= large_s[p - 1] - sp;
            let precomp = [
                0,
                0,
                p << 1,
                0,
                p << 2,
                0,
                (p << 2) + (p << 1),
                0,
                p << 3,
                0,
                (p << 1) + (p << 3),
            ];
            let mut incrs = wheel_incr.clone();

            let incr = usize::from(unsafe { incrs.next().unwrap_unchecked() });
            d += incr;
            dp += incr * p;
            //assert!(isqrt < xpp);
            while d <= isqrt {
                if unsafe { *bitmap.add(d / 30) } & MOD30_TO_MASK[d % 30] == 0 {
                    let xpd = xp / d;
                    large_s[d - 1] -= if xpd <= isqrt {
                        count += 1;
                        let mut ret = sieve.sum(xpd / 240) as usize
                            - (sieve_raw[xpd / 240] | REMOVE_LESS[xpd % 240]).count_zeros()
                                as usize;
                        ret += usize::from(xpd > 1) + usize::from(xpd > 2) + usize::from(xpd > 4);
                        ret
                    } else {
                        large_s[dp - 1]
                    } - sp;
                }
                let incr = usize::from(unsafe { incrs.next().unwrap_unchecked() });
                d += incr;
                dp += precomp[incr];
            }
            let mut multiple = pp;
            for incr in wheel_incr.clone() {
                if multiple > y {
                    break;
                }
                if unsafe { *bitmap.add(multiple / 30) } & MOD30_TO_MASK[multiple % 30] == 0 {
                    unsafe { *bitmap.add(multiple / 30) |= MOD30_TO_MASK[multiple % 30] };
                    sieve.dec(multiple / 240);
                }
                multiple += precomp[incr as usize];
            }
        }
        p += usize::from(unsafe { wheel_incr.next().unwrap_unchecked() });
    }
    // flatten tree here, no longer update small_s
    let small_s = sieve.flatten();
    dbg!(count, pi);
    let mut primes = Vec::with_capacity(small_s[isqrt / 240] as usize);
    //primes.extend([2, 3, 5]);
    let bitmap64: *const u64 = bitmap.cast();
    let mut base = 0;
    for k in 0..bitmap_size - 1 {
        let mut bitset = unsafe { *bitmap64.add(k) };
        bitset = !bitset;
        while bitset != 0 {
            let r = bitset.trailing_zeros() as usize;
            primes.push((base + BIT64TOVAL240[r] as usize) as u32);
            bitset &= bitset - 1;
        }

        base += 240;
    }
    let mut bitset = unsafe { *bitmap64.add(bitmap_size - 1) };
    bitset = !bitset;
    while bitset != 0 {
        let r = bitset.trailing_zeros() as usize;
        let prime_cand = base + BIT64TOVAL240[r] as usize;
        if prime_cand > isqrt {
            break;
        }
        primes.push(prime_cand as u32);
        bitset &= bitset - 1;
    }
    let pi_4th_root = pi - 3;
    let pi_cbrt = primes.partition_point(|&p| (p as usize).pow(3) <= x);
    dbg!(pi_4th_root, pi_cbrt);
    for (i, &p) in primes[pi_4th_root..pi_cbrt].iter().enumerate() {
        let p = p as usize;
        let xp = x / p;
        let xpp = xp / p;
        let sp = pi;
        pi += 1;
        large_s[0] -= large_s[p - 1] - sp;
        // each iteration does pi(x/(p*p)) - pi(p) work, x^1/4<=p<x^1/3
        for &d in &primes[pi_4th_root..][i + 1..] {
            let d = d as usize;
            if d > xpp {
                break;
            }
            let xpd = xp / d;
            let ret = 3 + small_s[xpd / 240] as usize
                - (sieve_raw[xpd / 240] | REMOVE_LESS[xpd % 240]).count_zeros() as usize;
            large_s[d - 1] -= ret - sp;
        }
    }
    let mut res = large_s[0];
    dbg!(res);
    // compute P2
    for &p in &primes[pi_cbrt..] {
        let p = p as usize;
        res -= large_s[p - 1] - pi;
        pi += 1;
    }
    res
}
*/
