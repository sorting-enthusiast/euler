use std::time::Instant;

use fastdivide::DividerU64;
use itertools::Itertools;

use crate::utils::{
    FIArray::{FIArray, FIArrayU64, FIArrayU128},
    bit_array::BitArray,
    fenwick::FenwickTree,
    multiplicative_function_summation::mobius_sieve,
    prime_sieves::{BIT64TOVAL240, WHEEL_2_3_5, WHEEL_2_3_5_7, sift},
};
const N: usize = 1e16 as usize;

// repeated convolution of the prefix sum representation of u with mu_p for p below sqrt(n)
// I guess this is essentially legendre's formula for prime counting, implemented using bottom-up dp
// not efficient, lucy is essentially a smarter version of this, reducing the complexity from O(n/logn) to O(n^0.75/logn)
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
                if small_s[d - 1] != small_s[d - 2] {
                    large_s[d - 1] -= if xp / d <= isqrt {
                        small_s[(xp / d) - 1]
                    } else {
                        large_s[dp - 1]
                    } - sp;
                }
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
                if small_s[d - 1] != small_s[d - 2] {
                    large_s[d - 1] -= if xp / d <= isqrt {
                        small_s[(xp / d) - 1]
                    } else {
                        large_s[dp - 1]
                    } - sp;
                }
                let incr = usize::from(unsafe { incrs.next().unwrap_unchecked() });
                d += incr;
                dp += incr * p;
            }
        }
        // todo: replace with fenwick tree
        for v in (pp..=isqrt).rev() {
            small_s[v - 1] -= small_s[(v / p) - 1] - sp;
        }
    }
    // flatten tree here, no longer update small_s
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
// 1e16: 236.6590369s
#[must_use]
pub fn prime_pi_fenwick(x: usize) -> usize {
    const WHEEL: u32 =
        (1 << 1) | (1 << 7) | (1 << 11) | (1 << 13) | (1 << 17) | (1 << 19) | (1 << 23) | (1 << 29);
    const LUT: [usize; 30] = [
        0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 8,
    ];
    const MOD30_TO_MASK: [u8; 30] = [
        0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 4, 4, 8, 8, 8, 8, 16, 16, 32, 32, 32, 32, 64, 64, 64, 64,
        64, 64, 128,
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
    let mut count = 0;
    let mut sieve = FenwickTree::new_with(isqrt + 1, |i| i64::from((WHEEL >> (i % 30)) & 1 == 1));
    sieve.add(1, -1);
    sieve.add(2, 1);
    sieve.add(3, 1);
    sieve.add(5, 1);

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
            sieve.sum(xp / d) as usize
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
            sieve.add(multiple, -1);
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
                        sieve.sum(xp / d) as usize
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
                    sieve.add(multiple, -1);
                }
                multiple += precomp[incr as usize];
            }
        }
        p += usize::from(unsafe { wheel_incr.next().unwrap_unchecked() });
    }
    // flatten tree here, no longer update small_s
    /* for i in 1..=isqrt {
        small_s[i - 1] = sieve.sum(i) as usize;
    } */
    let small_s = (1..=isqrt).map(|i| sieve.sum(i) as usize).collect_vec();
    count += isqrt;
    dbg!(count);
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
    let pi_4th_root = primes.partition_point(|p| p.pow(2) <= isqrt);
    let pi_cbrt = primes.partition_point(|p| p.pow(3) <= x);
    dbg!(pi_4th_root, pi_cbrt);
    for (i, &p) in primes[pi_4th_root..pi_cbrt].iter().enumerate() {
        let xp = x / p;
        let xpp = xp / p;
        let sp = small_s[p - 2];
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
        res -= large_s[p - 1] - small_s[p - 1] + 1;
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

const fn sum_n<const MOD: u128>(x: u128) -> u128 {
    let x = x % (MOD << 1);
    (if x & 1 == 0 {
        (x >> 1) * (x + 1)
    } else {
        x.div_ceil(2) * x
    }) % MOD
}
#[must_use]
pub fn lucy_sum<const MOD: u128>(x: u128) -> FIArrayU128 {
    let mut s = FIArrayU128::new(x);
    let keys = FIArrayU128::keys(x).collect_vec();

    unsafe { core::hint::assert_unchecked(s.arr.len() == keys.len()) };
    for (i, &v) in keys.iter().enumerate() {
        s.arr[i] = (sum_n::<MOD>(v) + MOD - 1) % MOD;
    }
    unsafe { core::hint::assert_unchecked(s.arr.len() as u128 > x.isqrt()) };
    for p in 2..=x.isqrt() {
        if s.arr[p as usize - 1] == s.arr[p as usize - 2] {
            continue;
        }
        let sp = s.arr[p as usize - 2];
        for (i, &v) in keys.iter().enumerate().rev() {
            if v < p * p {
                break;
            }
            s.arr[i] += MOD - ((p % MOD) * (s[v / p] + MOD - sp) % MOD) % MOD;
            s.arr[i] %= MOD;
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

    /*println!("prime counting using the logarithm of the zeta function:");
    let start = Instant::now();
    let count = log_zeta_reordered(N as _)[N as _]; // n^(2/3)
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    let start = Instant::now();
    let count = log_zeta(N as _)[N as _]; // n^(2/3)
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");*/

    /* println!("legendre:");
    let start = Instant::now();
    let count = legendre(N as _);
    let end = start.elapsed();
    println!("res = {count}, took {end:?}"); */

    println!("standard-ish lucy");

    let start = Instant::now();
    let count = prime_pi_fenwick(N as _);
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    let start = Instant::now();
    let count = prime_pi(N as _);
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");
    /* let start = Instant::now();
    let count = lucy_dumber(N as _)[N as _];
    let end = start.elapsed();
    println!("res = {count}, took {end:?}"); */

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

    let start = Instant::now();
    let count = lucy_wheel(N as _)[N as _];
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    let start = Instant::now();
    let count = lucy_wheel210(N)[N];
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    let start = Instant::now();
    let count = lucy_alt(N)[N];
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    let start = Instant::now();
    let count = lucy_non_fiarray_alt(N as _);
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    println!("strength reduced:");
    let start = Instant::now();
    let count = lucy_strengthreduce(N)[N];
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    let start = Instant::now();
    let count = lucy_strengthreduce_alt(N)[N];
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    println!("fast divide:");
    let start = Instant::now();
    let count = lucy_fastdivide(N as _)[N as _];
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    let start = Instant::now();
    let count = lucy_fastdivide_wheel(N as _)[N as _];
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    let start = Instant::now();
    let count = lucy_fastdivide_wheel210(N as _)[N as _];
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    let start = Instant::now();
    let count = lucy_fastdivide_alt(N as _)[N as _];
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

    let mut count_sum = 0usize;
    let mut count_add = 0usize;
    for p in 2..=sqrtx {
        if !sieve_raw.get(p) {
            count_sum += 1;
            let sp = sieve.sum(p - 1) as usize;
            let lim = (x / y).min(x / (p * p));
            for i in 1..=lim {
                let xip = x / (i * p);
                s.arr[len - i] -= if xip <= y {
                    count_sum += 1;
                    sieve.sum(xip) as usize
                } else {
                    s[xip]
                } - sp;
            }
            for j in (p * p..=y).step_by(p) {
                if !sieve_raw.get(j) {
                    sieve_raw.set(j);
                    sieve.add(j, -1);
                    count_add += 1;
                }
            }
        }
    }
    for (i, v) in FIArray::keys(x).take_while(|&v| v <= y).enumerate() {
        s.arr[i] = sieve.sum(v) as usize;
        count_sum += 1;
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

    let mut count_sum = 0usize;
    let mut count_add = 0usize;
    for p in 2..=sqrtx {
        if sieve_raw.get(p) {
            continue;
        }
        let sp = sieve.sum(p - 1) as usize;
        count_sum += 1;

        let lim = (x / y).min(x / (p * p));
        let xp = x / p;
        s.arr[len - 1] -= if xp <= y {
            count_sum += 1;
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
                count_sum += 1;
                sieve.sum(xip) as usize
            } else {
                s[xip]
            } - sp;
        }
        for j in (p * p..=y).step_by(p) {
            if !sieve_raw.get(j) {
                sieve_raw.set(j);
                sieve.add(j, -1);
                count_add += 1;
            }
        }
    }
    /* dbg!((count_sum, count_add));
       dbg!(y * (y as f64).ln().ln() as usize);
    */
    s[x]
}

// based on https://codeforces.com/blog/entry/91632?#comment-802482, https://codeforces.com/blog/entry/117783
// O(n^(2/3)) time, O(n^(1/2)) space. Pretty slow, despite noticeably superior time complexity.
// Likely due to repeated calls to dirichlet_mul, which is not particularly fast.
// Moreover, this function needs 2 times more memory than the O(n^0.75/log(n)) lucy_hedgehog based functions.
// O(n^0.75/log(n)) with low constant factors is better than O(n^(2/3)) with medium constant factors out to quite large n
// uses the fact that the coefficients of the logarithm of the DGF of u(n) = 1 (i.e. the zeta function)
// are exactly 1/k for p^k for some prime p, and 0 otherwise.
// Note: similarly to lucy_hedgehog, this code can be adapted to calculate the sum of totally multiplicative functions
// over the primes, though tbh you should probably just use lucy's algorithm for that.
// TODO: try to speed up the convolution steps more somehow, as they are the main bottleneck
const fn icbrt(x: usize) -> usize {
    let mut rt = 0;
    let mut rt_squared = 0;
    let mut rt_cubed = 0;
    while rt_cubed <= x {
        rt += 1;
        rt_squared += 2 * rt - 1;
        rt_cubed += 3 * rt_squared - 3 * rt + 1;
    }
    rt - 1
}
#[must_use]
pub fn log_zeta(n: usize) -> FIArray {
    const INVS: [usize; 6] = [0, 60, 30, 20, 15, 12];
    let rt = n.isqrt();
    let mut zeta = FIArray::unit(n as _);
    let len = zeta.arr.len();

    let mut buffer = zeta.clone();

    let mut ret = FIArray::new(n);
    let x = icbrt(rt) * (n as f64).ln() as usize;
    // remove contributions of small primes
    for p in 2..x {
        let val = zeta.arr[p - 1] - 1;
        if val == 0 {
            //not prime
            continue;
        }
        ret.arr[p - 1] = 1;
        for (i, nk) in buffer.arr.iter().enumerate().rev() {
            if i < p {
                break;
            }
            zeta.arr[i] -= zeta[nk / p];
        }
        zeta.arr[p - 1] = 1;
    }
    zeta.arr[..x - 1].fill(0);

    for i in x..=len {
        zeta.arr[i - 1] -= 1;
    }

    // zeta now equals zeta_t - 1
    // compute log(zeta_t) using log(x + 1) = x^5 / 5 - x^4 / 4 + x^3 / 3 - x^2 / 2 + x
    // x is zeta_t - 1.
    // in order to not have to deal with rational numbers, we compute 60 * log(zeta_t)
    // and adjust later
    // the contributions of x^4 and x^5 are 0 for essentially all reasonable n

    for i in x..=len {
        ret.arr[i - 1] = zeta.arr[i - 1] * INVS[1];
    }

    let mut pow_zeta = dirichlet_mul_zero_prefix(&zeta, &zeta, n, x - 1, x - 1);
    for i in x..=len {
        ret.arr[i - 1] -= pow_zeta.arr[i - 1] * INVS[2];
    }

    dirichlet_mul_zero_prefix_with_buffer(
        &zeta,
        &pow_zeta,
        n,
        &mut buffer,
        x - 1,
        pow_zeta.arr.iter().take_while(|&&e| e == 0).count(),
    );
    core::mem::swap(&mut pow_zeta.arr, &mut buffer.arr);

    for i in x..=len {
        ret.arr[i - 1] += pow_zeta.arr[i - 1] * INVS[3];
    }

    dirichlet_mul_zero_prefix_with_buffer(
        &zeta,
        &pow_zeta,
        n,
        &mut buffer,
        x - 1,
        pow_zeta.arr.iter().take_while(|&&e| e == 0).count(),
    );
    core::mem::swap(&mut pow_zeta.arr, &mut buffer.arr);

    for i in x..=len {
        ret.arr[i - 1] -= pow_zeta.arr[i - 1] * INVS[4];
    }

    dirichlet_mul_zero_prefix_with_buffer(
        &zeta,
        &pow_zeta,
        n,
        &mut buffer,
        x - 1,
        pow_zeta.arr.iter().take_while(|&&e| e == 0).count(),
    );
    core::mem::swap(&mut pow_zeta.arr, &mut buffer.arr);

    for i in x..=len {
        ret.arr[i - 1] += pow_zeta.arr[i - 1] * INVS[5];
    }

    // correction phase: get rid of contributions of prime powers
    for i in (x + 1..=len).rev() {
        ret.arr[i - 1] -= ret.arr[i - 2];
    }

    for x in x..=rt {
        let v = ret.arr[x - 1] / 60;
        let mut e = 1;
        let mut pv = v;
        let mut px = x;
        while px <= n / x {
            e += 1;
            px *= x;
            pv *= v;

            ret[px] -= pv * INVS[e];
        }
    }
    for i in 1..len {
        if i >= x - 1 {
            ret.arr[i] /= 60;
        }
        ret.arr[i] += ret.arr[i - 1];
    }
    ret
}

// identical to log_zeta, but convolutions are reordered in order to maximise shared 0 prefix, minor speedup for large n
// 1e17: 14411.9862525s
// 1e16: 2035.5664288s
// 1e15: 362.4408137s
// 1e14: 65.0882843s
// 1e13: 12.9792081s
// 1e12: 2.5241903s
// 1e11: 477.0335ms
// 1e10: 90.5864ms
// 1e9: 17.5497ms
// 1e8: 3.5963ms
// can try to write version which only computes final result:
#[must_use]
pub fn log_zeta_reordered(n: usize) -> FIArray {
    const INVS: [usize; 6] = [0, 60, 30, 20, 15, 12];
    let rt = n.isqrt();
    let mut zeta = FIArray::unit(n as _);
    let len = zeta.arr.len();

    let mut buffer = zeta.clone();

    let mut ret = FIArray::new(n);
    let x = {
        let mut x = 2;
        let mut x_cubed = 8;
        while x_cubed <= rt {
            x += 1;
            x_cubed += 3 * x * (x - 1) + 1;
        }
        x
    } * (n as f64).ln() as usize; // since primes are sparse, can afford to increase x by logarithmic factor without hurting complexity
    // remove contributions of small primes (first ~n^1/6 of them)
    for p in 2..x {
        let val = zeta.arr[p - 1] - 1;
        if val == 0 {
            //not prime
            continue;
        }
        ret.arr[p - 1] = 1;
        for (i, nk) in buffer.arr.iter().enumerate().rev() {
            if i < p {
                break;
            }
            zeta.arr[i] -= zeta[nk / p];
        }
        zeta.arr[p - 1] = 1;
    }
    //let prime_count = zeta.arr[..x - 1].iter().filter(|&&e| e == 1).count();
    zeta.arr[..x - 1].fill(0);

    for i in x..=len {
        zeta.arr[i - 1] -= 1;
    }

    // zeta now equals zeta_t - 1
    // compute log(zeta_t) using log(x + 1) = x^5 / 5 - x^4 / 4 + x^3 / 3 - x^2 / 2 + x
    // x is zeta_t - 1.
    // in order to not have to deal with rational numbers, we compute 60 * log(zeta_t)
    // and adjust later
    // the contributions of x^4 and x^5 are 0 for essentially all reasonable n

    for i in x..=len {
        ret.arr[i - 1] = zeta.arr[i - 1] * INVS[1];
    }
    //let start = std::time::Instant::now();
    let pow_zeta = dirichlet_mul_zero_prefix(&zeta, &zeta, n, x - 1, x - 1);
    //dbg!(start.elapsed());
    let z2_pref = pow_zeta.arr.iter().take_while(|&&e| e == 0).count();

    for i in z2_pref..=len {
        ret.arr[i - 1] -= pow_zeta.arr[i - 1] * INVS[2];
    }

    if z2_pref * z2_pref < n {
        dirichlet_mul_zero_prefix_with_buffer(
            &pow_zeta,
            &pow_zeta,
            n,
            &mut buffer,
            z2_pref,
            z2_pref,
        );

        for i in z2_pref..=len {
            ret.arr[i - 1] -= buffer.arr[i - 1] * INVS[4];
        }
    }
    if (x - 1) * z2_pref < n {
        dirichlet_mul_zero_prefix_with_buffer(&zeta, &pow_zeta, n, &mut buffer, x - 1, z2_pref);

        for i in z2_pref..=len {
            ret.arr[i - 1] += buffer.arr[i - 1] * INVS[3];
        }
        let z3_pref = buffer.arr.iter().take_while(|&&e| e == 0).count();

        if z2_pref * z3_pref < n {
            dirichlet_mul_zero_prefix_with_buffer(
                &pow_zeta, &buffer, n, &mut zeta, z2_pref, z3_pref,
            );
            for i in z3_pref..=len {
                ret.arr[i - 1] += zeta.arr[i - 1] * INVS[5];
            }
        }
    }
    //dbg!(prime_count + ret.arr[len - 1] / 60); // approximate final result
    // correction phase: get rid of contributions of prime powers
    for i in (x + 1..=len).rev() {
        ret.arr[i - 1] -= ret.arr[i - 2];
    }

    for x in x..=rt {
        let v = ret.arr[x - 1];
        if v == 0 {
            continue;
        }
        let mut e = 1;
        //let mut pv = v;
        let mut px = x;
        while px <= n / x {
            e += 1;
            px *= x;
            //pv *= v;

            ret[px] -= /* pv * */ INVS[e];
        }
    }
    for i in 1..len {
        if i >= x - 1 {
            ret.arr[i] /= 60;
        }
        ret.arr[i] += ret.arr[i - 1];
    }
    ret
}

// TODO: fix
/* pub fn log_zeta_reordered_single(n: usize) -> usize {
    const INVS: [usize; 6] = [0, 60, 30, 20, 15, 12];
    let rt = n.isqrt();
    let mut zeta = FIArray::unit(n as _);
    let len = zeta.arr.len();

    let mut buffer = zeta.clone();

    let mut ret = 0;
    let x = {
        let mut x = 2;
        let mut x_cubed = 8;
        while x_cubed <= rt {
            x += 1;
            x_cubed += 3 * x * (x - 1) + 1;
        }
        x
    } * (n as f64).ln() as usize; // since primes are sparse, can afford to increase x by logarithmic factor without hurting complexity
    // remove contributions of small primes (first ~n^1/6 of them)
    for p in 2..x {
        let val = zeta.arr[p - 1] - 1;
        if val == 0 {
            //not prime
            continue;
        }
        ret += 1;
        for (i, nk) in buffer.arr.iter().enumerate().rev() {
            if i < p {
                break;
            }
            zeta.arr[i] -= zeta[nk / p];
        }
        zeta.arr[p - 1] = 1;
    }
    //let prime_count = zeta.arr[..x - 1].iter().filter(|&&e| e == 1).count();
    zeta.arr[..x - 1].fill(0);

    for i in x..=len {
        zeta.arr[i - 1] -= 1;
    }

    // zeta now equals zeta_t - 1
    // compute log(zeta_t) using log(x + 1) = x^5 / 5 - x^4 / 4 + x^3 / 3 - x^2 / 2 + x
    // x is zeta_t - 1.
    // in order to not have to deal with rational numbers, we compute 60 * log(zeta_t)
    // and adjust later
    // the contributions of x^4 and x^5 are 0 for essentially all reasonable n

    ret += zeta.arr[len - 1] * INVS[1];

    //let start = std::time::Instant::now();
    let pow_zeta = dirichlet_mul_zero_prefix(&zeta, &zeta, n, x - 1, x - 1);
    //dbg!(start.elapsed());
    let z2_pref = pow_zeta.arr.iter().take_while(|&&e| e == 0).count();

    ret -= pow_zeta.arr[len - 1] * INVS[2];

    if z2_pref * z2_pref < n {
        ret -= INVS[4]
            * dirichlet_mul_single_zero_prefix_usize(&pow_zeta, &pow_zeta, n, z2_pref, z2_pref);
    }
    if (x - 1) * z2_pref < n {
        //if (x - 1) * z2_pref * z2_pref < n {
        println!("hello 1");
        dirichlet_mul_zero_prefix_with_buffer(&zeta, &pow_zeta, n, &mut buffer, x - 1, z2_pref);

        ret += buffer.arr[len - 1] * INVS[3];

        let z3_pref = buffer.arr.iter().take_while(|&&e| e == 0).count();

        if z2_pref * z3_pref < n {
            ret += dirichlet_mul_single_zero_prefix_usize(&pow_zeta, &buffer, n, z2_pref, z3_pref)
                * INVS[5];
        }
        /* } else {
        println!("hello 2");

        ret +=
            dirichlet_mul_single_zero_prefix_usize(&zeta, &pow_zeta, n, x - 1, z2_pref) * INVS[3];
        } */
    }
    //dbg!(prime_count + ret.arr[len - 1] / 60); // approximate final result
    // correction phase: get rid of contributions of prime powers
    let primes = sift(rt as _);

    for x in primes
        .into_iter()
        .filter_map(|p| (p as usize >= x).then_some(p as usize))
    {
        let mut e = 1;
        let mut px = x;
        while px <= n / x {
            e += 1;
            px *= x;
            ret -= INVS[e];
        }
    }

    ret / 60
}
 */
#[must_use]
pub fn dirichlet_mul_zero_prefix(
    F: &FIArray,
    G: &FIArray,
    n: usize,
    prefix_f: usize,
    prefix_g: usize,
) -> FIArray {
    assert!(prefix_f > 0);
    assert!(prefix_g > 0);
    assert!(prefix_f <= prefix_g);
    let mut H = FIArray::new(n as _);
    let len = H.arr.len();
    let rt_n = n.isqrt();

    let real_pref_f = if prefix_f <= rt_n {
        prefix_f
    } else {
        n / (len - prefix_f) - 1
    };
    let real_pref_g = if prefix_g <= rt_n {
        prefix_g
    } else {
        n / (len - prefix_g) - 1
    };

    if real_pref_f * real_pref_g >= n {
        return H;
    }

    let to_ord = |x| {
        if x <= rt_n { x } else { len + 1 - (n / x) }
    };
    let mut propogate = |(x0, x1), (y0, y1), (z0, z1)| {
        let f_x1 = F.arr[x1 - 1];
        let g_y1 = G.arr[y1 - 1];
        let f_x0_1 = F.arr.get(x0 - 2).copied().unwrap_or_default();
        let g_y0_1 = G.arr.get(y0 - 2).copied().unwrap_or_default();

        let t = (f_x1 - f_x0_1) * (g_y1 - g_y0_1);
        H.arr[z0 - 1] += t;
        if let Some(v) = H.arr.get_mut(z1) {
            *v -= t;
        }
    };

    let prefix_h = to_ord(real_pref_f * real_pref_g) - 1;

    for k in prefix_f..=len {
        let z = len + 1 - k;
        if k >= prefix_h {
            for x in prefix_f.. {
                let y_lo_ord = 1 + to_ord(x).max(to_ord(z));
                let y_hi_ord = to_ord(n / (x * z));
                if y_hi_ord < y_lo_ord {
                    break;
                }
                propogate((x, x), (y_lo_ord, y_hi_ord), (k, k));
                propogate((y_lo_ord, y_hi_ord), (x, x), (k, k));
            }
        }

        let x = k;
        for y in prefix_f..k {
            let z_lo_ord = to_ord(x * y);
            let z_hi_ord = to_ord(n / x);
            if z_hi_ord < z_lo_ord {
                break;
            }
            propogate((x, x), (y, y), (z_lo_ord, z_hi_ord));
            propogate((y, y), (x, x), (z_lo_ord, z_hi_ord));
        }

        if prefix_g <= x && x <= rt_n {
            propogate((x, x), (x, x), (to_ord(x * x), len));
        }
    }

    for i in 1..len {
        H.arr[i] += H.arr[i - 1];
    }
    H
}

pub fn dirichlet_mul_zero_prefix_with_buffer(
    F: &FIArray,
    G: &FIArray,
    n: usize,
    H: &mut FIArray,
    prefix_f: usize,
    prefix_g: usize,
) {
    assert!(prefix_f > 0);
    assert!(prefix_g > 0);
    assert!(prefix_f <= prefix_g);

    H.arr.fill(0);
    let len = H.arr.len();
    if prefix_f == len || prefix_g == len {
        return;
    }
    let rt_n = n.isqrt();

    let real_pref_f = if prefix_f <= rt_n {
        prefix_f
    } else {
        n / (len - prefix_f) - 1
    };
    let real_pref_g = if prefix_g <= rt_n {
        prefix_g
    } else {
        n / (len - prefix_g) - 1
    };

    if real_pref_f * real_pref_g >= n {
        return;
    }

    let to_ord = |x| {
        if x <= rt_n { x } else { len + 1 - (n / x) }
    };
    let mut propogate = |(x0, x1), (y0, y1), (z0, z1)| {
        let f_x1 = F.arr[x1 - 1];
        let g_y1 = G.arr[y1 - 1];
        let f_x0_1 = F.arr.get(x0 - 2).copied().unwrap_or_default();
        let g_y0_1 = G.arr.get(y0 - 2).copied().unwrap_or_default();

        let t = (f_x1 - f_x0_1) * (g_y1 - g_y0_1);
        H.arr[z0 - 1] += t;
        if let Some(v) = H.arr.get_mut(z1) {
            *v -= t;
        }
    };
    let prefix_h = to_ord(real_pref_f * real_pref_g) - 1;

    for k in prefix_f..=len {
        let z = len + 1 - k;
        if k >= prefix_h {
            for x in prefix_f.. {
                let y_lo_ord = 1 + to_ord(x).max(to_ord(z));
                let y_hi_ord = to_ord(n / (x * z));
                if y_hi_ord < y_lo_ord {
                    break;
                }
                propogate((x, x), (y_lo_ord, y_hi_ord), (k, k));
                propogate((y_lo_ord, y_hi_ord), (x, x), (k, k));
            }
        }
        let x = k;
        for y in prefix_f..k {
            let z_lo_ord = to_ord(x * y);
            let z_hi_ord = to_ord(n / x);
            if z_hi_ord < z_lo_ord {
                break;
            }
            propogate((x, x), (y, y), (z_lo_ord, z_hi_ord));
            propogate((y, y), (x, x), (z_lo_ord, z_hi_ord));
        }

        if prefix_g <= x && x <= rt_n {
            propogate((x, x), (x, x), (to_ord(x * x), len));
        }
    }

    for i in 1..len {
        H.arr[i] += H.arr[i - 1];
    }
}

/// Calculate the logarithmic integral using
/// Ramanujan's formula:
/// <https://en.wikipedia.org/wiki/Logarithmic_integral_function#Series_representation>
///
#[must_use]
pub fn li(x: f64) -> f64 {
    const GAMMA: f64 = 0.577215664901532860606512090082402431_f64;
    if x <= 1. {
        return 0.;
    }
    let mut sum = 0.;
    let mut inner_sum = 0.;
    let mut factorial = 1.;
    let mut p = -1.;
    let mut q;
    let mut power2 = 1.;
    let logx = x.ln();
    let mut k = 0;

    // The condition n < ITERS is required in case the computation
    // does not converge. This happened on Linux i386 where
    // the precision of the libc math functions is very limited.
    for n in 1..1000 {
        p *= -logx;
        factorial *= f64::from(n);
        q = factorial * power2;
        power2 *= 2.;

        while k <= (n - 1) / 2 {
            inner_sum += f64::from(2 * k + 1).recip();
            k += 1;
        }

        let old_sum = sum;
        sum += (p / q) * inner_sum;

        // Not converging anymore
        if (sum - old_sum).abs() <= f64::EPSILON {
            break;
        }
    }

    GAMMA + logx.ln() + x.sqrt() * sum
}

/// Calculate the Eulerian logarithmic integral which is a very
/// accurate approximation of the number of primes <= x.
/// Li(x) > pi(x) for 24 <= x <= ~ 10^316
///
#[must_use]
pub fn Li(x: usize) -> usize {
    const li2: f64 = 1.045163780117492784844588889194613136_f64;

    if x <= 2 {
        0
    } else {
        (li(x as f64) - li2) as usize
    }
}

#[must_use]
pub fn R(x: usize) -> usize {
    /* let Li = |x: f64| {
        const li2: f64 = 1.045163780117492784844588889194613136_f64;
        if x <= 2. { 0. } else { li(x as f64) - li2 }
    }; */
    let x = x as f64;
    let mobius = mobius_sieve(50);
    let mut sum = 0.;
    for i in 1..50 {
        let inv_i = (i as f64).recip();
        let li = li(x.powf(inv_i));
        if li == 0. {
            break;
        }
        sum += f64::from(mobius[i]) * inv_i * li;
    }
    sum as usize
}
