use std::{collections::HashSet, time::Instant};

use fastdivide::DividerU64;
use itertools::Itertools;

use crate::utils::{
    FIArray::{FIArray, FIArrayU64, FIArrayU128},
    bit_array::BitArray,
    fenwick::FenwickTree,
    multiplicative_function_summation::mobius_sieve,
    prime_sieves::{WHEEL_2_3_5, WHEEL_2_3_5_7, sift},
};

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
// 1e17: prime counting took 4042.9456095s: 2623557157654233
// 1e16: prime counting took 525.8536176s: 279238341033925
// 1e15: prime counting took 102.6765764s: 29844570422669
// 1e14: prime counting took 19.8688139s: 3204941750802
// 1e13: prime counting took 3.527102s: 346065536839
// 1e12: prime counting took 649.5083ms: 37607912018
// 1e11: prime counting took 117.2172ms: 4118054813
// 1e10: prime counting took 23.0519ms: 455052511
// 1e9: prime counting took 4.6073ms: 50847534
#[must_use]
pub fn lucy(x: usize) -> FIArray {
    let mut s = FIArray::new(x);
    let keys = FIArray::keys(x).collect_vec();
    unsafe { core::hint::assert_unchecked(s.arr.len() == keys.len()) };
    for (i, v) in keys.iter().enumerate() {
        s.arr[i] = (v + 1) >> 1;
    }
    s.arr[0] = 0;
    s.arr[2] = 2; // deal with 3 separately
    let mut pp = 1;
    unsafe { core::hint::assert_unchecked(s.arr.len() > x.isqrt()) };
    for p in (3..=x.isqrt()).step_by(2) {
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
pub fn lucy_alt(x: usize) -> FIArray {
    let primes = sift(x.isqrt() as u64);
    let mut s = FIArray::new(x);
    let keys = FIArray::keys(x).collect_vec();

    unsafe { core::hint::assert_unchecked(s.arr.len() == keys.len()) };
    for (i, v) in keys.iter().enumerate() {
        s.arr[i] = (v + 1) >> 1;
    }
    s.arr[0] = 0;
    s.arr[2] = 2;
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
pub fn lucy_wheel(x: usize) -> FIArray {
    let mut s = FIArray::new(x);
    let keys = FIArray::keys(x).collect_vec();

    unsafe { core::hint::assert_unchecked(s.arr.len() == keys.len()) };
    for (i, v) in keys.iter().enumerate() {
        s.arr[i] = (v + 1) >> 1;
    }
    s.arr[0] = 0;
    s.arr[2] = 2;
    unsafe { core::hint::assert_unchecked(s.arr.len() > x.isqrt()) };
    let lim = x.isqrt();
    assert!(lim >= 5);
    for p in [3, 5] {
        let sp = s.arr[p - 2];

        for (i, &v) in keys.iter().enumerate().rev() {
            if v < p * p {
                break;
            }
            s.arr[i] -= s[v / p] - sp;
        }
    }
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
    s.arr[2] = 2;
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
    const N: usize = 1e13 as usize;

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

    println!("lucy fenwick:");
    let start = Instant::now();
    let count = lucy_fenwick_trick(N as _);
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    let start = Instant::now();
    let count = lucy_trick(N as _);
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");
    let start = Instant::now();
    let count = lucy_fenwick(N as _)[N as _];
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

    /* println!("legendre:");
    let start = Instant::now();
    let count = legendre(N as _);
    let end = start.elapsed();
    println!("res = {count}, took {end:?}"); */

    println!("standard-ish lucy");
    let start = Instant::now();
    let count = lucy_dumber(N as _)[N as _];
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    let start = Instant::now();
    let count = lucy(N as _)[N as _];
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
            let mut j = p * p;
            while j <= y {
                if !sieve_raw.get(j) {
                    sieve_raw.set(j);
                    sieve.add(j, -1);
                    count_add += 1;
                }
                j += p;
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

#[must_use]
pub fn lucy_trick(x: usize) -> usize {
    let mut s = FIArray::new(x);
    let sqrtx = x.isqrt();
    let len = s.arr.len();
    let y = sqrtx;
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
        count_sum += 1;

        let sp = sieve.sum(p - 1) as usize;
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
    dbg!(y * (y as f64).ln().ln() as usize); */

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
#[must_use]
pub fn log_zeta(n: usize) -> FIArray {
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
    } * (n as f64).ln() as usize;
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
