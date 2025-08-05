use std::time::Instant;

use fastdivide::DividerU64;
use itertools::Itertools;

use crate::utils::{
    FIArray::{FIArray, FIArrayU64, FIArrayU128},
    bit_array::BitArray,
    fenwick::FenwickTree,
    multiplicative_function_summation::{dirichlet_mul_usize, dirichlet_mul_with_buffer_usize},
    prime_sieves::{WHEEL_2_3_5, WHEEL_2_3_5_7, sift},
};

// repeated convolution of the prefix sum representation of u with mu_p for p below sqrt(n)
// I guess this is essentially legendre's formula for prime counting, implemented using bottom-up dp
// not efficient, lucy is essentially a smarter version of this, reducing the complexity from O(n/logn) to O(n^0.75/logn)
pub fn test(x: usize) -> usize {
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
// 1e16: prime counting took 525.8536176s: 279238341033925
// 1e15: prime counting took 102.6765764s: 29844570422669
// 1e14: prime counting took 19.8688139s: 3204941750802
// 1e13: prime counting took 3.527102s: 346065536839
// 1e12: prime counting took 649.5083ms: 37607912018
// 1e11: prime counting took 117.2172ms: 4118054813
// 1e10: prime counting took 23.0519ms: 455052511
// 1e9: prime counting took 4.6073ms: 50847534
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
    const N: usize = 1e11 as usize;

    let start = Instant::now();
    let count = lucy_trick(N as _);
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    println!("prime counting using the logarithm of the zeta function:");
    let start = Instant::now();
    let count = log_zeta_alt(N as _)[N as _];
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    let start = Instant::now();
    let count = log_zeta(N as _)[N as _];
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    println!("lucy fenwick:");
    let start = Instant::now();
    let count = lucy_fenwick_trick(N as _);
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    let start = Instant::now();
    let count = lucy_fenwick(N as _)[N as _];
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    println!("standard-ish lucy");
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

//never faster for me, though asymptotically better: O(x^(2/3) (loglogx)^(1/3)) vs O(x^(3/4))
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
            let mut j = p * p;
            while j <= y {
                if !sieve_raw.get(j) {
                    sieve_raw.set(j);
                    sieve.add(j, -1);
                }
                j += p;
            }
        }
    }
    for (i, v) in FIArray::keys(x).take_while(|&v| v <= y).enumerate() {
        s.arr[i] = sieve.sum(v) as usize;
    }
    s
}

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
    s[x]
}

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
    s[x]
}

// based on https://codeforces.com/blog/entry/91632?#comment-802482
// O(n^(2/3)) time, O(n^(1/2)) space. Pretty slow, despite noticeably superior time complexity.
// Likely due to repeated calls to dirichlet_mul, which is not particularly fast.
// Moreover, this function needs 4 times more memory than the O(n^0.75/log(n)) lucy_hedgehog based functions.
// O(n^0.75/log(n)) with low constant factors is better than O(n^(2/3)) with medium constant factors out to quite large n
// uses the fact that the coefficients of the logarithm of the DGF of u(n) = 1 (i.e. the zeta function)
// are exactly 1/k for p^k for some prime p, and 0 otherwise.
// Note: similarly to lucy_hedgehog, this code can be adapted to calculate the sum of totally multiplicative functions
// over the primes, though tbh you should probably just use lucy's algorithm for that.
// TODO: try to speed up the convolution steps more, as each of the arrays has a 0 prefix of length ~ n^(1/6)
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
    };
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
    // note that almost the entirety of the time spent by this function is in the following 4 convolutions.
    // literally 95%+ of time taken for large n
    let start = std::time::Instant::now();
    for i in x..=len {
        ret.arr[i - 1] = zeta.arr[i - 1] * INVS[1];
    }

    let mut pow_zeta = dirichlet_mul_zero_prefix(&zeta, &zeta, n, x - 1);
    for i in x..=len {
        ret.arr[i - 1] -= pow_zeta.arr[i - 1] * INVS[2];
    }

    dirichlet_mul_zero_prefix_with_buffer(&pow_zeta, &zeta, n, &mut buffer, x - 1);
    core::mem::swap(&mut pow_zeta, &mut buffer);

    for i in x..=len {
        ret.arr[i - 1] += pow_zeta.arr[i - 1] * INVS[3];
    }

    dirichlet_mul_zero_prefix_with_buffer(&pow_zeta, &zeta, n, &mut buffer, x - 1);
    core::mem::swap(&mut pow_zeta, &mut buffer);

    for i in x..=len {
        ret.arr[i - 1] -= pow_zeta.arr[i - 1] * INVS[4];
    }

    dirichlet_mul_zero_prefix_with_buffer(&pow_zeta, &zeta, n, &mut buffer, x - 1);
    core::mem::swap(&mut pow_zeta, &mut buffer);

    for i in x..=len {
        ret.arr[i - 1] += pow_zeta.arr[i - 1] * INVS[5];
    }
    let end = start.elapsed();
    println!("second phase took {end:?}");

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

pub fn dirichlet_mul_zero_prefix(F: &FIArray, G: &FIArray, n: usize, prefix: usize) -> FIArray {
    assert!(prefix > 0);
    /* assert!(F.arr.iter().take_while(|&&e| e == 0).count() >= prefix);
       assert!(G.arr.iter().take_while(|&&e| e == 0).count() >= prefix);
    */
    let mut H = FIArray::new(n as _);
    let len = H.arr.len();

    let rt_n = n.isqrt();

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

    for k in 2..=len {
        let z = len + 1 - k;
        for x in prefix.. {
            let y_lo_ord = 1 + to_ord(x).max(to_ord(z));
            let y_hi_ord = to_ord(n / (x * z));
            if y_hi_ord < y_lo_ord {
                break;
            }
            propogate((x, x), (y_lo_ord, y_hi_ord), (k, k));
            propogate((y_lo_ord, y_hi_ord), (x, x), (k, k));
        }

        let x = k;
        for y in prefix..k {
            let z_lo_ord = to_ord(x * y);
            let z_hi_ord = to_ord(n / x);
            if z_hi_ord < z_lo_ord {
                break;
            }
            propogate((x, x), (y, y), (z_lo_ord, z_hi_ord));
            propogate((y, y), (x, x), (z_lo_ord, z_hi_ord));
        }

        if prefix <= x && x <= rt_n {
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
    prefix: usize,
) {
    assert!(prefix > 0);

    H.arr.fill(0);
    let len = H.arr.len();

    let rt_n = n.isqrt();

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

    for k in 2..=len {
        let z = len + 1 - k;
        for x in prefix.. {
            let y_lo_ord = 1 + to_ord(x).max(to_ord(z));
            let y_hi_ord = to_ord(n / (x * z));
            if y_hi_ord < y_lo_ord {
                break;
            }
            propogate((x, x), (y_lo_ord, y_hi_ord), (k, k));
            propogate((y_lo_ord, y_hi_ord), (x, x), (k, k));
        }

        let x = k;
        for y in prefix..k {
            let z_lo_ord = to_ord(x * y);
            let z_hi_ord = to_ord(n / x);
            if z_hi_ord < z_lo_ord {
                break;
            }
            propogate((x, x), (y, y), (z_lo_ord, z_hi_ord));
            propogate((y, y), (x, x), (z_lo_ord, z_hi_ord));
        }

        if prefix <= x && x <= rt_n {
            propogate((x, x), (x, x), (to_ord(x * x), len));
        }
    }

    for i in 1..len {
        H.arr[i] += H.arr[i - 1];
    }
}

// worse time complexity, better performance.
// could try using an O(n^0.75) approach to the convolution, wouldn't hurt the complexity but might grant speed up - nvm doesn't help
pub fn log_zeta_alt(n: usize) -> FIArray {
    const INVS: [usize; 4] = [0, 6, 3, 2];
    let rt = n.isqrt();
    let mut zeta = FIArray::unit(n as _);
    let len = zeta.arr.len();

    let mut buffer = zeta.clone();

    let mut ret = FIArray::new(n);
    let x = rt.isqrt() + 1;
    // remove contributions of small primes up to n^(1/4)
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
    // compute log(zeta_t) using log(x + 1) = x^3 / 3 - x^2 / 2 + x
    // x is zeta_t - 1.
    // in order to not have to deal with rational numbers, we compute 6*log(zeta_t)
    // and adjust later
    for i in x..=len {
        ret.arr[i - 1] = zeta.arr[i - 1] * INVS[1];
    }

    let mut pow_zeta = dirichlet_mul_zero_prefix(&zeta, &zeta, n, x - 1);
    for i in x..=len {
        ret.arr[i - 1] -= pow_zeta.arr[i - 1] * INVS[2];
    }

    dirichlet_mul_zero_prefix_with_buffer(&pow_zeta, &zeta, n, &mut buffer, x - 1);
    core::mem::swap(&mut pow_zeta, &mut buffer);

    for i in x..=len {
        ret.arr[i - 1] += pow_zeta.arr[i - 1] * INVS[3];
    }

    // correction phase: get rid of contributions of prime powers
    for i in (x + 1..=len).rev() {
        ret.arr[i - 1] -= ret.arr[i - 2];
    }

    for x in x..=rt {
        let v = ret.arr[x - 1] / 6;
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
    for i in 1..ret.arr.len() {
        if i >= x - 1 {
            ret.arr[i] /= 6;
        }
        ret.arr[i] += ret.arr[i - 1];
    }
    ret
}
pub fn conv_with_buffer_zero_prefix(
    F: &FIArray,
    G: &FIArray,
    n: usize,
    H: &mut FIArray,
    prefix_f: usize,
    prefix_g: usize,
) {
    assert!(prefix_f > 0);
    assert!(prefix_g > 0);
    H.arr.fill(0);
    for (i, v) in FIArray::keys(n).enumerate() {
        if v <= prefix_f.min(prefix_g) {
            continue;
        }
        let vsqrt = v.isqrt();
        let mut h = 0;
        for i in (1 + prefix_f)..=vsqrt {
            h += G[v / i] * (F.arr[i - 1] - F.arr[i - 2]);
        }
        for i in (1 + prefix_g)..=vsqrt {
            h += F[v / i] * (G.arr[i - 1] - G.arr[i - 2]);
        }
        h -= F[vsqrt] * G[vsqrt];
        H.arr[i] = h;
    }
}

pub fn conv_zero_prefix(
    F: &FIArray,
    G: &FIArray,
    n: usize,
    prefix_f: usize,
    prefix_g: usize,
) -> FIArray {
    let mut H = FIArray::new(n);
    conv_with_buffer_zero_prefix(F, G, n, &mut H, prefix_f, prefix_g);
    H
}
