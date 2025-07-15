use std::time::Instant;

use fastdivide::DividerU64;
use itertools::Itertools;

use crate::utils::{
    FIArray::{FIArray, FIArrayU64, FIArrayU128},
    prime_sieves::{WHEEL_2_3_5, WHEEL_2_3_5_7, sift},
};

// O(x^(3/4)) prime counting function
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
        s.arr[i] = v - 1;
    }
    unsafe { core::hint::assert_unchecked(s.arr.len() > x.isqrt()) };
    for p in 2..=x.isqrt() {
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
pub fn lucy_alt(x: usize) -> FIArray {
    let primes = sift(x.isqrt() as u64 + 1);
    let mut s = FIArray::new(x);
    let keys = FIArray::keys(x).collect_vec();

    unsafe { core::hint::assert_unchecked(s.arr.len() == keys.len()) };
    for (i, v) in keys.iter().enumerate() {
        s.arr[i] = v - 1;
    }
    for p in primes {
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
        s.arr[i] = v - 1;
    }
    unsafe { core::hint::assert_unchecked(s.arr.len() > x.isqrt()) };
    let lim = x.isqrt();
    assert!(lim >= 5);
    for p in [2, 3, 5] {
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
        s.arr[i] = v - 1;
    }
    unsafe { core::hint::assert_unchecked(s.arr.len() as u64 > x.isqrt()) };
    for p in 2..=x.isqrt() {
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
pub fn lucy_fastdivide_alt(x: u64) -> FIArrayU64 {
    let primes = sift(x.isqrt() + 1);

    let mut s = FIArrayU64::new(x);
    let keys = FIArrayU64::keys(x).collect_vec();

    unsafe { core::hint::assert_unchecked(s.arr.len() == keys.len()) };
    for (i, v) in keys.iter().enumerate() {
        s.arr[i] = v - 1;
    }
    unsafe { core::hint::assert_unchecked(s.arr.len() as u64 > x.isqrt()) };
    for p in primes {
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
        s.arr[i] = v - 1;
    }
    unsafe { core::hint::assert_unchecked(s.arr.len() as u64 > x.isqrt()) };
    let lim = x.isqrt();
    assert!(lim >= 5);
    for p in [2, 3, 5] {
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
        s.arr[i] = v - 1;
    }
    unsafe { core::hint::assert_unchecked(s.arr.len() as u64 > x.isqrt()) };
    let lim = x.isqrt();
    assert!(lim >= 7);
    for p in [2, 3, 5, 7] {
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
        s.arr[i] = v - 1;
    }
    unsafe { core::hint::assert_unchecked(s.arr.len() > x.isqrt()) };
    for p in 2..=x.isqrt() {
        if s.arr[p - 1] == s.arr[p - 2] {
            continue;
        }
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
pub fn lucy_strengthreduce_alt(x: usize) -> FIArray {
    let primes = sift(x.isqrt() as u64 + 1);
    let mut s = FIArray::new(x);
    let keys = FIArray::keys(x).collect_vec();

    unsafe { core::hint::assert_unchecked(s.arr.len() == keys.len()) };
    for (i, v) in keys.iter().enumerate() {
        s.arr[i] = v - 1;
    }
    unsafe { core::hint::assert_unchecked(s.arr.len() > x.isqrt()) };
    for p in primes {
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
fn lucy_dumber(x: usize) -> FIArray {
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
    const N: usize = 4e7 as _;
    let start = Instant::now();
    let s = lucy(N as _);
    let count = s[N as _];
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    let start = Instant::now();
    let s2 = lucy_wheel(N as _);
    let count = s2[N as _];
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
    dbg!(s == s2);

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
    let s = lucy_fastdivide(N as _);
    let count = s[N as _];
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    let start = Instant::now();
    let s2 = lucy_fastdivide_wheel(N as _);
    let count = s2[N as _];
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

    dbg!(s == s2);
}
