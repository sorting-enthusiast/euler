use itertools::Itertools;

use crate::utils::{
    FIArray::FIArray, multiplicative_function_summation::sum_n_usize, prime_sieves::sift,
};
const MOD: usize = 1e9 as _;
const N: usize = 1e12 as usize;
/*
pub fn sum_over_primes(
    x: usize,
    mut f: impl FnMut(usize) -> usize,
    mut F: impl FnMut(usize) -> usize,
) -> FIArray {
    let primes = sift(x.isqrt() as u64);
    let mut s = FIArray::new(x);
    let keys = FIArray::keys(x).collect_vec();

    unsafe { core::hint::assert_unchecked(s.arr.len() == keys.len()) };
    for (i, &v) in keys.iter().enumerate() {
        s.arr[i] = F(v);
        assert!(s.arr[i] < MOD);
    }
    unsafe { core::hint::assert_unchecked(s.arr.len() > x.isqrt()) };
    for p in primes {
        let p = p as usize;
        let minus_sp = MOD - s.arr[p - 2];
        for (i, &v) in keys.iter().enumerate().rev() {
            if v < p * p {
                break;
            }
            s.arr[i] += MOD - (f(p) * (s[v / p] + minus_sp) % MOD) % MOD;
            s.arr[i] %= MOD;
        }
    }
    s
}
 */
// finding the actual solution was quite easy,
// although I wasted an unbelievable amount of time debugging the code,
// turns out i was subtracting 1 from the sum in F wrong,
// not making sure the result mod MOD was nonnegative
pub fn main() {
    let start = std::time::Instant::now();

    let f = |p: usize| {
        let p = p % MOD;
        ((p * p) % MOD * p) % MOD
    };
    let F = |v| {
        let sn = sum_n_usize::<MOD>(v);
        ((sn * sn) % MOD + MOD - 1) % MOD
    };
    let primes = sift(N.isqrt() as u64);

    let mut s = FIArray::new(N);
    let keys = FIArray::keys(N).collect_vec();

    unsafe { core::hint::assert_unchecked(s.arr.len() == keys.len()) };
    for (i, &v) in keys.iter().enumerate() {
        s.arr[i] = F(v);
        assert!(s.arr[i] < MOD);
    }
    unsafe { core::hint::assert_unchecked(s.arr.len() > N.isqrt()) };
    for &p in &primes {
        let p = p as usize;
        let minus_sp = MOD - s.arr[p - 2];
        for (i, &v) in keys.iter().enumerate().rev() {
            if v < p * p {
                break;
            }
            s.arr[i] += MOD - (f(p) * (s[v / p] + minus_sp) % MOD) % MOD;
            s.arr[i] %= MOD;
        }
    }
    let mut sum = 0;
    for &p in &primes {
        let p = p as usize;
        sum += (f(p) * (s[N / p] + MOD - s.arr[p - 1]) % MOD) % MOD;
        sum %= MOD;
    }
    let end = start.elapsed();
    println!("{sum}, {end:?}");
}
//f(n) = sum over d|n, d<n of f(d)
//f(p^e) = 2^(e-1)
