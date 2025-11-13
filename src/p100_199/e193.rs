use std::time::Instant;

use itertools::Itertools;

use crate::utils::{
    FIArray::FIArrayU64, binary_search::binsearch, bit_array::BitArray,
    multiplicative_function_summation::mobius_sieve, sieve_of_pritchard::sift,
};

fn incl_excl(limit: u64, primes: &[u64]) -> u64 {
    let mut res = 0;
    for (i, &p) in primes.iter().enumerate() {
        let prod = p * p;
        if prod > limit {
            break;
        }
        res += limit / prod;
        res -= incl_excl(limit / prod, &primes[i + 1..]);
    }
    res
}
#[must_use]
fn count_squarefree(limit: u64) -> u64 {
    let first_primes = sift(limit.isqrt()).into_boxed_slice();
    limit - incl_excl(limit, &first_primes)
}

pub fn main() {
    const N: i64 = 1e16 as _;
    let start = Instant::now();
    let res = count_squarefree(N as _);
    let end = start.elapsed();
    println!("res = {res}, took {end:?}");
    let start = Instant::now();
    let res = count_squarefree_mobius_alt(N as _);
    let end = start.elapsed();
    println!("res = {res}, took {end:?}");
    let start = Instant::now();
    let res = count_squarefree_mobius(N as _);
    let end = start.elapsed();
    println!("res = {res}, took {end:?}");
    let start = Instant::now();
    let res = dirichlet_mul_based_opt(N as _);
    let end = start.elapsed();
    println!("res = {res}, took {end:?}");
    /* let start = Instant::now();
    let res = dirichlet_mul_based(N as _);
    let end = start.elapsed();
    println!("res = {res}, took {end:?}"); */
}
#[must_use]
fn count_squarefree_mobius(n: usize) -> usize {
    let m = mobius_sieve(n.isqrt() + 1).into_boxed_slice();
    (1..=n.isqrt())
        .map(|d| i64::from(m[d]) * (n / (d * d)) as i64)
        .sum::<i64>() as usize
}

#[must_use]
fn count_squarefree_mobius_alt(limit: usize) -> usize {
    let n = limit.isqrt() + 1;
    let mut sum = limit as i64;
    let mut res = vec![0i8; n];
    let mut composite = BitArray::zeroed(n);
    let mut primes = vec![];
    res[1] = 1;
    for i in 2..n {
        if !composite.get(i) {
            primes.push(i);
            res[i] = -1;
            sum -= (limit / (i * i)) as i64;
        }
        for &p in &primes {
            let ip = i * p;
            if ip >= n {
                break;
            }
            composite.set(ip);
            if i % p == 0 {
                res[ip] = 0;
                break;
            }
            res[ip] = -res[i];
            sum += i64::from(res[ip]) * (limit / (ip * ip)) as i64;
        }
    }
    sum as usize
}
// computation of zeta(2)/zeta(2s), but only for n/d^2 instead of for every n/d
fn dirichlet_mul_based_opt(x: u64) -> u64 {
    let xsqrt = x.isqrt();
    let primes = sift(xsqrt).into_boxed_slice();
    let xcbrt = (x as f64).cbrt() as u64;
    let mut keys = (1..=xcbrt)
        .map(|d| x / (d * d))
        .chain((1..=xcbrt).rev())
        .collect_vec();
    keys.reverse();
    assert!(keys.is_sorted());
    //dbg!(keys.len());
    keys.dedup();
    let keys = keys.into_boxed_slice();
    let rank = |e| binsearch(keys.len(), |i| keys[i], |v| v < e); //keys.partition_point(|&v| v < e);
    let mut s = keys.clone();
    for &p in &primes[1..] {
        for (i, &v) in keys.iter().enumerate().rev() {
            if v < p * p {
                break;
            }
            s[i] -= s[rank(v / (p * p))];
        }
    }
    s[rank(x)] - s[rank(x >> 2)]
}

// computation of zeta(2)/zeta(2s)
fn dirichlet_mul_based(x: u64) -> FIArrayU64 {
    let mut s = FIArrayU64::unit(x);
    let keys = FIArrayU64::keys(x).collect_vec().into_boxed_slice();
    let primes = sift(x.isqrt()).into_boxed_slice();

    for p in primes {
        for (i, &v) in keys.iter().enumerate().rev() {
            if v < p * p {
                break;
            }
            s.arr[i] -= s[v / (p * p)];
        }
    }
    s
}
