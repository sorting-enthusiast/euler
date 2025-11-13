use itertools::Itertools;

use crate::utils::{FIArray::FIArrayU64, binary_search::binsearch, prime_sieves::sift};

// Hilbert numbers are numbers that are products of any number of primes equivalent to 1 mod 4,
// and an even number of primes equivalent to 3 mod 4

// a hilbert number is h-squarefree if it is either squarefree, or if only 1 of its prime factors = 3 mod 4 appears with multiplicity 1< <4

pub fn main() {
    const N: u64 = 1e16 as _;
    println!("res = {}", 2_327_213_148_095_366_u64);
    let start = std::time::Instant::now();
    let res = solve(N as _);
    let end = start.elapsed();
    println!("res = {res}, took {end:?}");
    let start = std::time::Instant::now();
    let res = lucy_based(N as _);
    let end = start.elapsed();
    println!("res = {res}, took {end:?}");
    let start = std::time::Instant::now();
    let res = lucy_based_2(N as _);
    let end = start.elapsed();
    println!("res = {res}, took {end:?}");
}

fn solve(x: u64) -> u64 {
    let primes = sift(x.isqrt());
    let squarefree = |x| ((x + 3) >> 2) - incl_excl(x, &primes[1..]);
    primes
        .iter()
        .filter(|&&p| p & 3 == 3)
        .fold(squarefree(x), |acc, p| acc + squarefree(x / (p * p)))
}
fn incl_excl(limit: u64, primes: &[u64]) -> u64 {
    let mut res = 0;
    for (i, &p) in primes.iter().enumerate() {
        let prod = p * p;
        if prod > limit {
            break;
        }
        res += ((limit / prod) + 3) >> 2;
        res -= incl_excl(limit / prod, &primes[i + 1..]);
    }
    res
}

fn lucy_based_2(x: u64) -> u64 {
    let xsqrt = x.isqrt();
    let primes = sift(xsqrt);
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
    //dbg!(keys.len());
    let rank = |e| binsearch(keys.len(), |i| keys[i], |v| v < e); //keys.partition_point(|&v| v < e);
    let mut squarefree = vec![0; keys.len()].into_boxed_slice();
    for (i, &v) in keys.iter().enumerate() {
        squarefree[i] = (v + 3) >> 2;
    }
    for &p in &primes[1..] {
        for (i, &v) in keys.iter().enumerate().rev() {
            if v < p * p {
                break;
            }
            squarefree[i] -= squarefree[rank(v / (p * p))];
        }
    }
    primes
        .iter()
        .filter(|&&p| p & 3 == 3)
        .fold(squarefree[rank(x)], |acc, p| {
            acc + squarefree[rank(x / (p * p))]
        })
}

fn lucy_based(x: u64) -> u64 {
    let xsqrt = x.isqrt();
    let primes = sift(xsqrt);
    let xcbrt = (x as f64).cbrt() as u64;
    let mut keys = (1..=xcbrt)
        .map(|d| x / (d * d))
        .chain((1..=xcbrt).rev())
        .collect_vec();
    keys.reverse();
    assert!(keys.is_sorted());
    keys.dedup();
    let keys = keys.into_boxed_slice();
    let mut squarefree = FIArrayU64::new(x);
    for &v in &keys {
        squarefree[v] = (v + 3) >> 2;
    }
    for &p in &primes[1..] {
        for &v in keys.iter().rev() {
            if v < p * p {
                break;
            }
            squarefree[v] -= squarefree[v / (p * p)];
        }
    }
    primes
        .iter()
        .filter(|&&p| p & 3 == 3)
        .fold(squarefree[x], |acc, p| acc + squarefree[x / (p * p)])
}
