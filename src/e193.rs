use std::time::Instant;

use crate::utils::{powerful_numbers::PowerfulExt, sieve_of_pritchard::sift};

fn incl_excl(limit: u128, acc: u128, primes: &[u64]) -> u128 {
    let mut res = 0;
    for i in 0..primes.len() {
        let p = u128::from(primes[i]);
        let prod = acc * (p * p);
        if prod > limit {
            break;
        }
        res += limit / prod;
        res -= incl_excl(limit, prod, &primes[i + 1..]);
    }
    res
}
#[must_use]
fn count_squarefree(limit: u128) -> u128 {
    let first_primes = sift((limit as f64).sqrt() as u64);
    limit - incl_excl(limit, 1, &first_primes)
}

pub fn main() {
    const N: i64 = 5e15 as _;
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
    let mut sum = 0;
    let h = |p, e| {
        if e == 2 { -1 } else { 0 }
    };
    for (n, hn) in PowerfulExt::<_, { i64::MAX }>::new(N, h).filter(|&(_, hn)| hn != 0) {
        sum += hn * (N / n);
    }
    let end = start.elapsed();
    println!("res = {sum}, took {end:?}");
}
// done stupid: can probably use hyperbola method or powerful numbers to speed it up
fn mobius_until(n: usize) -> Vec<i8> {
    let mut m = vec![2; n + 1];
    m[1] = 1;
    for i in 1..=n >> 1 {
        m[i << 1] = -((i & 1) as i8);
    }
    for p in 3..=n >> 1 {
        if m[p] != 2 {
            continue;
        }
        m[p] = -1;

        let mut multiple = p << 1;
        let mut i = 2;
        while multiple <= n {
            m[multiple] = if i == p { 0 } else { -(m[multiple].min(1)) };
            i = if i == p { 0 } else { i } + 1;
            multiple += p;
        }
    }
    m
}

#[must_use]
fn count_squarefree_mobius(n: usize) -> usize {
    let m = mobius_until(n.isqrt() + 1);
    (1..=n.isqrt())
        .map(|d| { if m[d] == 2 { -1 } else { m[d] } } as isize * (n / (d * d)) as isize)
        .sum::<isize>() as usize
}

#[must_use]
fn count_squarefree_mobius_alt(limit: usize) -> usize {
    let n = limit.isqrt() + 1;
    let mut mu = vec![2; n + 1];
    mu[1] = 1;
    let mut sum = limit - (limit >> 2);
    for i in 1..=n >> 1 {
        mu[i << 1] = -((i & 1) as i8);
    }
    for p in 3..=n {
        if mu[p] == 2 {
            // prime
            sum -= limit / (p * p);
            if p * 2 > n {
                continue;
            }
            let mut multiple = p << 1;
            let mut i = 2;
            while multiple <= n {
                mu[multiple] = if i == p { 0 } else { -(mu[multiple].min(1)) };
                i = if i == p { 0 } else { i } + 1;
                multiple += p;
            }
        } else {
            // composite
            sum = (sum as isize + mu[p] as isize * (limit / (p * p)) as isize) as usize;
        }
    }
    sum
}
