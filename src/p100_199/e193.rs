use std::time::Instant;

use crate::utils::{
    bit_array::BitArray, multiplicative_function_summation::mobius_sieve, sieve_of_pritchard::sift,
};

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
pub fn count_squarefree(limit: u128) -> u128 {
    let first_primes = sift((limit as f64).sqrt() as u64);
    limit - incl_excl(limit, 1, &first_primes)
}

pub fn main() {
    const N: i64 = 1 << 62;
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
}
#[must_use]
fn count_squarefree_mobius(n: usize) -> usize {
    let m = mobius_sieve(n.isqrt() + 1);
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
