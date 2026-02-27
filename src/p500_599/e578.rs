use crate::utils::{FIArray::FIArray, primes::log_zeta::log_zeta_2};
const N: usize = 1e13 as _;
const SQRT: usize = N.isqrt();

pub fn main() {
    let start = std::time::Instant::now();
    let pi = log_zeta_2(N);
    let mut primes = vec![];
    for p in 2..=pi.isqrt {
        if pi.arr[p - 1] != pi.arr[p - 2] {
            primes.push(p);
        }
    }
    dbg!(start.elapsed());
    let res = rec(N, N.ilog2() as _, &primes, &pi);
    println!("res = {res}, took {:?}", start.elapsed());
}
fn rec(lim: usize, e: u8, primes: &[usize], pi: &FIArray) -> usize {
    if primes.is_empty() {
        return 1;
    }
    let rt = lim.isqrt();
    let mut ret = 1 + pi[lim.max(primes[0] - 1)] - pi[rt.max(primes[0] - 1)];
    for (i, &p) in primes.iter().enumerate() {
        if rt < p {
            break;
        }
        let mut new_lim = lim / p;
        ret += sqf(new_lim, &primes[i + 1..], pi);

        for ex in 2..=e {
            new_lim /= p;
            if new_lim == 0 {
                break;
            }
            ret += rec(new_lim, ex, &primes[i + 1..], pi);
        }
    }
    ret
}
fn sqf(lim: usize, primes: &[usize], pi: &FIArray) -> usize {
    if primes.is_empty() {
        return 1;
    }
    let rt = lim.isqrt();
    let mut ret = 1 + pi[lim.max(primes[0] - 1)] - pi[rt.max(primes[0] - 1)];
    for (i, &p) in primes.iter().enumerate() {
        if rt < p {
            break;
        }
        ret += sqf(lim / p, &primes[i + 1..], pi);
    }
    ret
}
