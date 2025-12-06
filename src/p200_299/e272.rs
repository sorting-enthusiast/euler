use itertools::Itertools;

use crate::utils::{FIArray::FIArrayU64, primes::prime_sieves::sift};

// f_3(n) counts # of cube roots of unity mod n, including 1
// f_3(n) is multiplicative
// if p = 3k+1: f_3(p^e) = 3;
// if p = 3k+2: f_3(p^e) = 1;
// f_3(3) = 1, f_3(3^e), e>=2 = 3
// sum all numbers with f_3(n) = 243 = 3^5
const N: u64 = 1e11 as _;

fn dfs(f: &impl Fn(u64) -> u64, depth: u8, acc: u64, lim: u64, primes: &[u64]) -> u64 {
    if depth == 0 {
        return f(acc);
    }
    let mut sum = 0;
    for (i, &p) in primes.iter().enumerate() {
        if p > lim {
            break;
        }
        let mut new_acc = acc * p;
        let mut new_lim = lim;
        loop {
            new_lim /= p;
            sum += dfs(f, depth - 1, new_acc, new_lim, &primes[i + 1..]);
            new_acc *= p;
            if p.pow(u32::from(depth) - 1) > new_lim {
                break;
            }
        }
    }
    sum
}
pub fn main() {
    let start = std::time::Instant::now();
    let primes = sift(N / (482_391 / 31))
        .into_iter()
        .filter(|p| p % 3 == 1)
        .collect_vec();
    dbg!(start.elapsed());
    let mut s = FIArrayU64::id::<{ u64::MAX >> 1 }>(N);
    let keys = FIArrayU64::keys(N)
        .take_while(|&v| v <= N / 482_391)
        .collect_vec()
        .into_boxed_slice();
    let lim = primes.partition_point(|&p| p <= N / 482_391);
    for &p in &primes[..lim] {
        for (i, &v) in keys.iter().enumerate().rev() {
            if p > v {
                break;
            }
            s.arr[i] -= p * s[v / p];
        }
    }
    dbg!(start.elapsed());

    let res = dfs(
        &|prod| prod * s[N / prod] - 9 * prod * s[N / (9 * prod)],
        5,
        1,
        N,
        &primes,
    ) + dfs(&|prod| prod * s[N / prod], 4, 9, N / 9, &primes);
    println!("res = {res}, took {:?}", start.elapsed());
}
