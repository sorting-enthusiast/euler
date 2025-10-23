use itertools::Itertools;

use crate::utils::prime_sieves::sift;

// Hilbert numbers are numbers that are products of any number of primes equivalent to 1 mod 4,
// and an even number of primes equivalent to 3 mod 4
fn incl_excl(limit: u128, acc: u128, primes: &[u64]) -> u128 {
    let mut res = 0;
    for i in 0..primes.len() {
        let p = u128::from(primes[i]);
        let prod = acc * (p * p);
        assert!(prod & 3 == 1);
        if prod > limit {
            break;
        }
        res += (limit / prod - 1) >> 2;
        res -= incl_excl(limit, prod, &primes[i + 1..]);
    }
    res
}
#[must_use]
pub fn count_squarefree(limit: u128) -> u128 {
    let mut first_primes = sift((limit as f64).sqrt() as u64);
    dbg!(first_primes.len());
    let mut primes = first_primes
        .iter()
        .copied()
        .filter(|&p| p & 3 == 1)
        .collect_vec();
    first_primes.retain(|p| p & 3 == 3);

    for (i, &p1) in first_primes.iter().enumerate() {
        //dbg!(i, p1);
        for &p2 in &first_primes[i..] {
            if (p1 * p2).pow(2) as u128 <= limit {
                primes.push(p1 * p2);
            }
        }
    }
    primes.sort();
    dbg!(primes.len());
    dbg!(&primes[..15]);
    ((limit - 1) >> 2) - incl_excl(limit, 1, &primes)
}

pub fn main() {
    const N: i64 = 1e16 as _;
    let start = std::time::Instant::now();
    let res = count_squarefree(N as _);
    let end = start.elapsed();
    println!("res = {res}, took {end:?}");
}
