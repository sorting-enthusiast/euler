use crate::sieve_of_pritchard::sift;

fn incl_excl(limit: u128, acc: u128, primes: &[u64]) -> u128 {
    let mut res = 0;
    for i in 0..primes.len() {
        let p = primes[i] as u128;
        let prod = acc * (p * p);
        if prod > limit {
            break;
        }
        res += limit / prod;
        res -= incl_excl(limit, prod, &primes[i + 1..]);
    }
    res
}
pub fn count_squarefree(limit: u128) -> u128 {
    let first_primes = sift((limit as f64).sqrt() as u64);
    limit - incl_excl(limit, 1, &first_primes)
}
