use crate::utils::{math, primes::prime_sieves::sift};
use std::time::Instant;
const N: u64 = 1e18 as u64;

pub fn main() {
    //302:
    //achilles number: every prime factor has multiplicity > 1 and gcd(multiplicities)=1
    //strong achilles number: both num and totient(num) are achilles numbers
    let start = Instant::now();
    //note: the multiplicity of the largest prime factor has to be at least 3 for num to be a strong achilles number,
    //and it has to have at least 2 prime factors -> 2^2 * p^3 <= 1e18 ->
    //p <= 1e6 / 2^(2/3)
    let mut achilles_candidate_factorization = vec![];
    let primes = sift((N as f64 / 4.0).cbrt() as u64);
    let end1 = start.elapsed();
    let ans = dfs(
        0,
        1,
        0,
        1,
        0,
        false,
        &primes,
        &mut achilles_candidate_factorization,
    );
    let end = start.elapsed();
    dbg!(end1);
    dbg!(end);
    dbg!(ans);
}
fn check_achilles_totient(
    phi: u64,
    primes: &[u64],
    candidate_factorization: &[(u64, u64)],
) -> bool {
    let mut gcd_multiplicities = 0;
    let mut x = phi;
    for &(prime_factor, multiplicity_minus_1) in candidate_factorization {
        let mut new_multiplicity = multiplicity_minus_1;
        while x % prime_factor == 0 {
            x /= prime_factor;
            new_multiplicity += 1;
        }
        if new_multiplicity == 1 {
            //if not powerful, not achilles
            return false;
        }
        gcd_multiplicities = math::gcd(gcd_multiplicities, new_multiplicity);
    }
    for &prime in primes {
        if prime * prime > x {
            break;
        }
        let mut multiplicity = 0;
        while x % prime == 0 {
            multiplicity += 1;
            x /= prime;
        }
        if multiplicity == 1 {
            //if not powerful, not achilles
            return false;
        }
        gcd_multiplicities = math::gcd(gcd_multiplicities, multiplicity);
    }

    gcd_multiplicities == 1 && x == 1 //is not perfect power
}
fn dfs(
    depth: usize,
    achilles_candidate: u64,
    gcd_multiplicities: u64,
    phi: u64,
    multiplicity_largest_prime_factor: u64,
    is_candidate: bool,
    primes: &[u64],
    candidate_factorization: &mut Vec<(u64, u64)>,
) -> u64 {
    let mut ans = 0;
    if gcd_multiplicities == 1 && multiplicity_largest_prime_factor > 2 && is_candidate {
        ans += u64::from(check_achilles_totient(phi, primes, candidate_factorization));
    }
    if depth > primes.len() - 1 {
        return ans;
    }
    let prime = primes[depth];
    let mut prime_power = prime * prime;
    let mut power = 2;
    if achilles_candidate > N / prime_power {
        return ans;
    }
    while achilles_candidate <= N / prime_power {
        candidate_factorization.push((prime, power - 1));
        ans += dfs(
            depth + 1,
            achilles_candidate * prime_power,
            math::gcd(gcd_multiplicities, power),
            phi * (prime - 1),
            power,
            true,
            primes,
            candidate_factorization,
        );
        candidate_factorization.pop();

        if prime > N / prime_power {
            break;
        }
        prime_power *= prime;
        power += 1;
    }
    ans + dfs(
        depth + 1,
        achilles_candidate,
        gcd_multiplicities,
        phi,
        multiplicity_largest_prime_factor,
        false,
        primes,
        candidate_factorization,
    )
}
