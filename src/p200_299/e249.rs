use crate::utils::primes::prime_sieves::erat_wheel210;
use std::time::Instant;
// calculate product of (1+x^p) for all p < 5000, i.e. the generating function of all subsets of sum i, sum coefficients of prime powers
pub fn main() {
    const MOD: usize = 1e16 as usize;

    let start = Instant::now();

    let n = 1_548_136; // sum of primes below 5000
    let primes = erat_wheel210(n);
    let mut dp = vec![0; n + 1];
    dp[0] = 1;

    let mut max = 0;
    for &p in &primes[..669] {
        max += p; // largest possible sum at current iteration
        for j in (p..=max).rev() {
            dp[j] += dp[j - p];
            if dp[j] >= MOD {
                dp[j] -= MOD;
            }
        }
    }
    let mut sum = 0;
    for p in primes {
        sum += dp[p];
        if sum >= MOD {
            sum -= MOD;
        }
    }
    println!("{:?}", start.elapsed());
    println!("{sum}");
}
