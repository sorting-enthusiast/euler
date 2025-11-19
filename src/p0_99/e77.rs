use crate::utils::primes::prime_sieves::sift;
//TODO
pub fn main() {
    const MAX: usize = 1e5 as _;

    let mut dp = vec![0; MAX + 1];
    let primes = sift(MAX as _);
    for &p in &primes {
        let p = p as usize;
        dp[p] = 1;
    }
    for n in 4..=MAX {
        for &p in &primes {
            let p = p as usize;
            if p > n - p {
                break;
            }
            dp[n] += dp[n - p];
        }
        dbg!((n, dp[n]));
        if dp[n] >= 5000 {
            break;
        }
    }
}
