use crate::utils::primes::{log_zeta::log_zeta, prime_sieves::sift};

pub fn main() {
    let primes = sift(1_991_010);
    let mut count = 0;
    let mut max_len = 40;
    let mut best_a = -1;
    let mut best_b = 41;
    for &p in primes.iter().take_while(|&&p| p <= 1000) {
        let b = p as i64;
        let lo = (-40 + (-b + 2) / 40) | 1;
        for a in (lo..1000).step_by(2) {
            if primes.binary_search(&((1 + a + b) as u64)).is_ok()
                && primes
                    .binary_search(&((40 * 40 + a * 40 + b) as u64))
                    .is_ok()
                && (2..40).all(|n| primes.binary_search(&((n * n + n * a + b) as u64)).is_ok())
            {
                let extra = (41..b)
                    .take_while(|n| primes.binary_search(&((n * n + n * a + b) as u64)).is_ok())
                    .count();
                if extra > max_len - 40 {
                    max_len = 40 + extra;
                    best_a = a;
                    best_b = b;
                }
                println!("n^2 + {a}n + {b}");
                count += 1;
            }
        }
    }
    dbg!(best_a, best_b, best_a * best_b, max_len);
    dbg!(count);
    dbg!(log_zeta(1000)[1000]);
}
