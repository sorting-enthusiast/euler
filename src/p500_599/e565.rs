use itertools::Itertools;

use crate::utils::primes::wheel_sieve;

// primitive root mod 2017 = 5
// find all prime powers p^e s.t. \frac{p^{e + 1} - 1}{p - 1} = 0 mod 2017
pub fn main() {
    let primes = wheel_sieve(1e7 as _);
    dbg!(primes.iter().filter(|&&p| p % 2017 == 1).count());
}
