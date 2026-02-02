use itertools::Itertools;

use crate::utils::multiplicative_function_summation::divisor_sieve;

pub fn main() {
    // only 1 and 5 are invertible mod 6,
    // n'th die is 1 iff it is a product of an even number of primes with multiplicity = 4 mod 6, and any number of primes with multiplicity = 0 mod 6
    const N: usize = 1e8 as _;
    let d = divisor_sieve(N + 1);
    for n in 1..=N {
        if d[n] % 6 == 1 {
            dbg!(n, d[n]);
        }
    }
}
