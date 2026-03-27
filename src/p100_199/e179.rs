use crate::utils::multiplicative_function_summation::divisor_sieve;

pub fn main() {
    dbg!(
        divisor_sieve(1e7 as usize + 1)
            .windows(2)
            .filter(|w| w[0] == w[1])
            .count()
    );
}
