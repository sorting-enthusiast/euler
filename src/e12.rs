use std::time::Instant;

use crate::utils::multiplicative_function_summation::divisor_sieve;

const N: usize = 15_000;
fn find() -> Option<usize> {
    let divs = divisor_sieve(N);
    for n in 20..=N {
        if n & 1 == 0 {
            if divs[n >> 1] * divs[n + 1] >= 500 {
                return Some((n >> 1) * (n + 1));
            }
        } else if divs[(n + 1) >> 1] * divs[n] >= 500 {
            return Some(((n + 1) >> 1) * n);
        }
    }
    None
}
pub fn main() {
    let start = Instant::now();
    let res = find();
    let end = start.elapsed();
    println!("res = {res:?}, took {end:?}");
}
