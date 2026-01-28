use crate::utils::multiplicative_function_summation::totient_sum_single;

pub fn main() {
    const N: usize = 1e8 as _;
    let start = std::time::Instant::now();
    let totient_sum = totient_sum_single::<0>(N);
    let sum = 6 * (((N * (N + 1)) >> 1) - totient_sum as usize);
    let end = start.elapsed();
    println!("res = {sum}, took {end:?}");
}
