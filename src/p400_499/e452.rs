use crate::utils::multiplicative_function_summation::general_divisor_summatory_alt;

// bruh
// sum of coefficients up to n in \zeta(s)^n
pub fn main() {
    const N: i64 = 1e9 as _;

    let start = std::time::Instant::now();
    let res = general_divisor_summatory_alt::<1_234_567_891>(N, N as _)[N];
    println!("res = {res}, took {:?}", start.elapsed());
}
