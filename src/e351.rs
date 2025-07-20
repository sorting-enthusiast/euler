use crate::utils::multiplicative_function_summation::totient_sum_alt;

pub fn main() {
    const N: i64 = 1e8 as i64;
    // original
    /* let start = std::time::Instant::now();
    let mut sum = (N - 1) * N;
    let mobius = mobius_sieve(N as usize + 1);
    for d in 1..=N {
        let nd = N / d;
        sum -= mobius[d as usize] * (nd - 1) * nd;
    }
    sum >>= 1;
    sum += N - 1;
    sum *= 6;
    let end = start.elapsed();
    println!("res = {sum}, took {end:?}"); */
    let start = std::time::Instant::now();
    let totient_sum = totient_sum_alt::<{ i64::MAX >> 1 }>(N)[N];
    let sum = 6 * (((N * (N + 1)) >> 1) - totient_sum);
    let end = start.elapsed();
    println!("res = {sum}, took {end:?}");
}
