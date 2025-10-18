use crate::utils::multiplicative_function_summation::totient_sum;

pub fn main() {
    let start = std::time::Instant::now();
    let mut n = 5e8 as i64;
    let sums = totient_sum::<{ i64::MAX >> 1 }>(n);
    let mut res = sums[n];
    n >>= 1;
    while n != 0 {
        res -= sums[n];
        n >>= 1;
    }
    let end = start.elapsed();
    println!("{res}, took {end:?}");
}
