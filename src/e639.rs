use std::time::Instant;

use crate::utils::{
    multiplicative_function_summation::{divisor_summatory, sum_n},
    powerful_numbers::PowerfulExt,
};
const N: i64 = 1e4 as i64;
const MOD: i64 = 1e9 as i64 + 7;
// powerful number trick?
// 50 lucy/unlucy iterations could work, but would be pretty ugly. could precompute coefficients in prefix sums of g, so it wouldn't be too bad, but nonetheless.

pub fn main() {
    let start = Instant::now();

    let h = |p, e| {
        unsafe { core::hint::assert_unchecked(e > 1) };
        p * (1 - p)
    };
    let mut sum = 0;
    let mut count = 0;
    let mut powerful = PowerfulExt::<_, MOD>::new(N, h);
    for (n, hn) in powerful.by_ref().filter(|&(_, hn)| hn != 0) {
        sum += MOD + (hn * sum_n::<MOD>(N / n)) % MOD;
        sum %= MOD;
        count += 1;
    }
    let end = start.elapsed();
    println!("res = {sum}, took {end:?}");
    dbg!((powerful.max_len, count));
}
