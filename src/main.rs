#![warn(clippy::pedantic)]
#![allow(long_running_const_eval)]
#![allow(dead_code)]
#![allow(clippy::too_many_lines)]
#![allow(clippy::cast_possible_truncation)]
#![allow(clippy::inline_always)]
#![allow(clippy::cast_sign_loss)]
#![allow(clippy::cast_precision_loss)]

use crate::utils::multiplicative_function_summation::{
    mertens, mertens_slow, totient_sum, totient_sum_single,
};

mod e107;
mod e153;
mod e156;
mod e169;
mod e175;
mod e193;
mod e214;
mod e233;
mod e245;
mod e249;
mod e250;
mod e255;
mod e258;
mod e302;
mod e351;
mod e355;
mod e401;
mod e432;
mod e448;
mod e484;
mod e501;
mod e508;
mod e512;
mod e521;
mod e530;
mod e625;
mod e639;
mod e708;
mod e745;
mod e810;
mod e955;
mod longest_collatz_chain;
mod utils;

pub fn main() {
    const N: i64 = 1e13 as i64;

    let start = std::time::Instant::now();
    let m1 = totient_sum_single::<{ 1e9 as _ }>(N);
    let end = start.elapsed();
    println!("res = {m1}, took {end:?}");

    let start = std::time::Instant::now();
    let m2 = totient_sum::<{ 1e9 as _ }>(N)[N];
    let end = start.elapsed();
    println!("res = {m2}, took {end:?}");

    /* let p = sift(1e6 as _);
    dbg!(p.len());
    dbg!(
        p.len()
            - p.iter()
                .filter(|&&p| powmod(5, (p - 1) >> 1, p) == p - 1)
                .count()
    );
    dbg!(
        p.iter()
            .filter(|&&p| p != 2 && powmod(5, (p - 1) >> 1, p) == 1)
            .take(10)
            .collect_vec()
    ); */
    // a^(phi(n)-2) = a^-1 mod n
    //dbg!(prime_modinv::<{ 1e9 as u128 + 7 }>(42));
}
