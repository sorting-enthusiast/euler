#![warn(clippy::pedantic)]
#![allow(long_running_const_eval)]
#![allow(dead_code)]
#![allow(clippy::too_many_lines)]
#![allow(clippy::cast_possible_truncation)]
#![allow(clippy::inline_always)]
#![allow(clippy::cast_sign_loss)]
#![allow(clippy::cast_precision_loss)]
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]
#![allow(clippy::large_stack_arrays)]
use chrono::Local;

use crate::utils::{
    FIArray::{DirichletFenwickU128, FIArray, FIArrayU128},
    fenwick::FenwickTreeU128,
    multiplicative_function_summation::{
        count_squarefree, inverse_pseudo_euler_transform, mertens, mertens_slow,
        pseudo_euler_transform, sum_n_u128,
    },
    primes::primecount::{lucy_fenwick, mertens_min25},
};

pub mod p0_99;
pub mod p100_199;
pub mod p200_299;
pub mod p300_399;
pub mod p400_499;
pub mod p500_599;
pub mod p600_699;
pub mod p700_799;
pub mod p800_899;
pub mod p900_999;
pub mod utils;

// digital root of n is just n mod 9 if n mod 9 != 0, otherwise 9
const fn is_target_little_endian() -> bool {
    u16::from_ne_bytes([1, 0]) == 1
}
// TODO: understand convex hull based lattice point counting, optimize dirichlet mul
pub fn main() {
    const { assert!(is_target_little_endian()) }; // some code relies on this
    println!("Started running at: {} ", Local::now().time());
    //p500_599::e580::main();
    //p800_899::e890::main();
    /* for i in 2..=4 {
        let n = 7u128.pow(1 << i) - 1;
        dbg!(n);
        let mut v =
            3 * sum_n_u128::<0>(n / 3) + 5 * sum_n_u128::<0>(n / 5) - 15 * sum_n_u128::<0>(n / 15);
        let mut sum = 0;
        while v != 0 {
            sum += v % 10;
            v /= 10;
        }
        print!("{i}:{sum}, ");
    }
    println!(); */
    //p300_399::e362::main();
    utils::primes::primecount::main();
    const N: i64 = 1e11 as _;
    let start = std::time::Instant::now();
    let s1 = mertens(N);
    let end = start.elapsed();
    dbg!(end, s1[N]);
    let start = std::time::Instant::now();
    let s1 = mertens_min25(N);
    let end = start.elapsed();
    dbg!(end, s1[N]);
    //utils::primes::prime_sieves::main();
    println!("Finished running at: {} ", Local::now().time());
}
