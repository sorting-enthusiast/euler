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

use crate::{
    p300_399::e362::mult,
    utils::{
        FIArray::FIArray,
        multiplicative_function_summation::mertens,
        primes::{
            log_zeta::log_zeta,
            primecount::{lucy_fenwick, mertens_min25},
        },
    },
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
    p300_399::e362::main();
    utils::primes::primecount::main();
    const N: i64 = 1e8 as _;
    assert_eq!(lucy_fenwick(N as _), log_zeta(N as _));
    let start = std::time::Instant::now();
    let s1 = mertens(N);
    let end = start.elapsed();
    dbg!(end, s1[N]);
    let start = std::time::Instant::now();
    let s1 = mertens_min25(N);
    let end = start.elapsed();
    dbg!(end, s1[N]);
    let pi = log_zeta(N as _);
    let mut pi2 = FIArray::new(N as _);
    for p in 2..=pi.isqrt {
        if pi.arr[p - 1] == pi.arr[p - 2] {
            continue;
        }
        pi2[p * p] += 1;
    }
    for i in 1..pi2.arr.len() {
        pi2.arr[i] += pi2.arr[i - 1];
    }
    let pi_squared = mult(&pi, &pi);
    for i in 0..pi2.arr.len() {
        pi2.arr[i] += pi_squared.arr[i];
    }
    for i in 0..pi2.arr.len() {
        pi2.arr[i] >>= 1;
    }
    dbg!(pi2[N as _]);
    p100_199::e187::main();
    //utils::primes::prime_sieves::main();
    println!("Finished running at: {} ", Local::now().time());
}
