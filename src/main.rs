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
    FIArray::FIArray,
    math::iroot,
    multiplicative_function_summation::{
        count_squarefree, inverse_pseudo_euler_transform, inverse_pseudo_euler_transform_fraction,
        mertens, pseudo_euler_transform, pseudo_euler_transform_fraction,
    },
    primes::{
        log_zeta::log_zeta,
        primecount::{lucy_fenwick, mertens_min25},
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
    //p300_399::e379::main();
    const N: i64 = 1e11 as _;
    let fsf = log_zeta(N as _);
    let start = std::time::Instant::now();
    let s2 = pseudo_euler_transform(&fsf);
    let end = start.elapsed();
    dbg!(end, s2[N as _]);
    let start = std::time::Instant::now();
    let s1 = pseudo_euler_transform_fraction(&fsf);
    let end = start.elapsed();
    dbg!(end, s1[N as _]);
    assert_eq!(s1, s2);
    println!("hello");
    assert_eq!(lucy_fenwick(N as _), log_zeta(N as _));
    let start = std::time::Instant::now();
    let s1 = mertens(N);
    let end = start.elapsed();
    dbg!(end, s1[N]);
    let start = std::time::Instant::now();
    let s1 = mertens_min25(N);
    let end = start.elapsed();
    dbg!(end, s1[N]);
    p300_399::e362::main();
    utils::primes::primecount::main();
    //utils::primes::prime_sieves::main();
    println!("Finished running at: {} ", Local::now().time());
}
