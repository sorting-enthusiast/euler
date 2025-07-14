#![warn(clippy::pedantic)]
#![allow(long_running_const_eval)]
#![allow(dead_code)]
#![allow(clippy::too_many_lines)]
#![allow(clippy::cast_possible_truncation)]
#![allow(clippy::inline_always)]
#![allow(clippy::cast_sign_loss)]
#![allow(clippy::cast_precision_loss)]

use std::time::Instant;

use crate::utils::{
    math::{modexp, prime_modinv},
    prime_sieves::{self, sieve_it},
    primecount::{
        lucy, lucy_alt, lucy_fastdivide, lucy_fastdivide_alt, lucy_strengthreduce,
        lucy_strengthreduce_alt, lucy_sum,
    },
};

mod e153;
mod e156;
mod e169;
mod e193;
mod e245;
mod e249;
mod e250;
mod e255;
mod e258;
mod e302;
mod e355;
mod e401;
mod e432;
mod e501;
mod e508;
mod e521;
mod e530;
mod e639;
mod e708;
mod e745;
mod e810;

pub mod longest_collatz_chain;
mod utils;
pub fn main() {
    //e153::main();
    //prime_sieves::main();
    const N: usize = 1e16 as _;
    let start = Instant::now();
    let count = lucy(N)[N];
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    let start = Instant::now();
    let count = lucy_alt(N)[N];
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    let start = Instant::now();
    let count = lucy_strengthreduce(N)[N];
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    let start = Instant::now();
    let count = lucy_strengthreduce_alt(N)[N];
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    let start = Instant::now();
    let count = lucy_fastdivide(N as _)[N as _];
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    let start = Instant::now();
    let count = lucy_fastdivide_alt(N as _)[N as _];
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    // a^(phi(n)-2) = a^-1 mod n
    //dbg!(prime_modinv::<{ 1e9 as u128 + 7 }>(42));
}
