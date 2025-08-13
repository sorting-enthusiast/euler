#![warn(clippy::pedantic)]
#![allow(long_running_const_eval)]
#![allow(dead_code)]
#![allow(clippy::too_many_lines)]
#![allow(clippy::cast_possible_truncation)]
#![allow(clippy::inline_always)]
#![allow(clippy::cast_sign_loss)]
#![allow(clippy::cast_precision_loss)]
use std::i64;

use crate::utils::{
    FIArray::FIArrayI64,
    multiplicative_function_summation::{
        dirichlet_div_i64, dirichlet_mulmod_i64, mertens, sum_n_i64, totient_sieve, totient_sum,
        totient_sum_single,
    },
    prime_sieves::sieve_it,
    primecount,
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
mod e273;
mod e302;
mod e351;
mod e355;
mod e401;
mod e415;
mod e432;
mod e439;
mod e448;
mod e484;
mod e501;
mod e508;
mod e512;
mod e521;
mod e530;
mod e606;
mod e625;
mod e639;
mod e668;
mod e708;
mod e745;
mod e810;
mod e955;
mod longest_collatz_chain;
mod utils;

const MOD: i64 = 1e9 as i64 + 7;
const fn powmod(mut x: i64, mut exp: i64) -> i64 {
    if exp == 0 {
        return 1;
    }
    let mut r = 1;
    x %= MOD;
    while exp > 1 {
        if exp & 1 == 1 {
            r = (r * x) % MOD;
        }
        x = (x * x) % MOD;
        exp >>= 1;
    }
    (r * x) % MOD
}
const fn modinv(x: i64) -> i64 {
    powmod(x, MOD - 2)
}

pub fn main() {
    dbg!(totient_sum_single::<{ i64::MAX >> 1 }>(16000));
    //e415::main();
}
