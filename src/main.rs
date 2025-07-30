#![warn(clippy::pedantic)]
#![allow(long_running_const_eval)]
#![allow(dead_code)]
#![allow(clippy::too_many_lines)]
#![allow(clippy::cast_possible_truncation)]
#![allow(clippy::inline_always)]
#![allow(clippy::cast_sign_loss)]
#![allow(clippy::cast_precision_loss)]

use crate::utils::{
    FIArray::FIArrayI64,
    multiplicative_function_summation::{
        dirichlet_mul, general_divisor_summatory, mertens, sum_n_i64, totient_sum,
    },
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
mod e668;
mod e708;
mod e745;
mod e810;
mod e955;
mod longest_collatz_chain;
mod utils;

pub fn main() {
    //e668::main();
    const N: i64 = 1e4 as _;
    const MOD: i64 = 1e9 as i64 + 7;

    let start = std::time::Instant::now();
    let divi = general_divisor_summatory(N, 63);
    let end = start.elapsed();
    println!("res = {}, took {end:?}", divi[N]);
}
