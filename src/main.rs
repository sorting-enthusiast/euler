#![warn(clippy::pedantic)]
#![allow(long_running_const_eval)]
#![allow(dead_code)]
#![allow(clippy::too_many_lines)]
#![allow(clippy::cast_possible_truncation)]
#![allow(clippy::inline_always)]
#![allow(clippy::cast_sign_loss)]
#![allow(clippy::cast_precision_loss)]

use crate::utils::primecount;

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
mod e484;
mod e501;
mod e508;
mod e521;
mod e530;
mod e639;
mod e708;
mod e745;
mod e810;
mod e955;

pub mod longest_collatz_chain;
mod utils;
pub fn main() {
    e508::main();
    // a^(phi(n)-2) = a^-1 mod n
    //dbg!(prime_modinv::<{ 1e9 as u128 + 7 }>(42));
}
