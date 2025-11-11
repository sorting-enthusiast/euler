#![warn(clippy::pedantic)]
#![allow(long_running_const_eval)]
#![allow(dead_code)]
#![allow(clippy::too_many_lines)]
#![allow(clippy::cast_possible_truncation)]
#![allow(clippy::inline_always)]
#![allow(clippy::cast_sign_loss)]
#![allow(clippy::cast_precision_loss)]

use chrono::Local;

use crate::utils::factorint::factor;

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
pub fn main() {
    println!("Started running at: {} ", Local::now().time());
    //p200_299::e269::main();
    //p500_599::e548::main();
    //utils::primecount::main();

    println!("Finished running at: {} ", Local::now().time());
}
