#![warn(clippy::pedantic)]
#![allow(long_running_const_eval)]
#![allow(dead_code)]
#![allow(clippy::too_many_lines)]
#![allow(clippy::cast_possible_truncation)]
#![allow(clippy::inline_always)]
#![allow(clippy::cast_sign_loss)]
#![allow(clippy::cast_precision_loss)]

use std::u64;

use chrono::Local;
use itertools::Itertools;

use crate::utils::{factorint::factor, multiplicative_function_summation::sum_n_u64};

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
    //p900_999::e967::main();
    //utils::primecount::main();
    dbg!(factor((1 << 15) + 1).iter().sorted_unstable());
    /*let n = 2e4 as i64;
    let prime_pi = utils::primecount::lucy(n as usize);
    let sqfree_1 = count_squarefree(n.isqrt() as _) as usize;
    let sqfree_2 = count_squarefree((n as f64).cbrt() as _) as usize;
    dbg!(dbg!(PowerfulExt::<_, 2>::new(n, |_, _| 1).count()) - dbg!(sqfree_1) - dbg!(sqfree_2) + 1); */
    println!("Finished running at: {} ", Local::now().time());
}
