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

use crate::utils::{FIArray::FIArray, fast_divisor_sums::sum_floors_fast};

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
// TODO: impl FIArrayFenwick as a struct of its own, rather than copy/paste the functionality every time
pub fn main() {
    const { assert!(is_target_little_endian()) }; // some code relies on this
    println!("Started running at: {} ", Local::now().time());
    /*let start = std::time::Instant::now();
    let d = sum_floors_fast(1e18 as _);
    let end = start.elapsed();
    dbg!(d, end);
    let start = std::time::Instant::now();
    let d = divisor_summatory_i64(1e18 as _);
    let end = start.elapsed();
    dbg!(d, end);*/
    //p200_299::e240::main();
    //p500_599::e580::main();
    p300_399::e379::main(); //132314136838185

    /* let n = 1e7 as i64;
    let mut s = FIArrayI64::eps(n);
    let keys = FIArrayI64::keys(n).collect_vec().into_boxed_slice();
    for p in sift(n as u64) {
        let p = p as i64;
        for (i, &v) in keys.iter().enumerate().rev() {
            if p > v {
                break;
            }
            s.arr[i] -= s[v / p];
        }
    }
    assert_eq!(s, mertens(n)); */
    //utils::primes::prime_sieves::main();
    //utils::primes::primecount::main();
    println!("Finished running at: {} ", Local::now().time());
}
