#![warn(clippy::pedantic)]
#![allow(long_running_const_eval)]
#![allow(dead_code)]
#![allow(clippy::too_many_lines)]
#![allow(clippy::cast_possible_truncation)]
#![allow(clippy::inline_always)]
#![allow(clippy::cast_sign_loss)]
#![allow(clippy::cast_precision_loss)]

use chrono::Local;

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
    //p200_299::e268::main();
    //p500_599::e580::main();
    //p100_199::e193::main();
    //p600_699::e634::main();
    //dbg!(opt(1e16 as _));
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
    /* let set = FIArray::keys(1e16 as _)
        .map(usize::isqrt)
        .collect::<HashSet<_>>();
    dbg!(set.len()); */
    //utils::prime_sieves::main();
    utils::primecount::main();
    println!("Finished running at: {} ", Local::now().time());
}
