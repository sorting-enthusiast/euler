#![warn(clippy::pedantic)]
#![allow(long_running_const_eval)]
#![allow(dead_code)]
#![allow(clippy::too_many_lines)]
#![allow(clippy::cast_possible_truncation)]
#![allow(clippy::inline_always)]
#![allow(clippy::cast_sign_loss)]
#![allow(clippy::cast_precision_loss)]

use crate::{e193::count_squarefree, utils::powerful_numbers::PowerfulExt};

mod e107;
mod e12;
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
mod e27;
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
mod e942;
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

// digital root of n is just n mod 9 if n mod 9 != 0, otherwise 9
pub fn main() {
    //e193::main();
    e942::main();
    //utils::primecount::main();
    /*let n = 2e4 as i64;
    let prime_pi = utils::primecount::lucy(n as usize);
    let sqfree_1 = count_squarefree(n.isqrt() as _) as usize;
    let sqfree_2 = count_squarefree((n as f64).cbrt() as _) as usize;
    dbg!(dbg!(PowerfulExt::<_, 2>::new(n, |_, _| 1).count()) - dbg!(sqfree_1) - dbg!(sqfree_2) + 1); */
}
