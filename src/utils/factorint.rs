use std::num::NonZeroU64;

use crate::utils::{
    math::gcd,
    primality::{is_prime, mulmod},
};

pub fn rho(n: u64) -> u64 {
    let mut x = 0;
    let mut y = 0;
    let mut t = 30;
    let mut prd = 2;
    let mut i = 1;
    let f = |x, i| mulmod(x, x, n) + i;
    while t % 40 != 0 || gcd(prd, n) == 1 {
        if x == y {
            i += 1;
            x = i;
            y = f(x, i);
        }
        prd = NonZeroU64::new(mulmod(prd, x.abs_diff(y), n)).map_or(prd, u64::from);
        x = f(x, i);
        y = f(f(y, i), i);

        t += 1;
    }
    gcd(prd, n)
}
pub fn print_factors(n: u64) {
    if n == 1 {
        return;
    }
    if is_prime(n) {
        print!("{n} ");
        return;
    }
    let d = rho(n);
    print_factors(d);
    print_factors(n / d);
}

pub fn factor(n: u64, factors: &mut Vec<u64>) {
    if n == 1 {
        return;
    }
    if is_prime(n) {
        factors.push(n);
        return;
    }
    let d = rho(n);
    factor(d, factors);
    factor(n / d, factors);
}
