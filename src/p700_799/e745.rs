use std::time::Instant;

use crate::utils::powerful_numbers::PowerfulExtSkipZero;
const N: i64 = 1e14 as i64;
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

// powerful number trick ftw: g(p) = 1(p) = 1
// g(p^e) = p^(2 * floor(e/2))
// h(p^e) = (f*g^-1)(p^e) = 0 if e odd, p^(e-2) * (p^2 - 1) if e even
pub fn main() {
    let start = Instant::now();

    let h = |p, e| {
        if e & 1 == 0 {
            (powmod(p, e - 2) * ((p * p - 1) % MOD)) % MOD
        } else {
            0
        }
    };
    let mut sum = 0;
    for (n, hn) in PowerfulExtSkipZero::<_, MOD>::new(N, h) {
        sum += hn * ((N / n) % MOD);
        sum %= MOD;
    }

    let end = start.elapsed();
    println!("res = {sum}, took {end:?}");
}
