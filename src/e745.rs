use std::time::Instant;

use crate::utils::{math::modexp, powerful_numbers::PowerfulExt};
const N: i64 = 1e14 as i64;
const MOD: i64 = 1e9 as i64 + 7;
// powerful number trick ftw: g(p) = 1(p) = 1
// g(p^e) = p^(2 * floor(e/2))
// h(p^e) = (f*g^-1)(p^e) = 0 if e odd, p^(e-2) * (p^2 - 1) if e even
pub fn main() {
    let start = Instant::now();

    let h = |p, e| {
        if e & 1 == 0 {
            (modexp::<{ MOD as u128 }>(p as _, e as u128 - 2) as i64 * ((p * p - 1) % MOD)) % MOD
        } else {
            0
        }
    };
    let mut sum = 0;
    for (n, hn) in PowerfulExt::<_, MOD>::new(N, h).filter(|&(_, hn)| hn != 0) {
        sum += hn * ((N / n) % MOD);
        sum %= MOD;
    }

    let end = start.elapsed();
    println!("res = {sum}, took {end:?}");
}
