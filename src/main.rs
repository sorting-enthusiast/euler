use std::time::Instant;

use crate::count_squarefree::count_squarefree;
pub mod achilles;
pub mod bit_array;
pub mod count_squarefree;
pub mod eratosthenes_variants;
pub mod longest_collatz_chain;
pub mod pandigital_products;
pub mod prime_sieves;
pub mod sieve_of_pritchard;
pub mod xorprimes;
const fn modular_exponentiation<const MODULO: u128>(mut x: u128, mut exp: u128) -> u128 {
    let mut r = 1;
    while exp > 1 {
        if exp & 1 == 1 {
            r = (r * x) % MODULO;
        }
        x = (x * x) % MODULO;
        exp >>= 1;
    }
    (r * x) % MODULO
}
fn gen_candidates_1(limit: u64, acc: u64, primes: &[u64]) -> Vec<u64> {
    let prime_len = primes.len();
    if prime_len == 0 {
        return vec![];
    }
    let mut res = vec![];
    for i in 1..=prime_len {
        let p = primes[prime_len - i];
        let mut prod = acc * (p * p);

        while prod <= limit {
            res.push(prod);
            res.extend(gen_candidates_1(limit, prod, &primes[..prime_len - i]));
            prod *= p;
        }
    }
    res
}
#[inline(never)]
pub fn gcd(mut u: u64, mut v: u64) -> u64 {
    if u == 0 || v == 0 {
        return u | v;
    }
    let shift = (u | v).trailing_zeros();
    u >>= u.trailing_zeros();
    loop {
        v >>= v.trailing_zeros();
        v -= u;
        let m = (v as i64 >> 63) as u64;
        u += v & m;
        v = (v + m) ^ m;
        if v == 0 {
            break;
        }
    }
    u << shift
}
pub fn main() {
    /* dbg!(modular_exponentiation::<10_000_000_000>(2, 7830457));
       dbg!(9700303872u128 * 28433 % 10_000_000_000);
       prime_sieves::main();
       longest_collatz_chain::main();
       const A: u64 = 51;
       const B: u64 = 257;
       dbg!(A);
       dbg!(B);

       dbg!(gcd(A, B));
       //193
       let start = Instant::now();
       dbg!(count_squarefree(1 << 50));
       let end = start.elapsed();
       println!("{:?}", end);
       //197
       let u_n = (2..513).fold(0.71, |acc, _| 1.42 * 2.0f64.powf(-acc * acc));
       dbg!(u_n + 1.42 * 2.0f64.powf(-u_n * u_n));

       //301: nim
       dbg!((1..=(1 << 30)).filter(|x| x ^ (x << 1) == 3 * x).count());
    */
    achilles::main();
    //xorprimes::main();
}
