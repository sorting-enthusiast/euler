pub mod achilles;
pub mod bit_array;
pub mod count_squarefree;
pub mod eratosthenes_variants;
pub mod longest_collatz_chain;
pub mod math;
pub mod pandigital_products;
pub mod prime_sieves;
pub mod sieve_of_pritchard;
pub mod xorprimes;

use std::time::Instant;

use crate::count_squarefree::count_squarefree;
pub fn main() {
    dbg!(math::modular_exponentiation::<10_000_000_000>(2, 7830457));
    dbg!(9700303872u128 * 28433 % 10_000_000_000);
    //prime_sieves::main();
    //longest_collatz_chain::main();
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

    //123: prime square remainder
    let start = Instant::now();
    let primes = sieve_of_pritchard::sift(1e6 as u64);
    for n in (7037..).step_by(2) {
        if (n << 1) * primes[n - 1] as usize > 1e10 as usize {
            dbg!(n);
            break;
        }
    }
    let end = start.elapsed();
    println!("{:?}", end);

    //120: square remainder
    let mut sum = 0;
    for a in 3..1001 {
        sum += a * ((a - 1) & !1);
    }
    dbg!(sum);

    //achilles::main();
    //xorprimes::main();
}
