use itertools::Itertools;

use crate::utils::{factorint::print_factors, primality::is_prime, prime_sieves::sieve_it};

// let f(n) count the # of gozinta chains for n:
// f(n) = sum over d|n, d<n of f(d)
// f(p^e) = 2^(e-1)
// up to 10^16 have at most 13 distinct prime factors
const N: usize = 1e9 as _;
pub fn main() {
    let mut goz = vec![1; N + 1];
    for i in 2..=N {
        if goz[i] == i {
            dbg!(i);
            dbg!(i.trailing_zeros());
            print_factors((i >> i.trailing_zeros()) as u64);
            println!();
        }
        for m in (2 * i..=N).step_by(i) {
            goz[m] += goz[i];
        }
    }
    dbg!(&goz[1..25]);
    let n: u64 = 296755361792;
    dbg!(n.trailing_zeros());
    print_factors(n >> n.trailing_zeros());
    println!();
}
// 4 1
// 8 1
// 6 1 1
// 12 1
// 4 4 1
// 12 1 1
// 14 1 1
// 20 1
// 18 1 1
// 24 1
