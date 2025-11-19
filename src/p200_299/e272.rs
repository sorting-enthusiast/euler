use crate::utils::primes::prime_sieves::{sieve_it, sift};

const fn powmod(mut x: u64, mut exp: u64, modulo: u64) -> u64 {
    if exp == 0 {
        return 1;
    }
    let mut r = 1;
    x %= modulo;
    while exp > 1 {
        if exp & 1 == 1 {
            r = mulmod(r, x, modulo);
        }
        x = mulmod(x, x, modulo);
        exp >>= 1;
    }
    mulmod(r, x, modulo)
}
const fn mulmod(mut x: u64, mut exp: u64, modulo: u64) -> u64 {
    if exp == 0 {
        return 0;
    }
    let mut r = 0;
    x %= modulo;
    while exp > 1 {
        if exp & 1 == 1 {
            r += x;
            if r >= modulo {
                r -= modulo;
            }
        }
        x <<= 1;
        if x >= modulo {
            x -= modulo;
        }
        exp >>= 1;
    }
    r += x;
    if r >= modulo {
        r -= modulo;
    }
    r
}

// f_3(n) counts # of cube roots of unity mod n, including 1
// f_3(n) is multiplicative
// if p = 3k+1: f_3(p^e) = 3;
// if p = 3k+2: f_3(p^e) = 1;
// f_3(3) = 1, f_3(3^e), e>=2 = 3
// sum all numbers with f_3(n) = 243 = 3^5
const N: u64 = 1e11 as _;
pub fn main() {
    let roots = |n: u64| {
        dbg!(n);
        let mut count = 0;
        for i in 2..n {
            if powmod(i, 3, n) == 1 {
                count += 1;
            }
        }
        dbg!(count + 1);
    };
    let mut n = 29;
    while n <= 1e6 as u64 {
        roots(n);
        n *= 29;
    }

    dbg!(
        sift(dbg!(
            1e11 as usize
                / (sieve_it()
                    .filter(|p| p % 3 == 1)
                    .take(3)
                    .inspect(|p| println!("{p}"))
                    .product::<usize>()
                    * 9)
        ) as u64)
        .len()
    );
}
