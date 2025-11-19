use crate::utils::primes::prime_sieves::sieve_it;

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
            r = (r + x) % modulo;
        }
        x = (x + x) % modulo;
        exp >>= 1;
    }
    (r + x) % modulo
}

// f_3(n) counts # of cube roots of unity mod n, including 1
// f_3(n) is multiplicative
// if p = 3k+1: f_3(p^e) = 3;
// if p = 3k+2: f_3(p^e) = 1;
// f_3(3) = 1, f_3(3^e), e>=2 = 3
// sum all numbers with f_3(n) = 243 = 3^5
pub fn main() {
    let roots = |n: u64| {
        dbg!(n);
        let mut count = 0;
        let mut sum = 1;
        for i in 2..n {
            if powmod(i, 3, n) == 1 {
                count += 1;
                sum += i;
                dbg!(i);
            }
        }
        dbg!(count + 1);
        dbg!((dbg!(sum / n) * n, sum));
    };
    roots(7 * 13 * 2);
    roots(7 * 13);
    roots(7 * 2 * 2 * 2 * 2);
    roots(7 * 2 * 2 * 2);
    roots(7 * 2 * 2);
    roots(7 * 2);
    roots(7);

    /* roots(
        sieve_it()
            .take(14)
            .filter(|p| p % 3 == 1)
            .product::<usize>() as u64,
    ); */
    // n * 365
    //dbg!(2u128 * 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23 * 29 * 31 * 37 * 41 * 43);
}
