use crate::utils::primes::{primality::powmod, prime_sieves::sieve_it};
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
