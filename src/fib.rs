use rug::{Complete, Integer};
pub fn main() {
    const N: usize = u32::MAX as _;
    println!("{N}");

    let start = std::time::Instant::now();
    let f2 = Integer::fibonacci(N as _).complete();
    let end = start.elapsed();
    println!("took {end:?}");

    let start = std::time::Instant::now();
    let f = fib_bostan_mori(N);
    let end = start.elapsed();
    println!("took {end:?}");

    assert_eq!(f, f2);
}
fn fib_bostan_mori(mut n: usize) -> Integer {
    if n < 2 {
        if n < 1 {
            return Integer::ZERO;
        }
        return Integer::ZERO + 1;
    }
    let mut c = Integer::from(3);
    let (mut a, mut b) = if n & 1 == 0 {
        (Integer::ZERO, Integer::ZERO + 1)
    } else {
        (Integer::ZERO + 1, Integer::ZERO - 1)
    };
    n >>= 1;
    while n > 1 {
        if n & 1 == 0 {
            b *= &c;
            b += &a;
        } else {
            a *= &c;
            a += &b;
        }
        c.square_mut();
        c -= 2;
        n >>= 1;
    }
    b + a * c
}
/* // copied from https://www.reddit.com/r/algorithms/comments/1im1kc5/20000000th_fibonacci_number_in_1_second/s
fn fib_luc(mut n: isize) -> (Integer, Integer) {
    if n == 0 {
        return (Integer::ZERO, Integer::from(2));
    }

    if n < 0 {
        n *= -1;
        let (fib, luc) = fib_luc(n);
        let k = n % 2 * 2 - 1;
        return (fib * k, luc * k);
    }

    if n & 1 == 1 {
        let (fib, luc) = fib_luc(n - 1);
        return (
            (&fib + &luc).complete() >> 1,
            ((5 * &fib).complete() + &luc).complete() >> 1,
        );
    }

    n >>= 1;
    let k = n % 2 * 2 - 1;
    let (fib, luc) = fib_luc(n);
    (fib * &luc, (&luc * &luc).complete() + 2 * k)
}
 */
