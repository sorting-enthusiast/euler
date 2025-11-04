use crate::utils::{math::gcd, primality::is_prime};

pub fn rho(n: u64, x0: u64, c: u64) -> u64 {
    let f = |x| (x * x + c) % n;
    let mut x: u64 = x0;
    let mut y = x0;
    let mut g = 1;
    while g == 1 {
        x = f(x);
        y = f(f(y));
        g = gcd(x.abs_diff(y), n);
    }
    g
}
pub fn print_factors(n: u64) {
    if n == 1 {
        return;
    }
    for (x0, c) in [(2, 1), (2, 2), (3, 1), (3, 2)] {
        let factor = rho(n, x0, c);
        if factor != n {
            print_factors(factor);
            print_factors(n / factor);
            return;
        }
    }
    if is_prime(n) {
        print!("{n} ");
    } else {
        print!("({n}) ");
    }
}
