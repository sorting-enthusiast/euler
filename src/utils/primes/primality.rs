#[must_use]
const fn powmod_(mut x: u128, mut exp: u128, modulo: u128) -> u128 {
    if exp == 0 {
        return 1;
    }
    let mut r = 1;
    x %= modulo;
    while exp > 1 {
        if exp & 1 == 1 {
            r = mulmod_(r, x, modulo);
        }
        x = mulmod_(x, x, modulo);
        exp >>= 1;
    }
    mulmod_(r, x, modulo)
}
#[must_use]
pub const fn mulmod_(mut x: u128, mut exp: u128, modulo: u128) -> u128 {
    /*if exp == 0 {
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
    (r + x) % modulo*/
    (x * exp) % modulo
}
#[must_use]
pub fn is_prime(n: u64) -> bool {
    assert!(n <= 7e18 as u64);
    let n = n as u128;
    if n < 2 || (n % 6) % 4 != 1 {
        return (n | 1) == 3;
    }
    let bases = [2u128, 325, 9375, 28178, 450_775, 97_805_044, 1_795_265_022];
    let s = (n - 1).trailing_zeros();
    let d = (n - 1) >> s;

    for a in bases {
        if a.is_multiple_of(n) {
            continue;
        }
        let mut p = powmod_(a, d, n);
        if p == 1 || p == n - 1 {
            continue;
        }
        for _ in 0..s {
            p = mulmod_(p, p, n);
            if p == 1 || p == n - 1 {
                break;
            }
        }

        if p != n - 1 {
            return false;
        }
    }
    true
}
#[must_use]
pub const fn powmod(mut x: u64, mut exp: u64, modulo: u64) -> u64 {
    powmod_(x as _, exp as _, modulo as _) as _
}
#[must_use]
pub const fn mulmod(mut x: u64, mut exp: u64, modulo: u64) -> u64 {
    mulmod_(x as _, exp as _, modulo as _) as _
}
