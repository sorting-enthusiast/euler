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
#[must_use]
pub fn is_prime(n: u64) -> bool {
    if n == 2 {
        return true;
    }
    if n & 1 == 0 {
        return false;
    }
    let bases = //[2, 325, 9375, 28178, 450775, 9780504, 1795265022]; //
    [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
    let s = (n - 1).trailing_zeros();
    let d = (n - 1) >> s;
    //dbg!(s, d);
    for a in bases {
        if a >= n {
            break;
        }
        let mut x = powmod(a, d, n);
        let mut y = 1;
        for _ in 0..s {
            y = mulmod(x, x, n);
            if y == 1 && x != 1 && x != n - 1 {
                return false;
            }
            x = y;
        }
        if y != 1 {
            return false;
        }
    }
    true
}
