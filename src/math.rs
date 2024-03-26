pub const fn modular_exponentiation<const MODULO: u128>(mut x: u128, mut exp: u128) -> u128 {
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
pub const fn gcd(mut u: u64, mut v: u64) -> u64 {
    if u == 0 || v == 0 {
        return u | v;
    }
    let shift = (u | v).trailing_zeros();
    u >>= u.trailing_zeros();
    loop {
        v >>= v.trailing_zeros();
        if u > v {
            let tmp = u;
            u = v;
            v = tmp;
        }
        v -= u;
        /* let m = (v as i64 >> 63) as u64;
        u += v & m;
        v = (v + m) ^ m; */
        if v == 0 {
            break;
        }
    }
    u << shift
}
pub const fn ipow(mut x: u64, mut exp: u64) -> u64 {
    let mut r = 1;
    while exp > 1 {
        if exp & 1 == 1 {
            r *= x;
        }
        x *= x;
        exp >>= 1;
    }
    r * x
}
