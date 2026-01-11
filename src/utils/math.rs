#[must_use]
pub const fn modmul<const MODULO: u128>(mut x: u128, mut y: u128) -> u128 {
    if x == 0 || y == 0 {
        return 0;
    }
    let mut r = 0;
    while y > 1 {
        if y & 1 == 1 {
            r = (r + x) % MODULO;
        }
        x = (x + x) % MODULO;
        y >>= 1;
    }
    (r + x) % MODULO
}
#[must_use]
pub const fn modexp<const MODULO: u128>(mut x: u128, mut exp: u128) -> u128 {
    if exp == 0 {
        return 1;
    }
    let mut r = 1;
    while exp > 1 {
        if exp & 1 == 1 {
            r = modmul::<MODULO>(r, x);
        }
        x = modmul::<MODULO>(x, x);
        exp >>= 1;
    }
    modmul::<MODULO>(r, x)
}
#[must_use]
pub const fn iroot<const k: usize>(x: usize) -> usize {
    let mut rt = 1usize << (1 + x.ilog2().div_ceil(k as _));
    let mut x_div_rtk1 = x / rt.pow(k as u32 - 1);
    while rt > x_div_rtk1 {
        rt = (rt * (k - 1) + x_div_rtk1) / k;
        x_div_rtk1 = x / rt.pow(k as u32 - 1);
    }
    rt
}
#[must_use]
pub const fn prime_modinv<const MODULO: u128>(x: u128) -> u128 {
    modexp::<MODULO>(x, MODULO - 2)
}
#[must_use]
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
#[must_use]
pub const fn mult(x: u64, y: u64) -> u64 {
    let mut res = 0;
    let (mut base, mut mul) = if x.count_ones() > y.count_ones() {
        (x, y)
    } else {
        (y, x)
    };
    if mul & 1 == 0 {
        let shift = mul.trailing_zeros();
        mul >>= shift;
        base <<= shift;
    }
    while mul > 0 {
        res += base;
        mul ^= 1;
        let shift = mul.trailing_zeros();
        mul >>= shift;
        base <<= shift;
    }
    res
}
