use std::{
    arch::x86_64::{
        _mm_abs_epi32, _mm_blendv_ps, _mm_cvtsi128_si32, _mm_min_epi32, _mm_set1_epi32,
        _mm_slli_epi32, _mm_srli_epi32, _mm_sub_epi32,
    },
    mem::transmute,
};
pub mod bit_array;
pub mod longest_collatz_chain;
pub mod pandigital_products;
pub mod prime_sieves;
pub mod xorprimes;
const fn modular_exponentiation<const MODULO: u128>(mut x: u128, mut exp: u128) -> u128 {
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
#[inline(never)]
fn gcd(mut a: i32, mut b: i32) -> i32 {
    if a == 0 || b == 0 {
        return a | b;
    }
    let az = a.trailing_zeros();
    let bz = b.trailing_zeros();
    let shift = if az < bz { az } else { bz };
    a >>= az;
    b >>= shift;
    unsafe {
        let mut ar = _mm_set1_epi32(a);
        let mut br = _mm_set1_epi32(b);

        loop {
            let b_even = _mm_srli_epi32(br, 1);

            let diff = _mm_sub_epi32(ar, br);
            let a_odd = _mm_min_epi32(ar, br);
            let b_odd = _mm_abs_epi32(diff);

            let t = _mm_slli_epi32(br, 31);

            ar = transmute(_mm_blendv_ps(transmute(ar), transmute(a_odd), transmute(t)));

            br = transmute(_mm_blendv_ps(
                transmute(b_even),
                transmute(b_odd),
                transmute(t),
            ));

            if _mm_cvtsi128_si32(br) == 0 {
                break;
            }
        }
        _mm_cvtsi128_si32(ar) << shift
    }
}
pub fn main() {
    dbg!(modular_exponentiation::<10_000_000_000>(2, 7830457));
    dbg!(9700303872u128 * 28433 % 10_000_000_000);
    prime_sieves::main();
    longest_collatz_chain::main();
    const A: i32 = 51;
    const B: i32 = 257;
    dbg!(A);
    dbg!(B);
    dbg!(gcd(A, B));

    //xorprimes::main();
}
