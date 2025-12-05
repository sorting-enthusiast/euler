use crate::utils::{
    FIArray::FIArray, math::gcd, multiplicative_function_summation::dirichlet_mul_usize,
};
const N: usize = 1e8 as _;
const SQRT_N: usize = N.isqrt();
pub fn main() {
    let start = std::time::Instant::now();
    let sum_divisors =
        dirichlet_mul_usize(&FIArray::unit(N), &FIArray::id::<{ usize::MAX >> 1 }>(N), N);
    let mut res = sum_divisors[N / 2];
    for a in 1..=SQRT_N {
        for b in a + 1..=SQRT_N {
            if gcd(a as _, b as _) == 1 {
                let n = a * a + b * b;
                if n > N {
                    break;
                }
                res += (a + b) * sum_divisors[N / n];
            }
        }
    }
    res <<= 1;
    res += sum_divisors[N];
    let end = start.elapsed();
    println!("res = {res}, took {end:?}");
}
