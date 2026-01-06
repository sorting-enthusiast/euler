use crate::utils::{multiplicative_function_summation::mobius_sieve, primes::primecount::prime_pi};

const N: usize = 9e18 as _;
const SQRT_N: usize = N.isqrt();
const CBRT_N: usize = icbrt(N);
const SIXTH_ROOT: usize = CBRT_N.isqrt();

const fn icbrt(x: usize) -> usize {
    let mut rt = 1 << (1 + x.ilog2().div_ceil(3));
    let mut x_div_rt2 = (x / rt) / rt;
    while rt > x_div_rt2 {
        rt = ((rt << 1) + x_div_rt2) / 3;
        x_div_rt2 = (x / rt) / rt;
    }
    rt
}
// all powerful numbers - squarefree numbers below cbrt - cubefree numbers below sqrt - pi(sixth root) + 1
pub fn main() {
    let start = std::time::Instant::now();
    //dbg!(sqrt, cbrt, sixth_root);
    let mut res = 1;
    let mobius = mobius_sieve(CBRT_N + 1);
    for b in 1..=CBRT_N {
        if mobius[b] == 0 {
            continue;
        }
        let bbb = b.pow(3);
        res += (N / bbb).isqrt() - 1;
        if bbb <= SQRT_N {
            if mobius[b] == 1 {
                res -= SQRT_N / bbb;
            } else {
                res += SQRT_N / bbb;
            }
        }
    }
    res -= prime_pi(SIXTH_ROOT);
    let end = start.elapsed();
    println!("res = {res}, took {end:?}");
}
