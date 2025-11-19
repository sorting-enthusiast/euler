use crate::utils::{multiplicative_function_summation::mobius_sieve, primes::primecount::prime_pi};

const N: usize = 9e18 as _;
fn icbrt(x: usize) -> usize {
    let mut rt = 0;
    let mut rt_squared = 0;
    let mut rt_cubed = 0;
    while rt_cubed <= x {
        rt += 1;
        rt_squared += 2 * rt - 1;
        rt_cubed += 3 * rt_squared - 3 * rt + 1;
    }
    rt - 1
}
// all powerful numbers - squarefree numbers below cbrt - cubefree numbers below sqrt - pi(sixth root) + 1
pub fn main() {
    let start = std::time::Instant::now();
    let sqrt = N.isqrt();
    let cbrt = icbrt(N);
    let sixth_root = cbrt.isqrt();
    //dbg!(sqrt, cbrt, sixth_root);
    let mut res = 1;
    let mobius = mobius_sieve(cbrt + 1);
    for b in 1..=cbrt {
        if mobius[b] == 0 {
            continue;
        }
        let bbb = b.pow(3);
        res += (N / bbb).isqrt() - 1;
        if bbb <= sqrt {
            if mobius[b] == 1 {
                res -= sqrt / bbb;
            } else {
                res += sqrt / bbb;
            }
        }
    }
    res -= prime_pi(sixth_root);
    let end = start.elapsed();
    println!("res = {res}, took {end:?}");
}
