use itertools::Itertools;

use crate::utils::{
    FIArray::FIArrayI64,
    multiplicative_function_summation::{divisor_summatory, mobius_sieve},
    primes::wheel_sieve,
};

// can/should replace to something based on dirichlet hyperbola
fn count_squarefree(x: i64) -> FIArrayI64 {
    let mut s = FIArrayI64::unit(x);
    let keys = FIArrayI64::keys(x).collect_vec().into_boxed_slice();
    let primes = wheel_sieve(x.isqrt() as _).into_boxed_slice();

    for p in primes {
        let p = p as i64;
        for (i, &v) in keys.iter().enumerate().rev() {
            if v < p * p {
                break;
            }
            s.arr[i] -= s[v / (p * p)];
        }
    }
    s
}
const N: i64 = 1e12 as _;
const SQRT_N: i64 = N.isqrt();
// A018892
// sum of (d(n^2) + 1) / 2
// let d2(n) = d(n^2)
// d2 = d * sqf
// the rest is just dirichlet hyperbola
pub fn main() {
    let start = std::time::Instant::now();
    let sqf = count_squarefree(N);
    dbg!(start.elapsed());
    let d = divisor_summatory(N);
    dbg!(start.elapsed());
    let len = d.arr.len();
    let mut d2 = sqf[N] + d[N] - d[SQRT_N] * sqf[SQRT_N];
    for i in 2..=SQRT_N as usize {
        d2 += (d.arr[i - 1] - d.arr[i - 2]) * sqf.arr[len - i];
        d2 += (sqf.arr[i - 1] - sqf.arr[i - 2]) * d.arr[len - i];
    }
    d2 += N;
    d2 >>= 1;
    println!("res = {d2}, took {:?}", start.elapsed());
    alt();
}
// can optimise even more, to O(n^5/9) time and O(n^1/3) memory,
// using a combination of a better D3 summation function and Jakub Pawlewicz's squarefree counting method (i.e. just smarter dirichlet hyperbola)
fn alt() {
    let start = std::time::Instant::now();
    let d = divisor_summatory(N);
    let d3 = |x: i64| {
        let xsqrt = x.isqrt();
        let mut ret = x + d[x] - xsqrt * d[xsqrt];
        for i in 2..=xsqrt {
            ret += d[x / i];
            ret += (d.arr[i as usize - 1] - d.arr[i as usize - 2]) * (x / i);
        }
        ret
    };
    dbg!(start.elapsed());
    let mob = mobius_sieve(SQRT_N as usize + 1);
    let mut res = N + d3(N);
    for i in 2..=SQRT_N as usize {
        res += i64::from(mob[i]) * d3(N / (i * i) as i64);
    }
    res >>= 1;
    println!("res = {res}, took {:?}", start.elapsed());
}
