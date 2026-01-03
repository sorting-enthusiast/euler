use std::collections::HashMap;

use crate::utils::fast_divisor_sums::icbrt;

const N: i64 = 3_141_592_653_589_793i64;
const SQRT_N: i64 = N.isqrt();
fn mobius_sieve(n: usize) -> Vec<i64> {
    let mut res = vec![0; n];
    if n < 2 {
        return res;
    }
    let mut composite = crate::utils::bit_array::BitArray::zeroed(n);
    let mut primes = vec![];
    res[1] = 1;
    for i in 2..n {
        if !composite.get(i) {
            primes.push(i);
            res[i] = -1;
        }
        for &p in &primes {
            if i * p >= n {
                break;
            }
            composite.set(i * p);
            if i % p == 0 {
                res[i * p] = 0;
                break;
            }
            res[i * p] = -res[i];
        }
    }
    res
}
// can optimize to O(n^1/3) space, but I'm lazy:
// all calls to coprime points is with input N/2^2k or (N/2)/2^2k
// can precompute M for all values of N/d^2 and (N/2)/d^2 in O(n^1/3) space
// same for all values of all_points
pub fn main() {
    let start = std::time::Instant::now();
    let mut mu = mobius_sieve(SQRT_N as usize + 1);
    for i in 1..=SQRT_N as usize {
        mu[i] += mu[i - 1];
    }
    let mut cache = HashMap::new();
    // |{(x,y): x > 0, y >= 0, x^2 + y^2 <= t}| a.k.a. sum of u * \chi_4
    let mut all_points = |t: i64| {
        *cache.entry(t).or_insert_with(|| {
            let mut res = 0;
            let tsqrt = t.isqrt();
            for k in 1..=tsqrt {
                let tk = t / k;
                res += [0, 1, 1, 0][tk as usize & 3];
                res += [0, 1, 0, -1][k as usize & 3] * tk;
            }
            res - [0, 1, 1, 0][tsqrt as usize & 3] * tsqrt
        })
    };
    // |{(x,y): x > 0, y >= 0, gcd(x, y) = 1, x^2 + y^2 <= t}|
    let mut coprime_points = |t: i64| {
        let I = icbrt(t);
        let D = (t / I).isqrt();
        let mut res = all_points(t);
        for d in 2..=D {
            let mu_d = mu[d as usize] - mu[d as usize - 1];
            if mu_d != 0 {
                res += mu_d * all_points(t / (d * d));
            }
        }
        for i in 1..I {
            let v = (t / i).isqrt() as usize;
            res += mu[v] * (all_points(i) - all_points(i - 1));
        }
        res -= mu[D as usize] * all_points(I - 1);
        res
    };
    let mut res = 0;
    let mut lim = N;
    let mut sign = 1;
    while lim != 1 {
        res += sign * coprime_points(lim);
        sign = -sign;
        lim >>= 1;
    }
    // res = |{(x,y): x > 0, y >= 0, gcd(x, y) = 1, (x ^ y) & 1 = 0, x^2 + y^2 <= t}|
    res -= 1; // remove (1, 0)
    res >>= 1; // remove order between remaining (x,y) pairs
    println!("res = {res}, took {:?}", start.elapsed());
}
