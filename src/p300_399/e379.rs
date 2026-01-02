use crate::utils::fast_divisor_sums::{d3, icbrt};

const N: i64 = 1e12 as _;
// A018892
// sum of (d(n^2) + 1) / 2
// let d2(n) = d(n^2)
// d2 = d * sqf = mu_sqrt * d3

// using d2 = mu_sqrt * d3
// O(n^2/3) time and O(n^1/3) space (can technically reduce to O(n^5/9) time, but won't necessarily be faster for required input)
// dirichlet hyperbola with n^2/3 - n^1/3 split
pub fn main() {
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

    const I: i64 = icbrt(N);
    const D: i64 = (N / I).isqrt();
    let start = std::time::Instant::now();
    let mut res = 0;
    let mut mertens_small = mobius_sieve(D as usize + 1);
    for d in 1..=D {
        if mertens_small[d as usize] != 0 {
            res += mertens_small[d as usize] * d3(N / (d * d));
        }
        mertens_small[d as usize] += mertens_small[d as usize - 1];
    }
    let mut small_diff = vec![0; I as usize - 1].into_boxed_slice();
    for i in 1..I {
        small_diff[i as usize - 1] = d3(i); // can use linear sieve, but the code is simpler this way and the time complexity doesn't change
    }
    res -= mertens_small[D as usize] * small_diff[I as usize - 2];

    for i in (2..I).rev() {
        small_diff[i as usize - 1] -= small_diff[i as usize - 2];
    }

    let mut mertens_big = vec![0; I as usize - 1].into_boxed_slice(); // indexed by denominator
    for i in (1..I).rev() {
        let v = (N / i).isqrt() as usize;
        let vsqrt = v.isqrt();
        let mut m = 1 - v as i64 + vsqrt as i64 * mertens_small[vsqrt];
        for d in 2..=vsqrt {
            m -= if v / d <= D as usize {
                mertens_small[v / d]
            } else {
                mertens_big[i as usize * d * d - 1]
            };
            m -= (mertens_small[d] - mertens_small[d - 1]) * (v / d) as i64;
        }
        mertens_big[i as usize - 1] = m;
        res += small_diff[i as usize - 1] * m;
    }
    res += N;
    res >>= 1;
    println!("res = {res}, took {:?}", start.elapsed());
}

// using d2 = d * sqf
// O(n^3/4) time, O(n^1/2) space
fn initial_solution() {
    use itertools::Itertools;

    use crate::utils::{
        FIArray::FIArrayI64, multiplicative_function_summation::divisor_summatory,
        primes::wheel_sieve,
    };
    const SQRT_N: i64 = N.isqrt();

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
}
