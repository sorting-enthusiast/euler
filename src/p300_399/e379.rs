use crate::{
    p300_399::e362::mult,
    utils::{
        FIArray::FIArray,
        fast_divisor_sums::{d3, icbrt},
        multiplicative_function_summation::{
            count_squarefree, dirichlet_mul_single_i64, dirichlet_mul_single_usize, divisor_summatory, sqf
        },
    },
};
const N: i64 = 1e12 as _;
// A018892
// sum of (d(n^2) + 1) / 2
// let d2(n) = d(n^2)
// d2 = d * sqf = mu_sqrt * d3

// using d(n^2) = (mu_sqrt * d3)(n)
// O(n^5/9) time and O(n^1/3) space
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
    dbg!(I,D);
    let start = std::time::Instant::now();
    let mut res = 0;
    let mut mertens_small = mobius_sieve(D as usize + 1);
    dbg!(start.elapsed());
    for d in 1..=D {
        if mertens_small[d as usize] != 0 {
            let v= d3((N / (d * d)) as _) as i64;
            res += mertens_small[d as usize] as i64 * v;
        }
        mertens_small[d as usize] += mertens_small[d as usize - 1];
    }
    dbg!(start.elapsed());
    let mut small_diff = vec![0; I as usize - 1].into_boxed_slice();
    for i in 1..I {
        small_diff[i as usize - 1] = d3(i as _) as i64; // can use linear sieve, but the code is simpler this way and the time complexity doesn't change
    }
    dbg!(start.elapsed());
    res -= mertens_small[D as usize] as i64 * small_diff[I as usize - 2];

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
        res += small_diff[i as usize - 1] * m as i64;
    }
    res += N as i64;
    res >>= 1;
    println!("res = {res}, took {:?}", start.elapsed());
    initial_solution();
}

// using d2 = d * sqf
// O(n^2/3 \log^-c n) time, O(n^1/2) space
fn initial_solution() {
    const N: usize = self::N as _;
    const SQRT_N: usize = N.isqrt();

    let start = std::time::Instant::now();
    let sqf = sqf(N);//count_squarefree(N);
    /* let u = FIArray::unit(N);
    let d = mult(&u, &u); */
 //dirichlet_mul_usize(&u, &u, N);
    let d = divisor_summatory(N);
    let res = (N + dirichlet_mul_single_i64(&sqf, &d) as usize) >> 1;
    println!("res = {res}, took {:?}", start.elapsed());
}
