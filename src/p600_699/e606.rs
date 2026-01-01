use crate::utils::{
    multiplicative_function_summation::{sum_n_i64, sum_over_primes},
    primes::prime_sieves::sift,
};
const MOD: i64 = 1e9 as _;
const N: i64 = 1e12 as i64;
// A034776

// finding the actual solution was quite easy,
// although I wasted an unbelievable amount of time debugging the code,
// turns out i was subtracting 1 from the sum in F wrong,
// not making sure the result mod MOD was nonnegative

// let f(n) count the # of gozinta chains for n:
// f(n) = sum over d|n, d<n of f(d)
// f(p^e) = 2^(e-1)
// only numbers for which f(n) = 252 are numbers of the form p^3q^3 = (pq)^3
// for 2 distinct primes p,q
// to sum all possible n, we iterate over p <= 10^6, adding p^3 * the sum of cubes of
// all primes between p + 1 and floor(10^12 / p)
// sums are calculated using lucy's algorithm.
pub fn main() {
    let start = std::time::Instant::now();

    let cubed = |p: i64| {
        let p = p % MOD;
        ((p * p) % MOD * p) % MOD
    };
    let sum_cubes = |v: i64| {
        // sum of cubes - 1
        let sn = sum_n_i64::<MOD>(v);
        ((sn * sn) % MOD + MOD - 1) % MOD
    };
    let s = sum_over_primes::<MOD>(N, cubed, sum_cubes);
    let mut sum = 0;
    for p in sift(N.isqrt() as u64) {
        let p = p as i64;
        sum += (cubed(p) * (s[N / p] + MOD - s.arr[p as usize - 1]) % MOD) % MOD;
        if sum >= MOD {
            sum -= MOD;
        }
    }
    let end = start.elapsed();
    println!("{sum}, {end:?}");
}
