use itertools::Itertools;

use crate::utils::{
    FIArray::{DirichletFenwick, FIArray, FIArrayI64},
    multiplicative_function_summation::sum_n_i64,
    primes::prime_sieves::sift,
};
// TODO: optimize to O(n^2/3)
const N: usize = 2018_2018_2018;
const SQRT_N: usize = N.isqrt();
const MOD: usize = 1e9 as _;
pub fn main() {
    let start = std::time::Instant::now();
    let primes = sift(SQRT_N as _);
    println!("Sieved primes: {:?}", start.elapsed());
    let mut pi = FIArray::new(N);
    let mut pi_sums = FIArrayI64::new(N);
    let keys = FIArray::keys(N).collect_vec();

    for (i, &v) in keys.iter().enumerate() {
        pi.arr[i] = v - 1;
        pi_sums.arr[i] = (sum_n_i64::<{ MOD as _ }>(v as _) + MOD as i64 - 1) % MOD as i64;
    }

    for &p in &primes {
        let p = p as usize;
        let sp = pi.arr[p - 2];
        let sp2 = pi_sums.arr[p - 2];

        for (i, &v) in keys.iter().enumerate().rev() {
            if v < p * p {
                break;
            }
            pi.arr[i] -= pi[v / p] - sp;
            pi_sums.arr[i] -= p as i64 * (pi_sums[v / p] - sp2);
            pi_sums.arr[i] %= MOD as i64;
            if pi_sums.arr[i] < 0 {
                pi_sums.arr[i] += MOD as i64;
            }
        }
    }
    let mut sum = ((1..SQRT_N).fold(0, |acc, k| acc + pi_sums[N / k] as usize)
        - (SQRT_N - 1) * pi_sums[const { N / SQRT_N }] as usize)
        % MOD;
    if sum >= MOD {
        sum -= MOD;
    }
    println!("Initialized sum: {:?}", start.elapsed());

    let mut smooth = DirichletFenwick::eps(N);
    let split = primes.partition_point(|p| p * p * p <= N as u64);
    for &p in &primes[..split] {
        let p = p as usize;
        //smooth.sparse_mul_unlimited(p, 1);
        {
            let x = p;
            let w = 1;
            let lim = smooth.x / x;
            let mut prev = 0;
            let mut i = 1;
            while i <= lim / i {
                let cur = smooth.bit.sum(i - 1);
                if cur != prev {
                    smooth.bit.add(smooth.get_index(i * x), w * (cur - prev));
                    prev = cur;
                }
                i += 1;
            }
            for j in (p..=lim / i).rev() {
                let cur = smooth.get_prefix(lim / j);
                if cur != prev {
                    smooth.bit.add(smooth.bit.0.len() - j, w * (cur - prev));
                    prev = cur;
                }
            }
        }
        sum += p * (smooth.get_prefix(N / p) % MOD);
        sum %= MOD;
    }
    for &p in &primes[split..] {
        // optimization, not required for 1-minute rule
        let p = p as usize;
        sum +=
            p * (N / p - (1..=(N / p) / p).fold(0, |acc, k| acc + pi[(N / p) / k] - pi.arr[p - 1]));
        sum %= MOD;
    }
    println!("res = {sum}, took {:?}", start.elapsed());
}
