use itertools::Itertools;

use crate::utils::{
    FIArray::{FIArray, FIArrayI64},
    multiplicative_function_summation::sum_n_i64,
    prime_sieves::sift,
};

const N: usize = 2018_2018_2018;
const SQRT_N: usize = N.isqrt();
const MOD: usize = 1e9 as _;
pub fn main() {
    let start = std::time::Instant::now();
    let primes = sift(SQRT_N as _);
    println!("Sieved primes: {:?}", start.elapsed());
    let mut pi = FIArray::new(N as _);
    let mut pi_sums = FIArrayI64::new(N as _);
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
            pi_sums.arr[i] -= p as i64 * (pi_sums[(v / p) as _] - sp2);
            pi_sums.arr[i] %= MOD as i64;
            if pi_sums.arr[i] < 0 {
                pi_sums.arr[i] += MOD as i64;
            }
        }
    }
    let mut sum = ((1..SQRT_N).fold(0, |acc, k| acc + pi_sums[(N / k) as _] as usize)
        - (SQRT_N - 1) * pi_sums[(const { N / SQRT_N }) as _] as usize)
        % MOD;
    if sum >= MOD {
        sum -= MOD;
    }
    println!("Initialized sum: {:?}", start.elapsed());

    let mut smooth = FIArray::new(N);
    for (i, &v) in keys.iter().enumerate() {
        smooth.arr[i] = 1 + v.ilog2() as usize;
    }

    sum += 2 * smooth[N / 2];
    sum %= MOD;
    let split = primes.partition_point(|p| p * p * p <= N as u64);
    for &p in &primes[1..] {
        if p < primes[split] {
            let p = p as usize;
            for (i, &v) in keys.iter().enumerate().skip(p - 1) {
                if v > N / p {
                    break; // will never read these values, can skip computing them
                }
                smooth.arr[i] += smooth[v / p];
                if smooth.arr[i] >= MOD {
                    smooth.arr[i] -= MOD;
                }
            }
            //dbg!(p);
            sum += p * smooth[N / p];
            sum %= MOD;
        } else {
            // optimization, not required for 1-minute rule
            let p = p as usize;
            sum += p
                * (N / p
                    - (1..=(N / p) / p).fold(0, |acc, k| acc + pi[(N / p) / k] - pi.arr[p - 1]));
            sum %= MOD;
        }
    }

    println!("res = {sum}, took {:?}", start.elapsed());
}
