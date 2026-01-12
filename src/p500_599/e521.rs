use itertools::Itertools;

use crate::utils::{FIArray::FIArray, multiplicative_function_summation::sum_n_usize};

const MOD: usize = 1e9 as usize;
const N: usize = 1e12 as _;
pub fn main() {
    let start = std::time::Instant::now();
    let mut sum = 0;
    let mut s = FIArray::new(N);
    let mut sums = s.clone();
    let keys = FIArray::keys(N).collect_vec();
    for (i, &v) in keys.iter().enumerate() {
        s.arr[i] = v - 1;
        sums.arr[i] = (sum_n_usize::<MOD>(v) - 1 + MOD) % MOD;
    }
    for p in 2..=s.isqrt {
        if s.arr[p - 1] == s.arr[p - 2] {
            continue;
        }
        let sp = s.arr[p - 2];
        let sp2 = sums.arr[p - 2];
        sum += (p * (s[N / p] - sp)) % MOD; // p = lpf
        if sum >= MOD {
            sum -= MOD;
        }
        for (i, &v) in keys.iter().enumerate().rev() {
            if v < p * p {
                break;
            }
            s.arr[i] -= s[v / p] - sp;
            sums.arr[i] += MOD - (p * (sums[v / p] + MOD - sp2)) % MOD;
            sums.arr[i] %= MOD;
        }
    }
    sum += sums[N] % MOD;
    if sum >= MOD {
        sum -= MOD;
    }
    println!("res = {sum}, took {:?}", start.elapsed());
}
