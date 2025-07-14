use crate::utils::prime_sieves::{MOD30_TO_MASK, WHEEL_2_3_5_7};
use std::time::Instant;

pub fn sieve(n: usize) -> Vec<usize> {
    if n < 7 {
        let mut ret = vec![2];
        if n > 2 {
            ret.push(3);
        }
        if n > 4 {
            ret.push(5);
        }
        return ret;
    }
    let bitmap_size = (((n / 30) >> 3) + 1) << 3;
    let mut sieve = vec![0u8; bitmap_size].into_boxed_slice();
    let bitmap = sieve.as_mut_ptr();
    let mut primes = Vec::with_capacity((n as f64 / (n as f64).log(3.0)) as usize);
    primes.extend([2, 3, 5, 7]);
    let mut num = 1;
    let mut wheel_incr = WHEEL_2_3_5_7.into_iter().cycle();

    loop {
        num += unsafe { wheel_incr.next().unwrap_unchecked() } as usize;
        if num * num > n {
            break;
        }
        if unsafe { *bitmap.add(num / 30) } & MOD30_TO_MASK[num % 30] == 0 {
            primes.push(num);
            let mut multiple = num * num;
            /*for incr in wheel_incr.clone() {
                if multiple > n {
                    break;
                }
                *bitmap.add(multiple / 30) |= MOD30_TO_MASK[num % 30];
                multiple += incr as usize * num;
            } */
            let precomp = [
                0,
                0,
                num << 1,
                0,
                num << 2,
                0,
                (num << 2) + (num << 1),
                0,
                num << 3,
                0,
                (num << 3) + (num << 1),
            ]; //to reduce bit complexity, only use adds and shifts instead of a multiplication every iteration of inner loop
            for incr in wheel_incr.clone() {
                if multiple > n {
                    break;
                }
                unsafe { *bitmap.add(multiple / 30) |= MOD30_TO_MASK[multiple % 30] };
                multiple += precomp[incr as usize];
            }
        }
    }
    for incr in wheel_incr {
        if num > n {
            break;
        }
        if unsafe { *bitmap.add(num / 30) & MOD30_TO_MASK[num % 30] } == 0 {
            primes.push(num);
        }
        num += incr as usize;
    }

    primes
}
// calculate product of (1+x^p) for all p < 5000, i.e. the generating function of all subsets of sum i, sum coefficients of prime powers
pub fn main() {
    const MOD: usize = 1e16 as usize;

    let start = Instant::now();

    let n = 1_548_136; // sum of primes below 5000
    let primes = sieve(n);
    let mut dp = vec![0; n + 1];
    dp[0] = 1;

    let mut max = 0;
    for &p in &primes[..669] {
        max += p; // largest possible sum at current iteration
        for j in (p..=max).rev() {
            dp[j] += dp[j - p];
            if dp[j] >= MOD {
                dp[j] -= MOD;
            }
        }
    }
    let mut sum = 0;
    for p in primes {
        sum += dp[p];
        if sum >= MOD {
            sum -= MOD;
        }
    }
    println!("{:?}", start.elapsed());
    println!("{sum}");
}
