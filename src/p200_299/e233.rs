use std::collections::BTreeMap;

use itertools::{Either, Itertools};

use crate::utils::primes::prime_sieves::sift;
const N: u64 = 1e11 as _;

// sum all numbers below lim s.t. all their prime factors are from the prime set given
fn dfs(acc: u64, lim: u64, primes: &[u64]) -> u64 {
    let mut sum = acc;
    for (i, &p) in primes.iter().enumerate() {
        if p * acc > lim {
            break;
        }
        let mut power = p;
        loop {
            sum += dfs(acc * power, lim, &primes[i + 1..]);
            power *= p;
            if power * acc > lim {
                break;
            }
        }
    }
    sum
}

pub fn main() {
    let start = std::time::Instant::now();

    let limit = N / (5u64.pow(3) * 13u64.pow(2));
    let primes = sift(limit + 1);
    let (primes_1mod4, rest): (Vec<_>, Vec<_>) = primes.into_iter().partition(|p| p & 3 == 1);
    let mut sum = 0;
    let mut cache = BTreeMap::new(); // few possible values for N / init, may as well cache them. 5ms speedup, 25%
    //let mut count = 0;
    // one prime contributes 7, one 5, one 3 for product of 105
    for p1 in &primes_1mod4 {
        let p1_cubed = p1.pow(3);
        if p1_cubed > N {
            break;
        }
        for p2 in &primes_1mod4 {
            if p2 == p1 {
                continue;
            }
            let p2_squared = p2.pow(2);
            if p1_cubed * p2_squared > N {
                break;
            }
            for p3 in &primes_1mod4 {
                if p3 == p1 || p3 == p2 {
                    continue;
                }
                let init = p1_cubed * p2_squared * p3;
                if init > N {
                    break;
                }
                sum += init
                    * *cache
                        .entry(N / init)
                        .or_insert_with(|| dfs(1, N / init, &rest));

                //count += 1;
            }
        }
    }

    // Reminder: dumb fucking mistake, forgot to consider other ways of partitioning 105 into divisors at first :(.
    // one prime contributes 15, one 7
    for p1 in &primes_1mod4 {
        let p1_seventh = p1.pow(7);
        if p1_seventh > N {
            break;
        }
        for p2 in &primes_1mod4 {
            if p2 == p1 {
                continue;
            }

            let init = p1_seventh * p2.pow(3);
            if init > N {
                break;
            }
            sum += init
                * *cache
                    .entry(N / init)
                    .or_insert_with(|| dfs(1, N / init, &rest));
            //count += 1;
        }
    }
    // one prime contributes 21, one 5
    for p1 in &primes_1mod4 {
        let p1_tenth = p1.pow(10);
        if p1_tenth > N {
            break;
        }
        for p2 in &primes_1mod4 {
            if p2 == p1 {
                continue;
            }

            let init = p1_tenth * p2.pow(2);
            if init > N {
                break;
            }
            sum += init
                * *cache
                    .entry(N / init)
                    .or_insert_with(|| dfs(1, N / init, &rest));
            //count += 1;
        }
    }
    // other combinations unnecessary, wont fit in range
    let end = start.elapsed();
    println!("res = {sum}, took {end:?}");
    //dbg!(count);
    dbg!(cache.len());
}
