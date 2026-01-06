use itertools::Itertools;

use crate::utils::{
    FIArray::FIArrayU64,
    primes::prime_sieves::{sieve_it, sift},
};

// Hilbert numbers are numbers that are products of any number of primes equivalent to 1 mod 4,
// and an even number of primes equivalent to 3 mod 4

// a hilbert number is h-squarefree if it is either squarefree, or if only 1 of its prime factors = 3 mod 4 appears with multiplicity 1< <4
const N: u64 = 1e18 as _;
const SQRT_N: u64 = N.isqrt();
pub fn main() {
    println!("res = {}", 2_327_213_148_095_366_u64);
    count_sqf(); // lol, can be made even faster using methods for counting primes in arithmetic progressions
    //solve_alt(); // case in point
    solve(); // can definitely optimize this more
    trivial();
}

fn trivial() {
    fn incl_excl(limit: u64, primes: &[u64]) -> u64 {
        let mut res = 0;
        for (i, p) in primes.iter().enumerate() {
            let prod = p * p;
            if prod > limit {
                break;
            }
            res += ((limit / prod) + 3) >> 2;
            res -= incl_excl(limit / prod, &primes[i + 1..]);
        }
        res
    }

    let start = std::time::Instant::now();
    let primes = sift(SQRT_N);
    println!("Generated primes up to {SQRT_N}: {:?}", start.elapsed());
    let squarefree = |x| ((x + 3) >> 2) - incl_excl(x, &primes[1..]);
    let res = primes
        .iter()
        .filter(|&&p| p & 3 == 3)
        .fold(squarefree(N), |acc, p| acc + squarefree(N / (p * p)));
    let end = start.elapsed();
    println!("res = {res}, took {end:?}");
}
fn solve_alt() {
    let start = std::time::Instant::now();
    let primes = sift(SQRT_N);
    let xcbrt = (N as f64).cbrt() as u64;
    let large_keys = (1..=xcbrt)
        .step_by(2)
        .map(|d| N / (d * d))
        .collect_vec()
        .into_boxed_slice();
    //let small_keys = (1..=xcbrt).collect_vec().into_boxed_slice();

    let mut large_sqf = vec![0; large_keys.len()].into_boxed_slice();
    let mut small_sqf = vec![0; 1 + ((xcbrt as usize - 1) >> 2)].into_boxed_slice();

    for v in (1..=xcbrt).step_by(4) {
        small_sqf[(v as usize - 1) >> 2] = (v + 3) >> 2;
    }
    for (i, &v) in large_keys.iter().enumerate() {
        large_sqf[i] = (v + 3) >> 2;
    }
    for &p in &primes[1..] {
        for (i, &v) in large_keys.iter().enumerate() {
            if v < p * p {
                break;
            }
            large_sqf[i] -= if v / (p * p) <= xcbrt {
                small_sqf[((v / (p * p)) as usize - 1) >> 2]
            } else {
                large_sqf[((2 * i + 1) * (p as usize)) >> 1]
            };
        }
        for v in (p * p..=(xcbrt & !3) | 1).rev().step_by(4) {
            small_sqf[(v as usize - 1) >> 2] -= small_sqf[((v / (p * p)) as usize - 1) >> 2];
        }
    }
    let mut res = large_sqf[0];
    for p in primes {
        if p & 3 != 3 {
            continue;
        }
        res += if N / (p * p) <= xcbrt {
            small_sqf[((N / (p * p)) as usize - 1) >> 2]
        } else {
            large_sqf[p as usize >> 1]
        };
    }
    let end = start.elapsed();
    println!("res = {res}, took {end:?}");
}
fn solve() {
    let start = std::time::Instant::now();
    let xsqrt = N.isqrt();
    let primes = sift(xsqrt);
    let xcbrt = (N as f64).cbrt() as u64;
    let mut keys = (1..=xcbrt)
        .step_by(2)
        .map(|d| N / (d * d))
        .chain((1..=xcbrt).rev())
        .collect_vec();
    keys.reverse();
    assert!(keys.is_sorted());
    keys.dedup();
    let keys = keys.into_boxed_slice();

    let rank = |e| {
        if e <= xcbrt {
            e as usize - 1
        } else {
            xcbrt as usize + keys[xcbrt as usize..].partition_point(|&v| v < e)
            /* xcbrt as usize
            + binsearch(
                keys.len() - xcbrt as usize,
                |i| keys[xcbrt as usize..][i],
                |v| v < e,
            ) */
        }
    };
    let mut squarefree = vec![0; keys.len()].into_boxed_slice();
    for (i, &v) in keys.iter().enumerate() {
        squarefree[i] = (v + 3) >> 2;
    }
    for &p in &primes[1..] {
        for (i, &v) in keys.iter().enumerate().rev() {
            if v < p * p {
                break;
            }
            squarefree[i] -= squarefree[rank(v / (p * p))];
        }
    }
    let res = primes
        .iter()
        .filter(|&&p| p & 3 == 3)
        .fold(squarefree[keys.len() - 1], |acc, p| {
            acc + squarefree[rank(N / (p * p))]
        });
    let end = start.elapsed();
    println!("res = {res}, took {end:?}");
}
fn solve_alt3() {
    let start = std::time::Instant::now();
    let xsqrt = N.isqrt();
    let primes = sift(xsqrt);
    let xcbrt = (N as f64).cbrt() as u64;
    /* let large_keys = (1..=xcbrt)
    .step_by(2)
    .map(|d| N / (d * d))
    .collect_vec()
    .into_boxed_slice(); */
    //let small_keys = (1..=xcbrt).collect_vec().into_boxed_slice();

    let mut large_sqf = vec![0; (xcbrt + 1) as usize >> 1].into_boxed_slice();
    let mut small_sqf = vec![0; xcbrt as usize].into_boxed_slice();

    for v in 1..=xcbrt {
        small_sqf[v as usize - 1] = (v + 3) >> 2;
    }
    for (i, v) in (1..=xcbrt).step_by(2).map(|d| N / (d * d)).enumerate() {
        large_sqf[i] = (v + 3) >> 2;
    }
    for &p in &primes[1..] {
        let pp = p * p;
        let npp = N / pp;
        large_sqf[0] -= if npp <= xcbrt {
            small_sqf[npp as usize - 1]
        } else {
            large_sqf[p as usize >> 1]
        };
        let mut dp = 3 * p;
        let mut dd = 1;
        // (d+2)^2 = d^2 + 4d+4 = d^2 + 4(d+2 - 1)
        for d in (3..=xcbrt).step_by(2) {
            dd += (d - 1) << 2;
            if d * d > npp {
                break;
            }
            large_sqf[d as usize >> 1] -= if npp / dd <= xcbrt {
                small_sqf[(npp / dd) as usize - 1]
            } else {
                large_sqf[dp as usize >> 1]
            };
            dp += p << 1;
        }
        for v in (pp..=xcbrt).rev() {
            small_sqf[v as usize - 1] -= small_sqf[(v / pp) as usize - 1];
        }
    }
    let mut res = large_sqf[0];
    for p in primes {
        if p & 3 != 3 {
            continue;
        }
        res += if N / (p * p) <= xcbrt {
            small_sqf[(N / (p * p)) as usize - 1]
        } else {
            large_sqf[p as usize >> 1]
        }
    }
    let end = start.elapsed();
    println!("res = {res}, took {end:?}");
}

fn solve_alt2() {
    let start = std::time::Instant::now();
    let xsqrt = N.isqrt();
    let xcbrt = (N as f64).cbrt() as u64;
    let large_keys = (1..=xcbrt)
        .step_by(2)
        .map(|d| N / (d * d))
        .collect_vec()
        .into_boxed_slice();
    //let small_keys = (1..=xcbrt).collect_vec().into_boxed_slice();

    let mut large_sqf = vec![0; large_keys.len()].into_boxed_slice();
    let mut small_sqf = vec![0; xcbrt as usize].into_boxed_slice();

    for v in 1..=xcbrt {
        small_sqf[v as usize - 1] = (v + 3) >> 2;
    }
    for (i, &v) in large_keys.iter().enumerate() {
        large_sqf[i] = (v + 3) >> 2;
    }
    for p in sieve_it().skip(1).take_while(|&p| p <= xsqrt as usize) {
        let p = p as u64;
        for (i, &v) in large_keys.iter().enumerate() {
            if v < p * p {
                break;
            }
            large_sqf[i] -= if v / (p * p) <= xcbrt {
                small_sqf[(v / (p * p)) as usize - 1]
            } else {
                large_sqf[((2 * i + 1) * (p as usize)) >> 1]
            };
        }
        // xcbrt < p^2 => p>x^1/6
        for v in (1..=xcbrt).rev() {
            if v < p * p {
                break;
            }
            small_sqf[v as usize - 1] -= small_sqf[(v / (p * p)) as usize - 1];
        }
    }
    let mut res = large_sqf[0];
    for p in sieve_it().take_while(|&p| p <= xsqrt as usize) {
        let p = p as u64;

        if p & 3 != 3 {
            continue;
        }
        res += if N / (p * p) <= xcbrt {
            small_sqf[(N / (p * p)) as usize - 1]
        } else {
            large_sqf[p as usize >> 1]
        }
    }
    let end = start.elapsed();
    println!("res = {res}, took {end:?}");
}

fn lucy_based_2(x: u64) -> u64 {
    let xsqrt = x.isqrt();
    let primes = sift(xsqrt);
    let xcbrt = (x as f64).cbrt() as u64;
    let mut keys = (1..=xcbrt)
        .map(|d| x / (d * d))
        .chain((1..=xcbrt).rev())
        .collect_vec();
    keys.reverse();
    assert!(keys.is_sorted());
    keys.dedup();
    let keys = keys.into_boxed_slice();
    let mut squarefree = FIArrayU64::new(x);
    for &v in &keys {
        squarefree[v] = (v + 3) >> 2;
    }
    for &p in &primes[1..] {
        for &v in keys.iter().rev() {
            if v < p * p {
                break;
            }
            squarefree[v] -= squarefree[v / (p * p)];
        }
    }
    primes
        .iter()
        .filter(|&&p| p & 3 == 3)
        .fold(squarefree[x], |acc, p| acc + squarefree[x / (p * p)])
}

fn count_sqf() {
    const N: usize = self::N as _;
    const SQRT_N: usize = N.isqrt();
    const fn icbrt(x: usize) -> usize {
        let mut rt = 1 << (1 + x.ilog2().div_ceil(3));
        let mut x_div_rt2 = (x / rt) / rt;
        while rt > x_div_rt2 {
            rt = ((rt << 1) + x_div_rt2) / 3;
            x_div_rt2 = (x / rt) / rt;
        }
        rt
    }
    const B: usize = icbrt(N);
    const A: usize = N / (B * B);

    let start = std::time::Instant::now();
    let primes = sift(SQRT_N as u64);
    dbg!(start.elapsed());
    let mut sqf_small = vec![1; 1 + ((A - 1) >> 2)].into_boxed_slice();
    for &p in &primes[1..] {
        let p = p as usize;
        let pp = p * p;
        if pp > A {
            break;
        }
        for m in (pp..=A).step_by(4 * pp) {
            sqf_small[(m - 1) >> 2] = 0;
        }
    }
    for i in 1..sqf_small.len() {
        sqf_small[i] += sqf_small[i - 1];
    }
    let sqrts = (1..=A)
        .step_by(4)
        .map(|i| (N / i).isqrt())
        .collect_vec()
        .into_boxed_slice(); // precompute once, used very often

    let mut sqf_big = vec![0; B >> 1].into_boxed_slice(); // indexed by denominator
    for d in (1..B).step_by(2).rev() {
        let v = N / (d * d);
        let b = icbrt(v);
        let a = v / (b * b);

        let mut sqf =
            ((v + 3) >> 2) + sqf_small[(a - 1) >> 2] * ((b + 1) >> 1) - (((SQRT_N / d) + 1) >> 1);

        /* for i in (5..=a).step_by(4) {
            if sqf_small[(i - 1) >> 2] != sqf_small[(i - 2) >> 2] {
                sqf -= ((sqrts[(i - 1) >> 2] / d) + 1) >> 1;
            }
        } */
        for i in 1..=(a - 1) >> 2 {
            if sqf_small[i] != sqf_small[i - 1] {
                sqf -= ((sqrts[i] / d) + 1) >> 1;
            }
        }
        for i in (3..=b).step_by(2) {
            sqf -= if i * d < B {
                sqf_big[(i * d) >> 1]
            } else {
                sqf_small[((v / (i * i)) - 1) >> 2]
            };
        }
        sqf_big[d >> 1] = sqf;
    }

    let mut res = sqf_big[0];
    let mut primes = primes;
    primes.retain(|&p| p & 3 == 3);
    let primes = primes;

    for &p in &primes {
        let p = p as usize;
        if p >= B {
            break;
        }
        res += sqf_big[p >> 1];
    }
    for (i, &v) in sqrts.iter().enumerate() {
        if i == 0 {
            res += primes.len();
            continue;
        }
        if sqf_small[i] != sqf_small[i - 1] {
            res += primes.partition_point(|&p| p <= v as u64);
        }
    }

    res -= sqf_small[(A - 1) >> 2] * primes.partition_point(|&p| p < B as u64);
    /* for &p in &primes[1..] {
        let p = p as usize;
        if p & 3 != 3 {
            continue;
        }
        res += if p < B {
            sqf_big[p >> 1]
        } else {
            sqf_small[((N / (p * p)) - 1) >> 2]
        };
    } */

    let end = start.elapsed();
    println!("res = {res}, took {end:?}");
}
