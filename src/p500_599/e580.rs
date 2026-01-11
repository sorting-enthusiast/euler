use itertools::Itertools;

use crate::utils::math::iroot;

// Hilbert numbers are numbers that are products of any number of primes equivalent to 1 mod 4,
// and an even number of primes equivalent to 3 mod 4

// a hilbert number is h-squarefree if it is either squarefree, or if only 1 of its prime factors = 3 mod 4 appears with multiplicity 1< <4
const N: usize = 1e16 as _;
const SQRT_N: usize = N.isqrt();
pub fn main() {
    println!("res = {}", 2_327_213_148_095_366_u64);
    solve();
}
/// O(n^1/2) time, O(n^1/3) space
/// Using the identity (n+3)/4 = sum_{d<=n^1/2, d odd} SQF'(n/d^2),
/// and the dirichlet hyperbola, we can compute all values of SQF'(n/d^2) in O(n^4/9) time, O(n^1/3) space.
/// The answer to the problem is then SQF'(n) + sum_{p<=n^1/2, p=3 mod 4} SQF'(n/p^2)
/// sum_{p<=n^1/2, p=3 mod 4} SQF'(n/p^2) can be computed with the dirichlet hyperbola method again.
/// Computing the # of primes = 3 mod 4 below each necessary threshold takes O(n^1/2) time and O(n^1/3) space.
fn solve() {
    const B: usize = iroot::<3>(N);
    const A: usize = N / (B * B);

    let start = std::time::Instant::now();
    let mut sqf_small = vec![1; (A + 3) >> 2].into_boxed_slice();
    for p in (3..=A.isqrt()).step_by(2) {
        let pp = p * p;
        if pp > A {
            break;
        }
        if sqf_small[(pp - 1) >> 2] == 0 {
            continue;
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
    let mut sqf_big = vec![0; (B + 1) >> 1].into_boxed_slice(); // indexed by denominator
    for d in (1..B + 1).step_by(2).rev() {
        let v = N / (d * d);
        let b = iroot::<3>(v);
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
            sqf -= if i * d <= B {
                sqf_big[(i * d) >> 1]
            } else {
                sqf_small[((v / (i * i)) - 1) >> 2]
            };
        }
        sqf_big[d >> 1] = sqf;
    }

    // Could try replacing this section with code computing prime counts and the sum of \chi_4 over primes,
    // though I don't see that improving the time complexity at all
    let mut pi1_big = sqrts.clone();
    let mut pi3_big = sqrts.clone();
    let mut pi1_small = vec![0; (B + 3) >> 2].into_boxed_slice();
    let mut pi3_small = vec![0; (B + 1) >> 2].into_boxed_slice();
    for (i, &v) in sqrts.iter().enumerate() {
        pi1_big[i] = ((v + 3) >> 2) - 1;
        pi3_big[i] = (v + 1) >> 2;
    }
    for i in (1..=B).step_by(2) {
        if i & 3 == 1 {
            pi1_small[(i - 1) >> 2] = ((i + 3) >> 2) - 1;
        } else {
            pi3_small[(i - 3) >> 2] = (i + 1) >> 2;
        }
    }
    for p in (3..=SQRT_N.isqrt()).step_by(2) {
        if p > 3 {
            if p & 3 == 1 {
                if pi1_small[(p - 1 - 1) >> 2] == pi1_small[(p - 1) >> 2] {
                    continue;
                }
            } else if pi3_small[(p - 1 - 3) >> 2] == pi3_small[(p - 3) >> 2] {
                continue;
            }
        }
        let pp = p * p;
        let k = (pp - 1) >> 2;
        let mut ind = k;
        if p & 3 == 1 {
            let sp1 = pi1_small[(p - 1 - 1) >> 2];
            let sp3 = pi3_small[(p - 1 - 3) >> 2];
            for (i, &v) in sqrts.iter().enumerate() {
                // v = sqrt(N / (4i+1))
                if v < pp {
                    break;
                }
                pi1_big[i] -= if ind < sqrts.len() {
                    pi1_big[ind]
                } else {
                    pi1_small[((v / p) - 1) >> 2]
                } - sp1;
                pi3_big[i] -= if ind < sqrts.len() {
                    pi3_big[ind]
                } else {
                    pi3_small[((v / p) - 3) >> 2]
                } - sp3;
                ind += (k << 2) | 1;
            }
            for i in (p * p..B + 1).step_by(2).rev() {
                if i & 3 == 1 {
                    pi1_small[(i - 1) >> 2] -= pi1_small[((i / p) - 1) >> 2] - sp1;
                } else {
                    pi3_small[(i - 3) >> 2] -= pi3_small[((i / p) - 3) >> 2] - sp3;
                }
            }
        } else {
            let sp1 = if p > 3 {
                pi3_small[(p - 1 - 3) >> 2]
            } else {
                0
            };
            let sp3 = pi1_small[(p - 1 - 1) >> 2];
            for (i, &v) in sqrts.iter().enumerate() {
                // v = sqrt(N / (4i+1))
                if v < pp {
                    break;
                }
                pi1_big[i] -= if ind < sqrts.len() {
                    pi3_big[ind]
                } else {
                    pi3_small[((v / p) - 3) >> 2]
                } - sp1;
                pi3_big[i] -= if ind < sqrts.len() {
                    pi1_big[ind]
                } else {
                    pi1_small[((v / p) - 1) >> 2]
                } - sp3;
                ind += (k << 2) | 1;
            }
            for i in (p * p..B + 1).step_by(2).rev() {
                if i & 3 == 1 {
                    pi1_small[(i - 1) >> 2] -= pi3_small[((i / p) - 3) >> 2] - sp1;
                } else {
                    pi3_small[(i - 3) >> 2] -= pi1_small[((i / p) - 1) >> 2] - sp3;
                }
            }
        }
    }

    let mut res = sqf_big[0] + pi3_big[0] + sqf_big[1];
    for i in 1..(B + 1) >> 2 {
        if pi3_small[i] != pi3_small[i - 1] {
            res += sqf_big[(i << 1) | 1];
        }
    }
    for i in 1..(A + 3) >> 2 {
        if sqf_small[i] != sqf_small[i - 1] {
            res += pi3_big[i];
        }
    }
    res -= sqf_small[((A + 3) >> 2) - 1] * pi3_small[((B + 1) >> 2) - 1];

    let end = start.elapsed();
    println!("res = {res}, took {end:?}");
}
/*
use crate::utils::{
    FIArray::FIArrayU64,
    multiplicative_function_summation::mobius_sieve,
    primes::prime_sieves::{sieve_it, sift},
};
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
fn solve_alt1() {
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
 */
fn test() {
    let mu = crate::utils::multiplicative_function_summation::mobius_sieve(SQRT_N as usize + 1);
    let mut sqf = (N + 3) >> 2;
    for d in (3..=SQRT_N).step_by(2) {
        if mu[d as usize] == -1 {
            sqf -= (N / (d * d) + 3) >> 2;
        } else if mu[d as usize] == 1 {
            sqf += (N / (d * d) + 3) >> 2;
        }
    }
    dbg!(sqf);
}
// can try reordering the convolutions, computing first 1_{p;4,3}^sqrt * (odd \mu^sqrt), and then convolving with (is_hilbert)
// what I currently do is first (odd \mu^sqrt) * (is_hilbert), i.e. count squarefree hilbert numbers, and then convolve with 1_{p;4,3}^sqrt
// the change would reduce the space complexity by a constant factor, but increase the time complexity a bit
