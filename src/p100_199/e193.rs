use std::time::Instant;

use itertools::Itertools;

use crate::utils::{
    FIArray::FIArrayU64, bit_array::BitArray, multiplicative_function_summation::mobius_sieve,
    sieve_of_pritchard::sift,
};

#[must_use]
fn count_squarefree(limit: u64) -> u64 {
    fn incl_excl(limit: u64, primes: &[u64]) -> u64 {
        let mut res = 0;
        for (i, &p) in primes.iter().enumerate() {
            let prod = p * p;
            if prod > limit {
                break;
            }
            res += limit / prod;
            res -= incl_excl(limit / prod, &primes[i + 1..]);
        }
        res
    }
    let first_primes = sift(limit.isqrt()).into_boxed_slice();
    limit - incl_excl(limit, &first_primes)
}

pub fn main() {
    const N: i64 = 1e18 as _;
    let start = Instant::now();
    let res = opt(N as _);
    let end = start.elapsed();
    println!("res = {res}, took {end:?}");
    let start = Instant::now();
    let res = opt2(N as _);
    let end = start.elapsed();
    println!("res = {res}, took {end:?}");
    let start = Instant::now();
    let res = dirimul_opt(N as _);
    let end = start.elapsed();
    println!("res = {res}, took {end:?}");
    let start = Instant::now();
    let res = dirichlet_mul_based_opt(N as _);
    let end = start.elapsed();
    println!("res = {res}, took {end:?}");
    let start = Instant::now();
    let res = count_squarefree(N as _);
    let end = start.elapsed();
    println!("res = {res}, took {end:?}");
    let start = Instant::now();
    let res = count_squarefree_mobius_alt(N as _);
    let end = start.elapsed();
    println!("res = {res}, took {end:?}");
    let start = Instant::now();
    let res = count_squarefree_mobius(N as _);
    let end = start.elapsed();
    println!("res = {res}, took {end:?}");
    /* let start = Instant::now();
    let res = dirichlet_mul_based(N as _);
    let end = start.elapsed();
    println!("res = {res}, took {end:?}"); */
}
#[must_use]
fn count_squarefree_mobius(n: usize) -> usize {
    let m = mobius_sieve(n.isqrt() + 1).into_boxed_slice();
    (1..=n.isqrt())
        .map(|d| i64::from(m[d]) * (n / (d * d)) as i64)
        .sum::<i64>() as usize
}

#[must_use]
fn count_squarefree_mobius_alt(limit: usize) -> usize {
    let n = limit.isqrt() + 1;
    let mut sum = limit as i64;
    let mut res = vec![0i8; n];
    let mut composite = BitArray::zeroed(n);
    let mut primes = vec![];
    res[1] = 1;
    for i in 2..n {
        if !composite.get(i) {
            primes.push(i);
            res[i] = -1;
            sum -= (limit / (i * i)) as i64;
        }
        for &p in &primes {
            let ip = i * p;
            if ip >= n {
                break;
            }
            composite.set(ip);
            if i % p == 0 {
                res[ip] = 0;
                break;
            }
            res[ip] = -res[i];
            sum += i64::from(res[ip]) * (limit / (ip * ip)) as i64;
        }
    }
    sum as usize
}
// computation of zeta(2)/zeta(2s)
#[must_use]
fn dirichlet_mul_based_opt(x: u64) -> u64 {
    let xsqrt = x.isqrt();
    let primes = sift(xsqrt).into_boxed_slice();
    let xcbrt = (x as f64).cbrt() as u64;
    let mut keys = (1..xcbrt)
        .map(|d| x / (d * d))
        .chain((xcbrt != x / (xcbrt * xcbrt)).then_some(x / (xcbrt * xcbrt)))
        .chain((1..=xcbrt).rev())
        .collect_vec(); // can probably reduce the number of keys
    keys.reverse();
    assert!(keys.is_sorted());
    //dbg!(keys.len());
    keys.dedup();
    let keys = keys.into_boxed_slice();
    let rank = |e| {
        if e <= xcbrt {
            e as usize - 1
        } else {
            xcbrt as usize + keys[xcbrt as usize..].partition_point(|&v| v < e)
        }
    };
    let mut s = keys.clone();
    let lim = primes.partition_point(|&p| p <= xsqrt.isqrt());
    for &p in &primes[..lim] {
        for (i, &v) in keys.iter().enumerate().rev() {
            if v < p * p {
                break;
            }
            s[i] -= s[rank(v / (p * p))];
        }
    }
    for &p in &primes[lim..] {
        s[keys.len() - 1] -= s[rank(x / (p * p))];
    }
    s[keys.len() - 1] //- s[keys.len() - 2]
}
fn dirimul_opt(N: u64) -> u64 {
    let xsqrt = N.isqrt();
    let primes = sift(xsqrt);
    let xcbrt = (N as f64).cbrt() as u64;
    let large_keys = (1..xcbrt)
        .map(|d| N / (d * d))
        .chain((xcbrt != N / (xcbrt * xcbrt)).then_some(N / (xcbrt * xcbrt)))
        .collect_vec()
        .into_boxed_slice();
    let mut large_sqf = vec![0; large_keys.len()].into_boxed_slice();
    let mut small_sqf = vec![0; xcbrt as usize].into_boxed_slice();

    for v in 1..=xcbrt {
        small_sqf[v as usize - 1] = v;
    }
    for (i, &v) in large_keys.iter().enumerate() {
        large_sqf[i] = v;
    }
    let lim = primes.partition_point(|&p| p <= xsqrt.isqrt());
    for &p in &primes[..lim] {
        for (i, &v) in large_keys.iter().enumerate() {
            if v < p * p {
                break;
            }
            large_sqf[i] -= if v / (p * p) <= xcbrt {
                small_sqf[(v / (p * p)) as usize - 1]
            } else {
                large_sqf[((i + 1) * (p as usize)) - 1]
            };
        }
        for v in (1..=xcbrt).rev() {
            if v < p * p {
                break;
            }
            small_sqf[v as usize - 1] -= small_sqf[(v / (p * p)) as usize - 1];
        }
    }
    for &p in &primes[lim..] {
        large_sqf[0] -= if N / (p * p) <= xcbrt {
            small_sqf[(N / (p * p)) as usize - 1]
        } else {
            large_sqf[p as usize - 1]
        };
    }
    large_sqf[0]
}
// computation of zeta(2)/zeta(2s)
#[must_use]
fn dirichlet_mul_based(x: u64) -> FIArrayU64 {
    let mut s = FIArrayU64::unit(x);
    let keys = FIArrayU64::keys(x).collect_vec().into_boxed_slice();
    let primes = sift(x.isqrt()).into_boxed_slice();

    for p in primes {
        for (i, &v) in keys.iter().enumerate().rev() {
            if v < p * p {
                break;
            }
            s.arr[i] -= s[v / (p * p)];
        }
    }
    s
}

// https://arxiv.org/pdf/1107.4890
// essentially just leverages the dirichlet hyperbola method
// becomes faster than incl-excl somewhere between 2^34 and 2^35
// O(n^0.4) space and time
#[must_use]
fn opt(x: usize) -> usize {
    fn mobius_sieve(n: usize) -> Vec<i64> {
        let mut res = vec![0; n];
        if n < 2 {
            return res;
        }
        let mut composite = BitArray::zeroed(n);
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

    let I = (x as f64).powf(0.2) as usize;
    let D = (x / I).isqrt();
    dbg!(I, D);
    let mut mertens_small = mobius_sieve(D + 1);
    let mut s1 = 0;
    for d in 1..=D {
        s1 += mertens_small[d] * (x / (d * d)) as i64;
        mertens_small[d] += mertens_small[d - 1];
    }

    let mut mertens_big = vec![0; I].into_boxed_slice(); // indexed by denominator
    for i in (1..I).rev() {
        let v = (x / i).isqrt();
        let vsqrt = v.isqrt();
        let mut m = 1 - v as i64 + vsqrt as i64 * mertens_small[vsqrt];
        for d in 2..=vsqrt {
            m -= if v / d <= D {
                mertens_small[v / d]
            } else {
                mertens_big[i * d * d]
            };
            m -= (mertens_small[d] - mertens_small[d - 1]) * (v / d) as i64;
        }
        mertens_big[i] = m;
    }
    let s2 = mertens_big.iter().sum::<i64>() - (I - 1) as i64 * mertens_small[D];

    (s1 + s2) as usize
}
// O(n^0.4375) time, O(n^0.375) space
#[must_use]
fn opt2(x: usize) -> usize {
    fn mobius_sieve(n: usize) -> Vec<i64> {
        let mut res = vec![0; n];
        if n < 2 {
            return res;
        }
        let mut composite = BitArray::zeroed(n);
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

    let I = (x as f64).powf(1. / 4.) as usize;
    let D = (x / I).isqrt();
    dbg!(I, D);
    let mut mertens_small = mobius_sieve(D + 1);
    let mut s1 = 0;
    for d in 1..=D {
        s1 += mertens_small[d] * (x / (d * d)) as i64;
        mertens_small[d] += mertens_small[d - 1];
    }

    let mut mertens_big = vec![0; I].into_boxed_slice(); // indexed by denominator
    for i in (1..I).rev() {
        let v = (x / i).isqrt();
        let vsqrt = v.isqrt();
        let mut m = 1 - v as i64 + vsqrt as i64 * mertens_small[vsqrt];
        for d in 2..=vsqrt {
            m -= if v / d <= D {
                mertens_small[v / d]
            } else {
                mertens_big[i * d * d]
            };
            m -= (mertens_small[d] - mertens_small[d - 1]) * (v / d) as i64;
        }
        mertens_big[i] = m;
    }
    let s2 = mertens_big.iter().sum::<i64>() - (I - 1) as i64 * mertens_small[D];

    (s1 + s2) as usize
}
