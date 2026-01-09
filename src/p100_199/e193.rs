use std::time::Instant;

use itertools::Itertools;

use crate::utils::{
    FIArray::FIArrayU64,
    bit_array::BitArray,
    multiplicative_function_summation::{self, mobius_sieve},
    primes::prime_sieves::sift,
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
    const N: u64 = 1 << 50;
    /*let start = Instant::now();
    let res = opt_blocked(N as _);
    let end = start.elapsed();
    println!("res = {res}, took {end:?}");*/
    let start = Instant::now();
    let res = opt(N as _);
    let end = start.elapsed();
    println!("res = {res}, took {end:?}");
    let start = Instant::now();
    let res = opt2(N as _);
    let end = start.elapsed();
    println!("res = {res}, took {end:?}");
    let start = Instant::now();
    let res = count_sqf(N as _);
    let end = start.elapsed();
    println!("res = {res}, took {end:?}");
    let start = Instant::now();
    let res = opt3(N as _);
    let end = start.elapsed();
    println!("res = {res}, took {end:?}");

    let start = Instant::now();
    let res = multiplicative_function_summation::count_squarefree(N as _)[N as _];
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
        large_sqf[0] -= if N / (p * p) <= xcbrt {
            small_sqf[(N / (p * p)) as usize - 1]
        } else {
            large_sqf[(p as usize) - 1]
        };
        for (i, &v) in large_keys.iter().enumerate().skip(p as usize) {
            if v < p * p {
                break;
            }
            large_sqf[i] -= if v / (p * p) <= xcbrt {
                small_sqf[(v / (p * p)) as usize - 1]
            } else {
                large_sqf[((i + 1) * (p as usize)) - 1]
            };
        }
        for v in (p * p..=xcbrt).rev() {
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

// O(n^4/9) time, O(n^1/3) space
fn count_sqf(x: usize) -> usize {
    const fn icbrt(x: usize) -> usize {
        let mut rt = 1 << (1 + x.ilog2().div_ceil(3));
        let mut x_div_rt2 = (x / rt) / rt;
        while rt > x_div_rt2 {
            rt = ((rt << 1) + x_div_rt2) / 3;
            x_div_rt2 = (x / rt) / rt;
        }
        rt
    }
    fn sqf_sieve(n: usize) -> Vec<usize> {
        unsafe { core::hint::assert_unchecked(n >= 1) };
        let mut sqf = vec![1; n];
        sqf[0] = 0;
        if n < 2 {
            return sqf;
        }
        let sqrtn = n.isqrt();
        let mut d2 = 1;
        for d in 2..=sqrtn {
            d2 += (d << 1) - 1;
            if sqf[d2] == 0 {
                continue;
            }
            for m in (d2..n).step_by(d2) {
                sqf[m] = 0;
            }
        }
        sqf
    }
    let B = icbrt(x);
    let A = x / (B * B);
    dbg!(A, B);
    let xsqrt = x.isqrt();

    let mut sqf_small = sqf_sieve(A + 1);
    for i in 2..=A {
        sqf_small[i] += sqf_small[i - 1];
    }
    let sqrts = (2..=A)
        .map(|i| (x / i).isqrt())
        .collect_vec()
        .into_boxed_slice(); // precompute once, used very often
    let mut sqf_big = vec![0; B - 1].into_boxed_slice(); // indexed by denominator
    for d in (1..B).rev() {
        let v = x / (d * d);
        let b = icbrt(v);
        let a = v / (b * b);

        let mut sqf = v + sqf_small[a] * b - xsqrt / d;
        for i in 2..=a {
            sqf -= (sqf_small[i] - sqf_small[i - 1]) * sqrts[i - 2] / d; //(v / i).isqrt();
        }
        for i in 2..=b {
            sqf -= if i * d < B {
                sqf_big[(i * d) - 1]
            } else {
                sqf_small[v / (i * i)]
            };
        }
        sqf_big[d - 1] = sqf;
    }
    sqf_big[0]
}

// https://arxiv.org/pdf/1107.4890
// essentially just leverages the dirichlet hyperbola method:
// sqf = mu_sqrt * u, \alpha = n^(4/5), \beta = n^(1/5)
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
                if i.is_multiple_of(p) {
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
    //dbg!(I, D);
    let mut mertens_small = mobius_sieve(D + 1).into_boxed_slice();
    let mut s1 = x as i64;
    for d in 2..=D {
        s1 += mertens_small[d] * (x / (d * d)) as i64;
        mertens_small[d] += mertens_small[d - 1];
    }
    let mut mertens_big = vec![0; I - 1]; // indexed by denominator
    for i in (1..I).rev() {
        let v = (x / i).isqrt();

        let vsqrt = v.isqrt();
        let mut m = 1 - v as i64 + vsqrt as i64 * mertens_small[vsqrt];
        for d in 2..=vsqrt {
            m -= if v / d <= D {
                mertens_small[v / d]
            } else {
                mertens_big[i * d * d - 1]
            };

            m -= (mertens_small[d] - mertens_small[d - 1]) * (v / d) as i64;
        }
        mertens_big[i - 1] = m;
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
    //dbg!(I, D);
    let mut mertens_small = mobius_sieve(D + 1);
    let mut s1 = 0;
    for d in 1..=D {
        s1 += mertens_small[d] * (x / (d * d)) as i64;
        mertens_small[d] += mertens_small[d - 1];
    }

    let mut mertens_big = vec![0; I - 1].into_boxed_slice(); // indexed by denominator
    for i in (1..I).rev() {
        let v = (x / i).isqrt();
        let vsqrt = v.isqrt();
        let mut m = 1 - v as i64 + vsqrt as i64 * mertens_small[vsqrt];
        for d in 2..=vsqrt {
            m -= if v / d <= D {
                mertens_small[v / d]
            } else {
                mertens_big[i * d * d - 1]
            };
            m -= (mertens_small[d] - mertens_small[d - 1]) * (v / d) as i64;
        }
        mertens_big[i - 1] = m;
    }
    let s2 = mertens_big.iter().sum::<i64>() - (I - 1) as i64 * mertens_small[D];

    (s1 + s2) as usize
}

// O(n^0.5) time, O(n^(1 / 3)) space
#[must_use]
fn opt3(x: usize) -> usize {
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
    const fn icbrt(x: usize) -> usize {
        let mut rt = 1 << x.ilog2().div_ceil(3);
        let mut x_div_rt2 = (x / rt) / rt;
        while rt > x_div_rt2 {
            rt = ((rt << 1) + x_div_rt2) / 3;
            x_div_rt2 = (x / rt) / rt;
        }
        rt
    }
    let I = icbrt(x >> 2);
    let D = (x / I).isqrt();
    let mut mertens_small = mobius_sieve(D + 1);
    let mut s1 = 0;
    for d in 1..=D {
        s1 += mertens_small[d] * (x / (d * d)) as i64;
        mertens_small[d] += mertens_small[d - 1];
    }

    let mut mertens_big = vec![0; I - 1].into_boxed_slice(); // indexed by denominator
    for i in (1..I).rev() {
        let v = (x / i).isqrt();
        let vsqrt = v.isqrt();
        let mut m = 1 - v as i64 + vsqrt as i64 * mertens_small[vsqrt];
        for d in 2..=vsqrt {
            m -= if v / d <= D {
                mertens_small[v / d]
            } else {
                mertens_big[i * d * d - 1]
            };
            m -= (mertens_small[d] - mertens_small[d - 1]) * (v / d) as i64;
        }
        mertens_big[i - 1] = m;
    }
    let s2 = mertens_big.iter().sum::<i64>() - (I - 1) as i64 * mertens_small[D];

    (s1 + s2) as usize
}

// idk it nearly works, can't figure out the issue
fn opt_blocked(x: usize) -> usize {
    let I = ((x as f64) * (x as f64).ln().ln().powi(4)).powf(0.2) as usize;
    let D = (x / I).isqrt();
    let B = D.isqrt();
    let L = D.div_ceil(B);
    //dbg!(I, D, B, L);
    let mut ilist = vec![const { vec![] }; L];
    let mut min_k = vec![1; I - 1];
    let mut Mx = vec![1; I - 1]; // indexed by denominator
    ilist[0].extend(1..I);

    let mut MxBlockUpdate = |a, b, i: usize, mut k, M: &[i64]| {
        let x_i = (x / i).isqrt();
        let mut d_a = x_i / k;
        loop {
            let d_b = x_i / (k + 1);
            Mx[i - 1] -= (d_a as i64 - d_b as i64) * M[k - a];
            k = x_i / d_b;
            d_a = d_b;
            if k > b {
                return k;
            }
        }
    };

    let mut s1 = 0;
    let mut mertens_D = 0;

    let primes = sift(B as u64);
    let mut mobius = vec![0; B + 1];
    let mut factor = vec![1; B + 1];
    for l in 0..L {
        mobius.fill(1);
        let a_l = l * B;
        let a_l1 = D.min(a_l + B);
        for k in a_l + 1..=a_l1 {
            factor[k - a_l] = k;
        }
        mobius[0] = mertens_D;
        for &p in &primes {
            let p = p as usize;
            let multiple = (a_l + 1).next_multiple_of(p * p);
            for k in (multiple..=a_l1).step_by(p * p) {
                mobius[k - a_l] = 0;
            }
            let multiple = (a_l + 1).next_multiple_of(p);
            for k in (multiple..=a_l1).step_by(p) {
                mobius[k - a_l] *= -1;
                factor[k - a_l] /= p;
            }
        }
        //assert_eq!(mobius[0], mertens_D);
        //assert_eq!(factor[0], 1);
        for k in a_l + 1..=a_l1 {
            if factor[k - a_l] != 1 {
                mobius[k - a_l] *= -1;
            }
            s1 += mobius[k - a_l] * (x / (k * k)) as i64;
            mobius[k - a_l] += mobius[k - a_l - 1];
        }
        for j in 0..ilist[l].len() {
            let i = ilist[l][j];
            min_k[i - 1] = MxBlockUpdate(a_l, a_l1, i, min_k[i - 1], &mobius);
            let l_prime = min_k[i - 1] / B;
            if l_prime < L && min_k[i - 1] < (x / i).isqrt() {
                ilist[l_prime].push(i);
            }
        }
        ilist[l].clear();
        mertens_D = mobius[a_l1 - a_l];
    }
    //assert!(ilist.iter().all(Vec::is_empty));
    for i in (1..I).rev() {
        let x_i = (x / i).isqrt();
        for d in 2.. {
            if x_i / d <= D {
                break;
            }
            Mx[i - 1] -= Mx[i * d * d - 1];
        }
    }

    let s2 = Mx.iter().sum::<i64>() - (I - 1) as i64 * mertens_D;
    //dbg!(s1, s2);

    (s1 + s2) as usize
}
