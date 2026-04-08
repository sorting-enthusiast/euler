#![warn(clippy::pedantic)]
#![allow(long_running_const_eval)]
#![allow(dead_code)]
#![allow(clippy::upper_case_acronyms)]
#![allow(clippy::too_many_lines)]
#![allow(clippy::cast_possible_truncation)]
#![allow(clippy::inline_always)]
#![allow(clippy::cast_sign_loss)]
#![allow(clippy::cast_precision_loss)]
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]
#![allow(clippy::large_stack_arrays)]
use chrono::Local;
use itertools::Itertools;

use crate::{
    fenwick_holes_test::{log_zeta_3, log_zeta_3_odd, log_zeta_fast_alt_2},
    incremental_flattening::{DynamicPrefixSumI64, count_squarefree_partial_flatten},
    p300_399::e362::{mult, mult_sparse},
    utils::{
        FIArray::{DirichletFenwick, DirichletFenwickI64, FIArray, FIArrayI64},
        fast_divisor_sums::{self, divisor_summatory},
        math::{iroot, sum_geometric_mod},
        multiplicative_function_summation::{
            count_squarefree, divisor_sieve, mertens, sqf, sqf_icy,
        },
        polymul::NTT,
        primes::{
            log_zeta::log_zeta_2,
            primecount::{lucy_fenwick, mertens_min25},
            wheel_sieve,
        },
    },
};
pub mod aebp;
pub mod fenwick_holes_test;
mod fib;
pub mod incremental_flattening;
pub mod p0_99;
pub mod p100_199;
pub mod p200_299;
pub mod p300_399;
pub mod p400_499;
pub mod p500_599;
pub mod p600_699;
pub mod p700_799;
pub mod p800_899;
pub mod p900_999;
pub mod test;
pub mod test2;
pub mod utils;
// digital root of n is just n mod 9 if n mod 9 != 0, otherwise 9
const fn is_target_little_endian() -> bool {
    u16::from_ne_bytes([1, 0]) == 1
}
// TODO: understand convex hull based lattice point counting, adapt https://github.com/dengtesla/acm/blob/master/acm%E6%A8%A1%E6%9D%BF/min25_new.cpp
// port icy's lattice point summer, optimize O(n^3/5) FIArray construction for the divisor function and O(n^5/8) prime counting
// Refactor utilities related to prime counting / multiplicative function summation / PET & IPET, + figure out exactly why/how current impl's of mult/mult_powerful work
// TODO: try alternate hybrid repr for FIArray: instead of partially flattened fenwick tree over all keys, use partially flattened fenwick tree for large keys and partially prefix-summed flat array for small keys
// might actually be what negiizhao meant
pub fn main() {
    const { assert!(is_target_little_endian()) }; // some code relies on this
    println!("Started running at: {} ", Local::now().time());
    //dbg!((10usize ^ 7 ^ 3).count_ones());
    p400_499::e487::main();
    p500_599::e556::main();
    p900_999::e989::main();
    p300_399::e362::main();

    //1.8656068s - 2^40, 61.5886751s - 2^48
    const N: usize = 1e13 as _;
    incremental_flattening::main();
    fenwick_holes_test::main();

    let start = std::time::Instant::now();
    let p = log_zeta_fast(N)[N];
    let end = start.elapsed();
    println!("\"fast\" prime counting output for {N}: {p} | {end:?}");

    let start = std::time::Instant::now();
    let p = log_zeta_fast_alt(N)[N];
    let end = start.elapsed();
    println!("\"fast\" prime counting output for {N}: {p} | {end:?}");
    assert_eq!(log_zeta_fast_alt_2(N), log_zeta_2(N));

    println!("counting sqf");
    let start = std::time::Instant::now();
    let s1 = count_squarefree(N);
    let end = start.elapsed();
    dbg!(end, s1[N]);
    let start = std::time::Instant::now();
    let s1 = count_squarefree_partial_flatten(N);
    let end = start.elapsed();
    dbg!(end, s1[N]);
    let start = std::time::Instant::now();
    let s1 = sqf(N);
    let end = start.elapsed();
    dbg!(end, s1[N]);
    let start = std::time::Instant::now();
    let s2 = sqf_icy(N);
    let end = start.elapsed();
    dbg!(end, s2[N]);
    assert_eq!(s1, s2);

    println!("summing mobius");
    let start = std::time::Instant::now();
    let s1 = mertens(N);
    let end = start.elapsed();
    println!("griff's blog: {end:?}, {}", s1[N]);

    let start = std::time::Instant::now();
    let s2 = div_i64(&FIArrayI64::eps(N), &FIArrayI64::unit(N));
    let end = start.elapsed();
    println!("Dirichlet division: {end:?}, {}", s2[N]);
    /*let start = std::time::Instant::now();
    let s2 = inv_i64(&FIArrayI64::unit(N));
    let end = start.elapsed();
    dbg!(end, s2[N]);
    assert_eq!(s1, s2);*/
    let start = std::time::Instant::now();
    let mut zeta = DirichletFenwickI64::zeta(N);
    let lim = iroot::<7>(N) + 1;
    let mut primes = vec![];
    println!("started removal of primes < {lim}: {:?}", start.elapsed());
    for p in 2..lim {
        if zeta.get_bucket_prefix(p - 1) == 1 {
            continue;
        }
        primes.push(p);
        zeta.sparse_mul_at_most_one(p, 1);
    }
    let zeta_lim = FIArrayI64::from(zeta);
    println!(
        "Finished removal of primes < {lim}, started convolution: {:?}",
        start.elapsed()
    );
    let mut zeta_2 = DirichletFenwickI64::from(div_i64(&FIArrayI64::eps(N), &zeta_lim));
    println!(
        "Finished convolution, started adding back primes < {lim}: {:?}",
        start.elapsed()
    );
    for &p in primes.iter().rev() {
        zeta_2.sparse_mul_at_most_one(p, 1);
    }
    let approx = FIArrayI64::from(zeta_2);
    println!(
        "Finished adding back primes < {lim}, started correction: {:?}",
        start.elapsed()
    );
    let end = start.elapsed();
    println!("Sparse division (fenwick): {end:?}, {}", approx[N]);

    let start = std::time::Instant::now();
    let mut zeta = FIArrayI64::unit(N);
    let lim = iroot::<7>(N) + 1;
    let mut primes = vec![];
    println!("started removal of primes < {lim}: {:?}", start.elapsed());
    let keys = FIArray::keys(N).collect_vec();

    for p in 2..lim {
        if zeta.arr[p - 1] == 1 {
            continue;
        }
        primes.push(p);
        for (i, &v) in keys.iter().enumerate().rev() {
            if v < p {
                break;
            }
            zeta.arr[i] -= zeta[v / p];
        }
    }
    println!(
        "Finished removal of primes < {lim}, started convolution: {:?}",
        start.elapsed()
    );
    let mut mu_lim = div_i64(&FIArrayI64::eps(N), &zeta);
    println!(
        "Finished convolution, started adding back primes < {lim}: {:?}",
        start.elapsed()
    );
    for &p in primes.iter().rev() {
        for (i, &v) in keys.iter().enumerate().rev() {
            if v < p {
                break;
            }
            mu_lim.arr[i] -= mu_lim[v / p];
        }
    }
    println!(
        "Finished adding back primes < {lim}, started correction: {:?}",
        start.elapsed()
    );
    let end = start.elapsed();
    println!("Sparse division (naive): {end:?}, {}", mu_lim[N]);

    let start = std::time::Instant::now();
    let s1 = mertens_min25(N);
    let end = start.elapsed();
    println!("Min_25 sieve based: {end:?}, {}", s1[N]);
    let start = std::time::Instant::now();
    let mut pi = inverse_pseudo_euler_transform_fraction_i64(FIArrayI64::unit(N));
    let mut primes = vec![];
    for p in 2..=pi.isqrt {
        if pi.arr[p - 1] != pi.arr[p - 2] {
            primes.push(p);
        }
    }
    for e in &mut pi.arr {
        *e *= -1;
    }

    let s1 = mult_correction_i64(
        &pseudo_euler_transform_fraction_i64(pi),
        &primes,
        |_pp, _p, e| -i64::from(e < 2),
    );
    let end = start.elapsed();
    println!("Pseudo-euler transform based: {end:?}, {}", s1[N]);

    let start = std::time::Instant::now();

    let mut zeta = FIArrayI64::unit(N);
    let lim = iroot::<8>(N) + 1;
    let primes = wheel_sieve(zeta.isqrt as u64);
    let keys = FIArray::keys(N).collect_vec();

    for &p in &primes {
        let p = p as usize;
        if p >= lim {
            break;
        }
        for (i, &v) in keys.iter().enumerate().rev() {
            if v < p {
                break;
            }
            zeta.arr[i] -= zeta[v / p];
        }
    }
    let zeta_lim = zeta;

    let mut mu0 = FIArrayI64::new(N);
    mu0.arr[0] = 1;
    for i in 1..=mu0.isqrt {
        for &p in &primes {
            let p = p as usize;
            if i * p > mu0.isqrt {
                break;
            }
            if p < lim {
                mu0.arr[i * p - 1] = 0;
                if i % p == 0 {
                    break;
                }
            } else {
                if i % p != 0 {
                    mu0.arr[i * p - 1] = -mu0.arr[i - 1];
                } else {
                    mu0.arr[i * p - 1] = 0;
                    break;
                }
            }
        }
    }
    mu0.partial_sum();

    let mut mu = mult_i64(&zeta_lim, &mu0);
    for e in &mut mu.arr {
        *e -= 1;
    }
    let tmp = mult_sparse_i64(&mu0, &mu);
    for i in 0..mu.arr.len() {
        mu.arr[i] = mu0.arr[i] - tmp.arr[i];
    }
    let i = primes.partition_point(|&p| p < lim as u64);
    for &p in primes[..i].iter().rev() {
        let p = p as usize;
        for (i, &v) in keys.iter().enumerate().rev() {
            if v < p {
                break;
            }
            mu.arr[i] -= mu[v / p];
        }
    }
    let end = start.elapsed();
    println!("Sparse inv (naive + newton iter): {end:?}, {}", mu[N]);

    println!("hello and goodbye");
    println!("summing divisor counts");
    dbg!(fast_divisor_sums::divisor_summatory(N));
    /* {
        let start = std::time::Instant::now();
        let mut zeta = DirichletFenwickI64::zeta(N);
        let lim = iroot::<9>(N) + 1;
        let mut primes = vec![];
        println!("started removal of primes < {lim}: {:?}", start.elapsed());
        for p in 2..lim {
            if zeta.get_bucket_prefix(p - 1) == 1 {
                continue;
            }
            primes.push(p);
            zeta.sparse_mul_at_most_one(p, 1);
        }
        let zeta_lim = FIArrayI64::from(zeta);
        println!(
            "Finished removal of primes < {lim}, started convolution: {:?}",
            start.elapsed()
        );
        let mut zeta_2 = DirichletFenwickI64::from(mult_i64(&zeta_lim, &zeta_lim));
        println!(
            "Finished convolution, started adding back primes < {lim}: {:?}",
            start.elapsed()
        );
        for &p in primes.iter().rev() {
            zeta_2.sparse_mul_unlimited(p, 2);
        }
        let approx = FIArrayI64::from(zeta_2);
        println!(
            "Finished adding back primes < {lim}, started correction: {:?}",
            start.elapsed()
        );
        let accurate = mult_correction(&approx, &primes, |_, _, e| e as i64 + 1);
        let end = start.elapsed();
        dbg!(end, accurate[N]);
    }
     */

    let start = std::time::Instant::now();
    let mut zeta = DirichletFenwickI64::zeta(N);
    let lim = iroot::<8>(N) + 1;
    let mut primes = vec![];
    println!("started removal of primes < {lim}: {:?}", start.elapsed());
    for p in 2..lim {
        if zeta.get_bucket_prefix(p - 1) == 1 {
            continue;
        }
        primes.push(p);
        zeta.sparse_mul_at_most_one(p, 1);
    }
    let zeta_lim = FIArrayI64::from(zeta);
    println!(
        "Finished removal of primes < {lim}, started convolution: {:?}",
        start.elapsed()
    );
    let mut zeta_2 = DirichletFenwickI64::from(mult_i64(&zeta_lim, &zeta_lim));
    println!(
        "Finished convolution, started adding back primes < {lim}: {:?}",
        start.elapsed()
    );
    for &p in primes.iter().rev() {
        zeta_2.sparse_mul_unlimited(p, 2);
    }
    let approx = FIArrayI64::from(zeta_2);
    println!(
        "Finished adding back primes < {lim}, started correction: {:?}",
        start.elapsed()
    );
    let accurate = mult_correction_i64(&approx, &primes, |_, _, e| e as i64 + 1);
    let end = start.elapsed();
    dbg!(end, accurate[N]); // 1e16: 185.5163734s, 1e15: 43.4229948s, 1e14: 9.991973s

    let start = std::time::Instant::now();
    let mut zeta = FIArrayI64::unit(N);
    let lim = iroot::<8>(N) + 1;
    let mut primes = vec![];
    println!("started removal of primes < {lim}: {:?}", start.elapsed());
    let keys = FIArray::keys(N).collect_vec();

    for p in 2..lim {
        if zeta.arr[p - 1] == 1 {
            continue;
        }
        primes.push(p);
        for (i, &v) in keys.iter().enumerate().rev() {
            if v < p {
                break;
            }
            zeta.arr[i] -= zeta[v / p];
        }
    }

    let zeta_lim = zeta;
    println!(
        "Finished removal of primes < {lim}, started convolution: {:?}",
        start.elapsed()
    );
    let mut zeta_2 = mult_i64(&zeta_lim, &zeta_lim);
    println!(
        "Finished convolution, started adding back primes < {lim}: {:?}",
        start.elapsed()
    );
    for &p in primes.iter().rev() {
        for (i, &v) in keys.iter().enumerate().skip(p - 1) {
            zeta_2.arr[i] += 2 * zeta_2[v / p];
        }
    }
    let approx = zeta_2;
    println!(
        "Finished adding back primes < {lim}, started correction: {:?}",
        start.elapsed()
    );
    let accurate = mult_correction_i64(&approx, &primes, |_, _, e| e as i64 + 1);
    let end = start.elapsed();
    dbg!(end, accurate[N]);

    /* {
           let start = std::time::Instant::now();
           let mut zeta = DirichletFenwickI64::zeta(N);
           let lim = iroot::<7>(N) + 1;
           let mut primes = vec![];
           println!("started removal of primes < {lim}: {:?}", start.elapsed());
           for p in 2..lim {
               if zeta.get_bucket_prefix(p - 1) == 1 {
                   continue;
               }
               primes.push(p);
               zeta.sparse_mul_at_most_one(p, 1);
           }
           let zeta_lim = FIArrayI64::from(zeta);
           println!(
               "Finished removal of primes < {lim}, started convolution: {:?}",
               start.elapsed()
           );
           let mut zeta_2 = DirichletFenwickI64::from(mult_i64(&zeta_lim, &zeta_lim));
           println!(
               "Finished convolution, started adding back primes < {lim}: {:?}",
               start.elapsed()
           );
           for &p in primes.iter().rev() {
               zeta_2.sparse_mul_unlimited(p, 2);
           }
           let approx = FIArrayI64::from(zeta_2);
           println!(
               "Finished adding back primes < {lim}, started correction: {:?}",
               start.elapsed()
           );
           let accurate = mult_correction(&approx, &primes, |_, _, e| e as i64 + 1);
           let end = start.elapsed();
           dbg!(end, accurate[N]);
       }
       {
           let start = std::time::Instant::now();
           let mut zeta = DirichletFenwickI64::zeta(N);
           let lim = iroot::<6>(N) + 1;
           let mut primes = vec![];
           println!("started removal of primes < {lim}: {:?}", start.elapsed());
           for p in 2..lim {
               if zeta.get_bucket_prefix(p - 1) == 1 {
                   continue;
               }
               primes.push(p);
               zeta.sparse_mul_at_most_one(p, 1);
           }
           let zeta_lim = FIArrayI64::from(zeta);
           println!(
               "Finished removal of primes < {lim}, started convolution: {:?}",
               start.elapsed()
           );
           let mut zeta_2 = DirichletFenwickI64::from(mult_i64(&zeta_lim, &zeta_lim));
           println!(
               "Finished convolution, started adding back primes < {lim}: {:?}",
               start.elapsed()
           );
           for &p in primes.iter().rev() {
               zeta_2.sparse_mul_unlimited(p, 2);
           }
           let approx = FIArrayI64::from(zeta_2);
           println!(
               "Finished adding back primes < {lim}, started correction: {:?}",
               start.elapsed()
           );
           let accurate = mult_correction(&approx, &primes, |_, _, e| e as i64 + 1);
           let end = start.elapsed();
           dbg!(end, accurate[N]);
       }
       {
           let start = std::time::Instant::now();
           let mut zeta = DirichletFenwickI64::zeta(N);
           let lim = iroot::<5>(N) + 1;
           let mut primes = vec![];
           println!("started removal of primes < {lim}: {:?}", start.elapsed());
           for p in 2..lim {
               if zeta.get_bucket_prefix(p - 1) == 1 {
                   continue;
               }
               primes.push(p);
               zeta.sparse_mul_at_most_one(p, 1);
           }
           let zeta_lim = FIArrayI64::from(zeta);
           println!(
               "Finished removal of primes < {lim}, started convolution: {:?}",
               start.elapsed()
           );
           let mut zeta_2 = DirichletFenwickI64::from(mult_i64(&zeta_lim, &zeta_lim));
           println!(
               "Finished convolution, started adding back primes < {lim}: {:?}",
               start.elapsed()
           );
           for &p in primes.iter().rev() {
               zeta_2.sparse_mul_unlimited(p, 2);
           }
           let approx = FIArrayI64::from(zeta_2);
           println!(
               "Finished adding back primes < {lim}, started correction: {:?}",
               start.elapsed()
           );
           let accurate = mult_correction(&approx, &primes, |_, _, e| e as i64 + 1);
           let end = start.elapsed();
           dbg!(end, accurate[N]);
       }
       {
           let start = std::time::Instant::now();
           let mut zeta = DirichletFenwickI64::zeta(N);
           let lim = iroot::<4>(N) + 1;
           let mut primes = vec![];
           println!("started removal of primes < {lim}: {:?}", start.elapsed());
           for p in 2..lim {
               if zeta.get_bucket_prefix(p - 1) == 1 {
                   continue;
               }
               primes.push(p);
               zeta.sparse_mul_at_most_one(p, 1);
           }
           let zeta_lim = FIArrayI64::from(zeta);
           println!(
               "Finished removal of primes < {lim}, started convolution: {:?}",
               start.elapsed()
           );
           let mut zeta_2 = DirichletFenwickI64::from(mult_i64(&zeta_lim, &zeta_lim));
           println!(
               "Finished convolution, started adding back primes < {lim}: {:?}",
               start.elapsed()
           );
           for &p in primes.iter().rev() {
               zeta_2.sparse_mul_unlimited(p, 2);
           }
           let approx = FIArrayI64::from(zeta_2);
           println!(
               "Finished adding back primes < {lim}, started correction: {:?}",
               start.elapsed()
           );
           let accurate = mult_correction(&approx, &primes, |_, _, e| e as i64 + 1);
           let end = start.elapsed();
           dbg!(end, accurate[N]);
       }
    */
    let start = std::time::Instant::now();
    let mut pi = inverse_pseudo_euler_transform_fraction_i64(FIArrayI64::unit(N));
    let mut primes = vec![];
    for p in 2..=pi.isqrt {
        if pi.arr[p - 1] != pi.arr[p - 2] {
            primes.push(p);
        }
    }
    for i in 0..pi.arr.len() {
        pi.arr[i] *= 2;
    }

    let s2 = mult_correction_i64(
        &pseudo_euler_transform_fraction_i64(pi),
        &primes,
        |_, _, e| e as i64 + 1,
    );
    let end = start.elapsed();
    dbg!(end, s2[N]); // 1e16: 412.19307s, 1e15: 94.1795311s, 1e14: 23.9598301s
    assert_eq!(accurate, s2);

    let start = std::time::Instant::now();
    let u = FIArray::unit(N);
    let s1 = mult(&u, &u);
    let end = start.elapsed();
    dbg!(end, s1[N]); // 1e16: 930.8709046s, 1e15: 163.2886445s, 1e14: 31.1254918s

    for i in 0..s1.arr.len() {
        assert_eq!(s1.arr[i], accurate.arr[i] as usize);
    }
    let start = std::time::Instant::now();
    let u = FIArray::unit(N);
    let s1 = icy_mult(&u, &u); // 1e11 168.5048ms, 1e12 872.4128ms, 1e15 149.5834006s
    let end = start.elapsed();
    dbg!(end, s1[N]);
    println!("hello and goodbye");
    //utils::primes::prime_sieves::main();
    //utils::primes::primecount::main();
    println!("Finished running at: {} ", Local::now().time());
}
pub fn divisor_sums(n: usize) -> FIArray {
    let mut ret = FIArray::new(n);
    let sqrt = ret.isqrt;
    //let start = std::time::Instant::now();
    // Guard tiny n to avoid ln(0/1), step_by(0), sqrt-1 underflow, etc.
    if n <= 1 {
        if n == 1 {
            ret.arr[0] = 1; // sum_{m<=1} tau(m) = 1
        }
        return ret;
    }

    let mut c = ((n as f64).powf(3.0 / 5.0)) as usize;
    //dbg!(c, n, sqrt);
    // Keep c in a sane range; the algorithm expects c >= sqrt+1 to do any block work.
    if c <= sqrt {
        c = sqrt + 1;
    }
    if c > n {
        c = n;
    }

    let mut keys = FIArray::keys(n).enumerate().skip(sqrt);

    // exact tau(m) for m <= sqrt
    ret.arr[..sqrt].fill(1);
    for i in 2..=sqrt {
        for m in (i..=sqrt).step_by(i) {
            ret.arr[m - 1] += 1;
        }
    }
    ret.partial_sum();
    let mut acc = ret.arr[sqrt - 1];
    //dbg!(start.elapsed());
    // ---- FIXED BLOCK SIEVE PART ----
    let mut block_len = c.isqrt();
    if block_len == 0 {
        block_len = 1;
    }
    //dbg!(block_len);
    let mut block = vec![0usize; block_len];

    let mut iv = keys.next();

    // We compute tau(m) for m in [sqrt+1, c) (c is excluded, matching your original range).
    // Important: last block may be shorter than block_len.
    for b in (sqrt + 1..c).step_by(block_len) {
        let len = (c - b).min(block_len);
        block[..len].fill(0);

        let end = b + len - 1;
        let dmax = end.isqrt(); // only need d up to sqrt(end)

        for d in 1..=dmax {
            let Some(dd) = d.checked_mul(d) else { break };

            // First multiple of d in [b, end]
            let mut m = b.next_multiple_of(d);

            // Enforce m >= d*d so we only count each divisor-pair once.
            if m < dd {
                m = dd;
            }
            if m > end {
                continue;
            }

            for mm in (m..=end).step_by(d) {
                // count d and mm/d (two divisors) unless it's a square
                block[mm - b] += 1 + usize::from(dd != mm);
            }
        }

        // turn tau(m) into prefix sums and write out requested keys
        for i in 0..len {
            acc += block[i];
            block[i] = acc;

            // your original “match exact v” logic, but safe for multiple keys at same v
            while let Some((ind, v)) = iv
                && v == b + i
            {
                ret.arr[ind] = acc;
                iv = keys.next();
            }
        }
    }
    // ---- END FIXED BLOCK SIEVE PART ----
    //dbg!(start.elapsed());

    // Remaining keys (v >= c) via divisor_summatory
    if let Some((i, v)) = iv {
        ret.arr[i] = divisor_summatory(v);
        for (i, v) in keys {
            ret.arr[i] = divisor_summatory(v);
        }
    }
    //dbg!(start.elapsed());

    ret
}

fn log_zeta_fast(n: usize) -> FIArray {
    const INVS: [usize; 4] = [0, 6, 3, 2];
    let start = std::time::Instant::now();
    let mut zeta = DirichletFenwick::zeta(n);
    let mut zeta_2 = DirichletFenwick::from(divisor_sums(n));
    dbg!(start.elapsed());
    let rt = zeta.isqrt;
    let len = zeta.bit.0.len();

    let mut ret = FIArray::new(n);

    let x = iroot::<4>(n) + 1;
    // remove contributions of small primes
    for p in 2..x {
        if zeta.bit.sum(p - 1) == 1 {
            //not prime
            continue;
        }
        ret.arr[p - 1] = 1;
        zeta.sparse_mul_at_most_one(p, 1);
        zeta_2.sparse_mul_at_most_one(p, 1);
        zeta_2.sparse_mul_at_most_one(p, 1);
    }
    zeta.bit.dec(0);
    let mut zeta = FIArray::from(zeta);
    let mut zeta_2 = FIArray::from(zeta_2);
    for i in 0..len {
        zeta_2.arr[i] -= 2 * zeta.arr[i] + 1;
    }
    dbg!(start.elapsed());
    // zeta now equals zeta_t - 1, and zeta_2 (zeta_t - 1)^2
    // compute log(zeta_t) using log(x + 1) = x^3 / 3 - x^2 / 2 + x
    // in order to not have to deal with rational numbers, we compute 6 * log(zeta_t)
    // and adjust later

    for i in x..=len {
        ret.arr[i - 1] = zeta.arr[i - 1] * INVS[1];
    }

    let pa = {
        let mut vec = vec![];
        for i in x..=rt {
            if zeta.arr[i - 1] != zeta.arr[i - 2] {
                vec.push(i);
            }
        }
        vec.push(rt + 1);
        vec
    };
    let va = &pa[..pa.len() - 1];

    /* let ind = zeta_2.get_index(x.pow(2));
    assert!(pow_zeta.arr[..ind].iter().all(|&e| e == 0));
    dbg!(ind, x.pow(2), rt - 1, len, len - ind); */

    for i in rt + 1..=len {
        ret.arr[i - 1] -= zeta_2.arr[i - 1] * INVS[2];
    }

    //pow_zeta = mult_sparse(&zeta, &pow_zeta);

    zeta.arr[rt..].fill(0);
    for &i in va {
        for j in 1..=rt / i {
            //zeta[n / j] += pow_zeta[n / (i * j)];
            zeta.arr[len - j] += zeta_2.arr[len - i * j];
        }
    }
    //zeta.arr[..rt].fill(0);
    let zeta_3 = zeta;

    let ind = zeta_3.get_index(x.pow(3));
    //dbg!(ind, len, len - ind);
    //assert!(pow_zeta.arr[..ind].iter().all(|&e| e == 0)); */
    for i in ind + 1..=len {
        ret.arr[i - 1] += zeta_3.arr[i - 1] * INVS[3];
    }
    dbg!(start.elapsed());
    // correction phase: get rid of contributions of prime powers
    for i in (x..len).rev() {
        ret.arr[i] -= ret.arr[i - 1];
    }

    for x in x..=rt {
        let v = ret.arr[x - 1] / 6;
        if v == 0 {
            continue;
        }
        let mut e = 1;
        let mut px = x;
        while px <= n / x {
            e += 1;
            px *= x;

            ret[px] -= INVS[e];
        }
    }
    for i in 1..x - 1 {
        ret.arr[i] += ret.arr[i - 1];
    }
    for i in x - 1..len {
        ret.arr[i] /= 6;
        ret.arr[i] += ret.arr[i - 1];
    }
    ret
}

fn log_zeta_fast_alt(n: usize) -> FIArray {
    fn mult_correction(d: &FIArray, primes: &[usize]) -> FIArray {
        struct Correction(FIArray, usize);
        impl Correction {
            fn fill(&mut self, primes: &[usize], lim: usize, x: usize, y: usize) {
                self.0[x] += y;
                self.1 += 1;
                for (i, &p) in primes.iter().enumerate() {
                    if p > lim / p {
                        break;
                    }
                    let mut pp = p * p;
                    let mut new_lim = lim / pp;
                    for e in 2.. {
                        let hp = 1 << (e - 2);
                        if hp != 0 {
                            self.fill(&primes[i + 1..], new_lim, x * pp, y * hp);
                        }
                        if p > new_lim {
                            break;
                        }
                        pp *= p;
                        new_lim /= p;
                    }
                }
            }
        }
        let mut correction = Correction(FIArray::new(d.x), 0);
        correction.fill(primes, d.x, 1, 1);
        for i in 1..correction.0.arr.len() {
            correction.0.arr[i] += correction.0.arr[i - 1];
        }
        dbg!(correction.1);
        mult_sparse(d, &correction.0)
    }
    const INVS: [usize; 4] = [0, 6, 3, 2];
    let start = std::time::Instant::now();
    let mut zeta = DirichletFenwick::zeta(n);
    let mut zeta_2 = DirichletFenwick::from(divisor_sums(n));
    dbg!(start.elapsed());
    let rt = zeta.isqrt;
    let len = zeta.bit.0.len();

    let mut ret = FIArray::new(n);

    let mut primes = vec![];
    let x = iroot::<4>(n) + 1;
    // remove contributions of small primes
    for p in 2..x {
        if zeta.bit.sum(p - 1) == 1 {
            //not prime
            continue;
        }
        primes.push(p);
        ret.arr[p - 1] = 1;
        zeta.sparse_mul_at_most_one(p, 1);
        zeta_2.sparse_mul_at_most_one(p, 2);
    }
    zeta.bit.dec(0);
    let mut zeta = FIArray::from(zeta);
    let mut zeta_2 = FIArray::from(zeta_2);
    dbg!(start.elapsed());

    zeta_2 = mult_correction(&zeta_2, &primes);
    for i in 0..len {
        zeta_2.arr[i] -= 2 * zeta.arr[i] + 1;
    }
    dbg!(start.elapsed());
    // zeta now equals zeta_t - 1, and zeta_2 (zeta_t - 1)^2
    // compute log(zeta_t) using log(x + 1) = x^3 / 3 - x^2 / 2 + x
    // in order to not have to deal with rational numbers, we compute 6 * log(zeta_t)
    // and adjust later

    for i in x..=len {
        ret.arr[i - 1] = zeta.arr[i - 1] * INVS[1];
    }

    let pa = {
        let mut vec = vec![];
        for i in x..=rt {
            if zeta.arr[i - 1] != zeta.arr[i - 2] {
                vec.push(i);
            }
        }
        vec.push(rt + 1);
        vec
    };
    let va = &pa[..pa.len() - 1];

    /* let ind = zeta_2.get_index(x.pow(2));
    assert!(pow_zeta.arr[..ind].iter().all(|&e| e == 0));
    dbg!(ind, x.pow(2), rt - 1, len, len - ind); */

    for i in rt + 1..=len {
        ret.arr[i - 1] -= zeta_2.arr[i - 1] * INVS[2];
    }

    //pow_zeta = mult_sparse(&zeta, &pow_zeta);

    zeta.arr[rt..].fill(0);
    for &i in va {
        for j in 1..=rt / i {
            //zeta[n / j] += pow_zeta[n / (i * j)];
            zeta.arr[len - j] += zeta_2.arr[len - i * j];
        }
    }
    //zeta.arr[..rt].fill(0);
    let zeta_3 = zeta;

    let ind = zeta_3.get_index(x.pow(3));
    //dbg!(ind, len, len - ind);
    //assert!(pow_zeta.arr[..ind].iter().all(|&e| e == 0)); */
    for i in ind + 1..=len {
        ret.arr[i - 1] += zeta_3.arr[i - 1] * INVS[3];
    }
    dbg!(start.elapsed());
    // correction phase: get rid of contributions of prime powers
    for i in (x..len).rev() {
        ret.arr[i] -= ret.arr[i - 1];
    }

    for x in x..=rt {
        let v = ret.arr[x - 1] / 6;
        if v == 0 {
            continue;
        }
        let mut e = 1;
        let mut px = x;
        while px <= n / x {
            e += 1;
            px *= x;

            ret[px] -= INVS[e];
        }
    }
    for i in 1..x - 1 {
        ret.arr[i] += ret.arr[i - 1];
    }
    for i in x - 1..len {
        ret.arr[i] /= 6;
        ret.arr[i] += ret.arr[i - 1];
    }
    ret
}

fn icy_mult(F: &FIArray, G: &FIArray) -> FIArray {
    unsafe { core::hint::assert_unchecked(F.x == G.x) };
    unsafe { core::hint::assert_unchecked(F.isqrt == G.isqrt) };
    unsafe { core::hint::assert_unchecked(F.arr.len() == G.arr.len()) };

    let f = |i: usize| {
        if i == 1 {
            F.arr[0]
        } else {
            F.arr[i - 1] - F.arr[i - 2]
        }
    };
    let g = |i: usize| {
        if i == 1 {
            G.arr[0]
        } else {
            G.arr[i - 1] - G.arr[i - 2]
        }
    };
    let N = F.x;
    let c_N = iroot::<3>(N);

    let mut H = FIArray::new(N);
    let len = H.arr.len();
    for x in 1..=c_N {
        let (fx, gx) = (f(x), g(x));
        let max_y = (N / x).isqrt();
        H[x * x] += fx * gx;
        for y in x + 1..=max_y {
            H[x * y] += fx * g(y) + gx * f(y);
        }
        H.arr[len - max_y] -= fx * G.arr[max_y - 1] + gx * F.arr[max_y - 1];
    }
    H.arr[len - c_N] += F.arr[c_N - 1] * G.arr[c_N - 1];
    H.partial_sum();
    for x in 1..=c_N {
        let Nx = N / x;
        let (fx, gx) = (f(x), g(x));
        let max_y = Nx.isqrt();
        H.arr[len - x] += fx * G[Nx / x] + gx * F[Nx / x];
        for y in x + 1..=max_y {
            H.arr[len - x] += f(y) * G[Nx / y] + g(y) * F[Nx / y];
            H.arr[len - y] += fx * G[Nx / y] + gx * F[Nx / y];
        }
        H.arr[len - x] -= G.arr[max_y - 1] * F.arr[max_y - 1];
    }
    H
}
#[must_use]
pub fn pseudo_euler_transform_i64(a: FIArrayI64) -> FIArrayI64 {
    let n = a.x;
    let rt = a.isqrt;
    let x = iroot::<3>(n) + 1;
    let mut r = a;
    let a = r.clone();
    for e in &mut r.arr[x - 1..] {
        *e -= a.arr[x - 2];
    }
    r.arr[..x - 1].fill(0);
    for i in x..=rt {
        let vi = a.arr[i - 1] - a.arr[i - 2];
        if vi == 0 {
            continue;
        }
        let n_i = n / i;
        for z in 1..=n_i / i {
            r[(n / z) as _] += vi * (a[n_i / z] - a.arr[i - 2]);
        }
    }
    let mut r = DirichletFenwickI64::from(r);
    r.bit.add(0, 1);
    for i in (2..x).rev() {
        let cur = a.arr[i - 1] - a.arr[i - 2];
        if cur == 0 {
            continue;
        }
        r.sparse_mul_unlimited(i, cur);
    }
    r.into()
}

#[must_use]
pub fn inverse_pseudo_euler_transform_i64(a: FIArrayI64) -> FIArrayI64 {
    let n = a.x;
    let rt = a.isqrt;
    let len = a.arr.len();
    let mut r = FIArrayI64::new(n);
    let mut a_bit = DirichletFenwickI64::from(a);
    let x = iroot::<3>(n) + 1;
    for i in 2..x {
        let cur = a_bit.bit.sum(i - 1) - 1;
        if cur == 0 {
            continue;
        }
        r.arr[i - 1] = cur;
        a_bit.sparse_mul_at_most_one(i, cur);
    }
    let a = FIArrayI64::from(a_bit);
    for i in x..=len {
        r.arr[i - 1] = a.arr[i - 1] - a.arr[i - 2];
    }
    for i in x..=rt {
        let vi = r.arr[i - 1];
        if vi == 0 {
            continue;
        }
        let n_i = n / i;
        for z in 1..=n_i / i {
            let v = vi * (a[n_i / z] - a.arr[i - 2]);
            r.arr[len - z] -= v;
            if z > 1 {
                r.arr[len - z + 1] += v;
            }
        }
    }
    for i in 1..len {
        r.arr[i] += r.arr[i - 1];
    }
    r
}

#[must_use]
pub fn mult_correction_single_i64(
    d: &FIArrayI64,
    primes: &[usize],
    f: impl Fn(usize, usize, u8) -> i64,
) -> i64 {
    fn dfs(
        d: &FIArrayI64,
        primes: &[usize],
        lim: usize,
        x: usize,
        hx: i64,
        f: &impl Fn(usize, usize, u8) -> i64,
    ) -> i64 {
        let mut res = hx * d[d.x / x];
        for (i, &p) in primes.iter().enumerate() {
            if p > lim / p {
                break;
            }
            let fp = f(p, p, 1);
            let mut prev = fp;
            let mut pp = p * p;
            let mut new_lim = lim / pp;
            for e in 2.. {
                let cur = f(pp, p, e);
                let hp = cur - fp * prev;
                if hp != 0 {
                    res += dfs(d, &primes[i + 1..], new_lim, x * pp, hx * hp, f);
                }
                prev = cur;
                if p > new_lim {
                    break;
                }
                pp *= p;
                new_lim /= p;
            }
        }
        res
    }

    dfs(d, primes, d.x, 1, 1, &f)
}

#[must_use]
pub fn mult_correction(
    d: &FIArray,
    primes: &[usize],
    f: impl Fn(usize, usize, u8) -> usize,
) -> FIArray {
    struct Correction(FIArray, usize);
    impl Correction {
        fn fill(
            &mut self,
            primes: &[usize],
            lim: usize,
            x: usize,
            y: usize,
            f: &impl Fn(usize, usize, u8) -> usize,
        ) {
            self.0[x] += y;
            self.1 += 1;
            for (i, &p) in primes.iter().enumerate() {
                if p > lim / p {
                    break;
                }
                let fp = f(p, p, 1);
                let mut prev = fp;
                let mut pp = p * p;
                let mut new_lim = lim / pp;
                for e in 2.. {
                    let cur = f(pp, p, e);
                    let hp = cur - fp * prev;
                    dbg!((p, e, hp));
                    if hp != 0 {
                        self.fill(&primes[i + 1..], new_lim, x * pp, y * hp, f);
                    }
                    prev = cur;
                    if p > new_lim {
                        break;
                    }
                    pp *= p;
                    new_lim /= p;
                }
            }
        }
    }
    let mut correction = Correction(FIArray::new(d.x), 0);
    correction.fill(primes, d.x, 1, 1, &f);
    for i in 1..correction.0.arr.len() {
        correction.0.arr[i] += correction.0.arr[i - 1];
    }
    dbg!(correction.1);
    mult_sparse(d, &correction.0)
}

#[must_use]
pub fn mult_correction_i64(
    d: &FIArrayI64,
    primes: &[usize],
    f: impl Fn(usize, usize, u8) -> i64,
) -> FIArrayI64 {
    struct Correction(FIArrayI64, usize);
    impl Correction {
        fn fill(
            &mut self,
            primes: &[usize],
            lim: usize,
            x: usize,
            y: i64,
            f: &impl Fn(usize, usize, u8) -> i64,
        ) {
            self.0[x as _] += y;
            self.1 += 1;
            for (i, &p) in primes.iter().enumerate() {
                if p > lim / p {
                    break;
                }
                let fp = f(p, p, 1);
                let mut prev = fp;
                let mut pp = p * p;
                let mut new_lim = lim / pp;
                for e in 2.. {
                    let cur = f(pp, p, e);
                    let hp = cur - fp * prev;
                    if hp != 0 {
                        self.fill(&primes[i + 1..], new_lim, x * pp, y * hp, f);
                    }
                    prev = cur;
                    if p > new_lim {
                        break;
                    }
                    pp *= p;
                    new_lim /= p;
                }
            }
        }
    }
    let mut correction = Correction(FIArrayI64::new(d.x), 0);
    correction.fill(primes, d.x, 1, 1, &f);
    for i in 1..correction.0.arr.len() {
        correction.0.arr[i] += correction.0.arr[i - 1];
    }
    dbg!(correction.1);
    mult_sparse_i64(d, &correction.0)
}
pub fn mult_sparse_with_buffer_i64(a: &FIArrayI64, b: &FIArrayI64, res: &mut FIArrayI64) {
    unsafe { core::hint::assert_unchecked(a.x == b.x && a.x == res.x) };
    res.arr.fill(0);
    let R2 = a.isqrt;
    let n = a.x;
    let s1 = |ds: &FIArrayI64| {
        let mut vec = vec![];
        if ds.arr[0] != 0 {
            vec.push((1, ds.arr[0]));
        }
        for i in 1..R2 {
            if ds.arr[i] != ds.arr[i - 1] {
                vec.push((i + 1, ds.arr[i] - ds.arr[i - 1]));
            }
        }
        vec.push((R2 + 1, 0));
        vec
    };
    let pa = s1(a);
    let pb = s1(b);
    let va = &pa[..pa.len() - 1];
    let vb = &pb[..pb.len() - 1];
    let len = res.arr.len();

    let mut r = va.len();
    for &(x, y) in vb {
        let Nx = n / x;
        while r != 0 && Nx / pa[r - 1].0 < r - 1 {
            r -= 1;
        }
        let mut i = 0;
        let X = R2 / x;
        let Nx = n / x;
        while i != r && pa[i].0 <= X {
            res.arr[x * pa[i].0 - 1] += y * pa[i].1;
            i += 1;
        }
        while i != r {
            res.arr[len - Nx / pa[i].0] += y * pa[i].1;
            i += 1;
        }
        if r != 0 && pa[r].0 <= Nx {
            res[(x * pa[r].0) as _] -= y * a.arr[pa[r - 1].0 - 1];
        }
    }
    for &(x, y) in va {
        res.arr[len - R2 / x] -= y * b.arr[R2 - 1];
    }
    for i in 1..len {
        res.arr[i] += res.arr[i - 1];
    }
    let mut r = va.len();
    for &(x, y) in vb {
        let Nx = n / x;
        while r != 0 && Nx / pa[r - 1].0 < r - 1 {
            r -= 1;
        }
        let mut i = Nx / pa[r].0;
        let X = R2 / x;
        while i > X {
            res.arr[len - i] += y * a.arr[Nx / i - 1];
            i -= 1;
        }
        while i > 0 {
            res.arr[len - i] += y * a.arr[len - x * i];
            i -= 1;
        }
    }
    for &(x, y) in va {
        for j in 1..=R2 / x {
            res.arr[len - j] += y * b.arr[len - x * j];
        }
    }
}
#[must_use]
pub fn mult_sparse_i64(a: &FIArrayI64, b: &FIArrayI64) -> FIArrayI64 {
    let mut res = a.clone();
    mult_sparse_with_buffer_i64(a, b, &mut res);
    res
}
pub fn mult_with_buffer_i64(a: &FIArrayI64, b: &FIArrayI64, res: &mut FIArrayI64) {
    res.arr.fill(0);
    unsafe { core::hint::assert_unchecked(a.x == b.x) };
    let R2 = a.isqrt;
    let n = a.x;
    let mut res = FIArrayI64::new(n);

    let s1 = |ds: &FIArrayI64| {
        let mut vec = vec![];
        if ds.arr[0] != 0 {
            vec.push((1, ds.arr[0]));
        }
        for i in 1..R2 {
            if ds.arr[i] != ds.arr[i - 1] {
                vec.push((i + 1, ds.arr[i] - ds.arr[i - 1]));
            }
        }
        vec.push((R2 + 1, 0));
        vec
    };
    let len = res.arr.len();
    if a == b {
        let pa = s1(a);
        let va = &pa[..pa.len() - 1];
        let mut r = va.len();
        let mut l = 0;
        for &(x, fx) in va {
            res[x * x] += fx * fx;

            l += 1;
            let Nx = n / x;
            while r > l && Nx / pa[r - 1].0 < r - 1 {
                r -= 1;
            }
            if r < l {
                r = l;
            }
            let X = R2 / x;
            let Nx = n / x;

            let mut i = l;

            while i != r && pa[i].0 <= X {
                res.arr[x * pa[i].0 - 1] += 2 * fx * pa[i].1;
                i += 1;
            }
            while i != r {
                res.arr[len - Nx / pa[i].0] += 2 * fx * pa[i].1;
                i += 1;
            }
            if r != 0 && pa[r].0 <= Nx {
                res[x * pa[r].0] -= fx * a.arr[pa[r - 1].0 - 1] * 2;
            }
        }
        for i in 1..len {
            res.arr[i] += res.arr[i - 1];
        }
        r = va.len();
        l = 0;
        for &(x, fx) in va {
            l += 1;
            let Nx = n / x;
            while r > l && Nx / pa[r - 1].0 < r - 1 {
                r -= 1;
            }
            if r < l {
                r = l;
            }

            let mut i = Nx / pa[r].0;
            let X = R2 / x;
            while i > X {
                res.arr[len - i] += 2 * fx * a.arr[Nx / i - 1];
                i -= 1;
            }
            while i > 0 {
                res.arr[len - i] += 2 * fx * a.arr[len - x * i];
                i -= 1;
            }
        }
    } else {
        let pa = s1(a);
        let va = &pa[..pa.len() - 1];
        let pb = s1(b);
        let vb = &pb[..pb.len() - 1];
        res.arr[0] += a.arr[0] * b.arr[0];
        for i in 2..=R2 {
            res[i * i] += (a.arr[i - 1] - a.arr[i - 2]) * (b.arr[i - 1] - b.arr[i - 2]);
        }
        let mut r = vb.len();
        let mut l = 0;
        for &(x, y) in va {
            while pb[l].0 <= x {
                l += 1;
            }
            let Nx = n / x;
            while r > l && Nx / pb[r - 1].0 < r - 1 {
                r -= 1;
            }
            if r < l {
                r = l;
            }
            let X = R2 / x;
            let Nx = n / x;

            let mut i = l;

            while i != r && pb[i].0 <= X {
                res.arr[x * pb[i].0 - 1] += y * pb[i].1;
                i += 1;
            }
            while i != r {
                res.arr[len - Nx / pb[i].0] += y * pb[i].1;
                i += 1;
            }
            if r != 0 && pb[r].0 <= Nx {
                res[(x * pb[r].0) as _] -= y * b.arr[pb[r - 1].0 - 1];
            }
        }
        r = va.len();
        l = 0;
        for &(x, y) in vb {
            while pa[l].0 <= x {
                l += 1;
            }
            let Nx = n / x;
            while r > l && Nx / pa[r - 1].0 < r - 1 {
                r -= 1;
            }
            if r < l {
                r = l;
            }
            let X = R2 / x;
            let Nx = n / x;

            let mut i = l;

            while i != r && pa[i].0 <= X {
                res.arr[x * pa[i].0 - 1] += y * pa[i].1;
                i += 1;
            }
            while i != r {
                res.arr[len - Nx / pa[i].0] += y * pa[i].1;
                i += 1;
            }
            if r != 0 && pa[r].0 <= Nx {
                res[(x * pa[r].0) as _] -= y * a.arr[pa[r - 1].0 - 1];
            }
        }
        for i in 1..len {
            res.arr[i] += res.arr[i - 1];
        }
        r = vb.len();
        l = 0;
        for &(x, y) in va {
            while pb[l].0 <= x {
                l += 1;
            }
            let Nx = n / x;
            while r > l && Nx / pb[r - 1].0 < r - 1 {
                r -= 1;
            }
            if r < l {
                r = l;
            }

            let mut i = Nx / pb[r].0;
            let X = R2 / x;
            while i > X {
                res.arr[len - i] += y * b.arr[Nx / i - 1];
                i -= 1;
            }
            while i > 0 {
                res.arr[len - i] += y * b.arr[len - x * i];
                i -= 1;
            }
        }
        r = va.len();
        l = 0;
        for &(x, y) in vb {
            while pa[l].0 <= x {
                l += 1;
            }
            let Nx = n / x;
            while r > l && Nx / pa[r - 1].0 < r - 1 {
                r -= 1;
            }
            if r < l {
                r = l;
            }

            let mut i = Nx / pa[r].0;
            let X = R2 / x;
            while i > X {
                res.arr[len - i] += y * a.arr[Nx / i - 1];
                i -= 1;
            }
            while i > 0 {
                res.arr[len - i] += y * a.arr[len - x * i];
                i -= 1;
            }
        }
    }
}
#[must_use]
pub fn mult_i64(a: &FIArrayI64, b: &FIArrayI64) -> FIArrayI64 {
    unsafe { core::hint::assert_unchecked(a.x == b.x) };
    let R2 = a.isqrt;
    let n = a.x;
    let mut res = FIArrayI64::new(n);

    let s1 = |ds: &FIArrayI64| {
        let mut vec = vec![];
        if ds.arr[0] != 0 {
            vec.push((1, ds.arr[0]));
        }
        for i in 1..R2 {
            if ds.arr[i] != ds.arr[i - 1] {
                vec.push((i + 1, ds.arr[i] - ds.arr[i - 1]));
            }
        }
        vec.push((R2 + 1, 0));
        vec
    };
    let len = res.arr.len();
    if a == b {
        let pa = s1(a);
        let va = &pa[..pa.len() - 1];
        let mut r = va.len();
        let mut l = 0;
        for &(x, fx) in va {
            res[x * x] += fx * fx;

            l += 1;
            let Nx = n / x;
            while r > l && Nx / pa[r - 1].0 < r - 1 {
                r -= 1;
            }
            if r < l {
                r = l;
            }
            let X = R2 / x;
            let Nx = n / x;

            let mut i = l;

            while i != r && pa[i].0 <= X {
                res.arr[x * pa[i].0 - 1] += 2 * fx * pa[i].1;
                i += 1;
            }
            while i != r {
                res.arr[len - Nx / pa[i].0] += 2 * fx * pa[i].1;
                i += 1;
            }
            if r != 0 && pa[r].0 <= Nx {
                res[x * pa[r].0] -= fx * a.arr[pa[r - 1].0 - 1] * 2;
            }
        }
        for i in 1..len {
            res.arr[i] += res.arr[i - 1];
        }
        r = va.len();
        l = 0;
        for &(x, fx) in va {
            l += 1;
            let Nx = n / x;
            while r > l && Nx / pa[r - 1].0 < r - 1 {
                r -= 1;
            }
            if r < l {
                r = l;
            }

            let mut i = Nx / pa[r].0;
            let X = R2 / x;
            while i > X {
                res.arr[len - i] += 2 * fx * a.arr[Nx / i - 1];
                i -= 1;
            }
            while i > 0 {
                res.arr[len - i] += 2 * fx * a.arr[len - x * i];
                i -= 1;
            }
        }
    } else {
        let pa = s1(a);
        let va = &pa[..pa.len() - 1];
        let pb = s1(b);
        let vb = &pb[..pb.len() - 1];
        res.arr[0] += a.arr[0] * b.arr[0];
        for i in 2..=R2 {
            res[i * i] += (a.arr[i - 1] - a.arr[i - 2]) * (b.arr[i - 1] - b.arr[i - 2]);
        }
        let mut r = vb.len();
        let mut l = 0;
        for &(x, y) in va {
            while pb[l].0 <= x {
                l += 1;
            }
            let Nx = n / x;
            while r > l && Nx / pb[r - 1].0 < r - 1 {
                r -= 1;
            }
            if r < l {
                r = l;
            }
            let X = R2 / x;
            let Nx = n / x;

            let mut i = l;

            while i != r && pb[i].0 <= X {
                res.arr[x * pb[i].0 - 1] += y * pb[i].1;
                i += 1;
            }
            while i != r {
                res.arr[len - Nx / pb[i].0] += y * pb[i].1;
                i += 1;
            }
            if r != 0 && pb[r].0 <= Nx {
                res[(x * pb[r].0) as _] -= y * b.arr[pb[r - 1].0 - 1];
            }
        }
        r = va.len();
        l = 0;
        for &(x, y) in vb {
            while pa[l].0 <= x {
                l += 1;
            }
            let Nx = n / x;
            while r > l && Nx / pa[r - 1].0 < r - 1 {
                r -= 1;
            }
            if r < l {
                r = l;
            }
            let X = R2 / x;
            let Nx = n / x;

            let mut i = l;

            while i != r && pa[i].0 <= X {
                res.arr[x * pa[i].0 - 1] += y * pa[i].1;
                i += 1;
            }
            while i != r {
                res.arr[len - Nx / pa[i].0] += y * pa[i].1;
                i += 1;
            }
            if r != 0 && pa[r].0 <= Nx {
                res[(x * pa[r].0) as _] -= y * a.arr[pa[r - 1].0 - 1];
            }
        }
        for i in 1..len {
            res.arr[i] += res.arr[i - 1];
        }
        r = vb.len();
        l = 0;
        for &(x, y) in va {
            while pb[l].0 <= x {
                l += 1;
            }
            let Nx = n / x;
            while r > l && Nx / pb[r - 1].0 < r - 1 {
                r -= 1;
            }
            if r < l {
                r = l;
            }

            let mut i = Nx / pb[r].0;
            let X = R2 / x;
            while i > X {
                res.arr[len - i] += y * b.arr[Nx / i - 1];
                i -= 1;
            }
            while i > 0 {
                res.arr[len - i] += y * b.arr[len - x * i];
                i -= 1;
            }
        }
        r = va.len();
        l = 0;
        for &(x, y) in vb {
            while pa[l].0 <= x {
                l += 1;
            }
            let Nx = n / x;
            while r > l && Nx / pa[r - 1].0 < r - 1 {
                r -= 1;
            }
            if r < l {
                r = l;
            }

            let mut i = Nx / pa[r].0;
            let X = R2 / x;
            while i > X {
                res.arr[len - i] += y * a.arr[Nx / i - 1];
                i -= 1;
            }
            while i > 0 {
                res.arr[len - i] += y * a.arr[len - x * i];
                i -= 1;
            }
        }
    }
    res
}
#[must_use]
pub fn div_i64(a: &FIArrayI64, b: &FIArrayI64) -> FIArrayI64 {
    unsafe { core::hint::assert_unchecked(a.x == b.x) };
    let R2 = a.isqrt;
    let n = a.x;
    let mut res = FIArrayI64::new(n);

    let s1 = |ds: &FIArrayI64| {
        let mut vec = vec![];
        if ds.arr[0] != 0 {
            vec.push((1, ds.arr[0]));
        }
        for i in 1..R2 {
            if ds.arr[i] != ds.arr[i - 1] {
                vec.push((i + 1, ds.arr[i] - ds.arr[i - 1]));
            }
        }
        vec.push((R2 + 1, 0));
        vec
    };
    let len = res.arr.len();

    let mut pa = vec![];
    let pb = s1(b);
    let vb = &pb[..pb.len() - 1];
    res.arr[0] = a.arr[0];
    for i in 1..R2 {
        res.arr[i] = a.arr[i] - a.arr[i - 1];
    }
    let mut sum = 0;
    for i in 1..=R2 {
        let val = res.arr[i - 1];
        sum += val;
        res.arr[i - 1] = sum;
        //dbg!(sum);
        if val == 0 {
            continue;
        }
        pa.push((i, val));
        for (y, fy) in &vb[1..] {
            if y * i > R2 {
                break;
            }
            res.arr[i * y - 1] -= val * fy;
        }
    }
    //dbg!(&res.arr[..R2]);
    pa.push((R2 + 1, 0));
    let va = &pa[..pa.len() - 1];

    let mut r = vb.len();
    let mut l0 = r;
    let mut l = 0;
    for &(x, fx) in va {
        while pb[l].0 <= x {
            l += 1;
        }
        let Nx = n / x;
        let X = R2 / x;

        while l0 > l && X < pb[l0 - 1].0 {
            l0 -= 1;
        }
        while r > l && Nx < (r - 1) * pb[r - 1].0 {
            r -= 1;
        }

        if l.max(l0) < r {
            for (y, fy) in &pb[l.max(l0)..r] {
                res.arr[len - Nx / y] += fx * fy;
            }
        }
        r = r.max(l);

        if r > 0 && pb[r].0 <= Nx {
            res.arr[len - Nx / pb[r].0] -= fx * b.arr[pb[r - 1].0 - 1];
        }
    }
    r = va.len();
    l0 = r;
    l = 0;
    let mut bound_z = n / (R2 + 1);
    for &(y, fy) in vb {
        while pa[l].0 < y {
            l += 1;
        }
        let Ny = n / y;
        let bound_y = Ny / (bound_z + 1);
        while l0 > l && R2 < y * pa[l0 - 1].0 {
            l0 -= 1;
        }
        while r > l && pa[r - 1].0 > bound_y && Ny < (r - 1) * pa[r - 1].0 {
            r -= 1;
        }
        if l.max(l0) < r {
            for (x, fx) in &va[l.max(l0)..r] {
                res.arr[len - Ny / x] += fy * fx;
            }
        }

        r = r.max(l);

        if r > 0 && pa[r].0 <= Ny {
            res.arr[len - Ny / pa[r].0] -= fy * res.arr[pa[r - 1].0 - 1];
        }
        bound_z = Ny / pa[r].0;
    }

    res.arr[R2] += a.arr[R2 - 1];
    for i in R2 + 1..len {
        res.arr[i] += res.arr[i - 1];
    }
    l = 0;
    r = vb.len();
    for &(x, fx) in va {
        while pb[l].0 <= x {
            l += 1;
        }
        let Nx = n / x;
        while r > l && Nx < (r - 1) * pb[r - 1].0 {
            r -= 1;
        }
        //assert!(r >= l);
        //mul_sparse_large(x, Nx, fx, pairs2[std::max(l, r)].first, block2, block_out);
        if r < l {
            r = l;
        }

        let mut i = Nx / pb[r].0;
        let X = R2 / x;
        while i > X {
            res.arr[len - i] += fx * b.arr[Nx / i - 1];
            i -= 1;
        }
        while i > 0 {
            res.arr[len - i] += fx * b.arr[len - x * i];
            i -= 1;
        }
    }
    let mut c_y = 0;
    l = 0;
    r = va.len();
    bound_z = n / (R2 + 1);
    for i in (1..=bound_z).rev() {
        while bound_z >= i {
            c_y += 1;
            let y = pb[c_y].0;
            while pa[l].0 < y {
                l += 1;
            }
            let Ny = n / y;
            let bound_y = Ny / (bound_z + 1);
            while r > l && pa[r - 1].0 > bound_y && (r - 1) * (pa[r - 1].0) > Ny {
                r -= 1;
            }
            r = r.max(l);
            bound_z = Ny / pa[r].0;
        }
        let Nz = n / i;
        let mut ans = a.arr[len - i] - res.arr[len - i];
        for (y, fy) in &pb[1..c_y] {
            ans -= fy * res[Nz / y];
        }
        res.arr[len - i] = ans;
    }
    res
}
// TODO: finish
#[must_use]
pub fn inv_i64(b: &FIArrayI64) -> FIArrayI64 {
    todo!();
    let R2 = b.isqrt;
    let n = b.x;
    let mut res = FIArrayI64::new(n);

    let s1 = |ds: &FIArrayI64| {
        let mut vec = vec![];
        if ds.arr[0] != 0 {
            vec.push((1, ds.arr[0]));
        }
        for i in 1..R2 {
            if ds.arr[i] != ds.arr[i - 1] {
                vec.push((i + 1, ds.arr[i] - ds.arr[i - 1]));
            }
        }
        vec.push((R2 + 1, 0));
        vec
    };
    let len = res.arr.len();

    let pa = [(1, 1), (R2 + 1, 0)];
    let va = &pa[..pa.len() - 1];

    let pb = s1(b);
    let vb = &pb[..pb.len() - 1];
    res.arr[0] = 1;
    let mut sum = 0;
    for i in 1..=R2 {
        let val = res.arr[i - 1];
        sum += val;
        res.arr[i - 1] = sum;
        dbg!(sum);
        if val == 0 {
            continue;
        }
        for (y, fy) in &vb[1..] {
            if y * i > R2 {
                break;
            }
            res.arr[i * y - 1] -= val * fy;
        }
    }
    dbg!(&res.arr[..R2]);
    let mut r = vb.len();
    let mut l0 = r;
    let mut l = 0;
    for &(x, fx) in va {
        while pb[l].0 <= x {
            l += 1;
        }
        let Nx = n / x;
        let X = R2 / x;

        while l0 > l && X < pb[l0 - 1].0 {
            l0 -= 1;
        }
        while r > l && Nx < (r - 1) * pb[r - 1].0 {
            r -= 1;
        }

        if l.max(l0) < r {
            for (y, fy) in &pb[l.max(l0)..r] {
                res.arr[len - Nx / y] += fx * fy;
            }
        }
        r = r.max(l);

        if r > 0 && pb[r].0 <= Nx {
            res.arr[len - Nx / pb[r].0] -= fx * b.arr[pb[r - 1].0 - 1];
        }
    }
    r = va.len();
    l0 = r;
    l = 0;
    let mut bound_z = n / (R2 + 1);
    for &(y, fy) in vb {
        while pa[l].0 < y {
            l += 1;
        }
        let Ny = n / y;
        let bound_y = Ny / (bound_z + 1);
        while l0 > l && R2 < y * pa[l0 - 1].0 {
            l0 -= 1;
        }
        while r > l && pa[r - 1].0 > bound_y && Ny < (r - 1) * pa[r - 1].0 {
            r -= 1;
        }
        if l.max(l0) < r {
            for (x, fx) in &va[l.max(l0)..r] {
                res.arr[len - Ny / x] += fy * fx;
            }
        }

        r = r.max(l);

        if r > 0 && pa[r].0 <= Ny {
            res.arr[len - Ny / pa[r].0] -= fy * res.arr[pa[r - 1].0 - 1];
        }
        bound_z = Ny / pa[r].0;
    }

    res.arr[R2] += 1;
    for i in R2 + 1..len {
        res.arr[i] += res.arr[i - 1];
    }
    l = 0;
    r = vb.len();
    for &(x, fx) in va {
        while pb[l].0 <= x {
            l += 1;
        }
        let Nx = n / x;
        while r > l && Nx < (r - 1) * pb[r - 1].0 {
            r -= 1;
        }
        //assert!(r >= l);
        //mul_sparse_large(x, Nx, fx, pairs2[std::max(l, r)].first, block2, block_out);
        if r < l {
            r = l;
        }

        let mut i = Nx / pb[r].0;
        let X = R2 / x;
        while i > X {
            res.arr[len - i] += fx * b.arr[Nx / i - 1];
            i -= 1;
        }
        while i > 0 {
            res.arr[len - i] += fx * b.arr[len - x * i];
            i -= 1;
        }
    }
    let mut c_y = 0;
    l = 0;
    r = va.len();
    bound_z = n / (R2 + 1);
    for i in (1..=bound_z).rev() {
        while bound_z >= i {
            c_y += 1;
            let y = pb[c_y].0;
            while pa[l].0 < y {
                l += 1;
            }
            let Ny = n / y;
            let bound_y = Ny / (bound_z + 1);
            while r > l && pa[r - 1].0 > bound_y && (r - 1) * (pa[r - 1].0) > Ny {
                r -= 1;
            }
            r = r.max(l);
            bound_z = Ny / pa[r].0;
        }
        let Nz = n / i;
        let mut ans = 1 - res.arr[len - i];
        for (y, fy) in &pb[1..c_y] {
            ans -= fy * res[Nz / y];
        }
        res.arr[len - i] = ans;
    }
    res
}

/* // faster by a log factor, but much more susceptible to overflow - multiplies input by a factor of 1296 before reducing
#[must_use]
pub fn pseudo_euler_transform_fraction_i64(a: FIArrayI64) -> FIArrayI64 {
    const INVS: [i64; 4] = [0, 6, 3, 2];

    let mut a = a;
    let len = a.arr.len();
    let n = a.x;
    let rt = a.isqrt;
    let x = iroot::<4>(n) + 1;

    for i in (1..len).rev() {
        a.arr[i] -= a.arr[i - 1];
        a.arr[i] *= INVS[1];
    }
    a.arr[0] *= INVS[1]; // kinda pointless tbh
    for i in (x..=rt).rev() {
        let v = a.arr[i - 1] / INVS[1];
        if v == 0 {
            continue;
        }
        let mut e = 1;
        let mut pi = i;
        let mut pv = v;
        while pi <= n / i {
            e += 1;
            pi *= i;
            pv *= v;
            a[pi] += pv * INVS[e];
        }
    }

    let mut v = FIArrayI64::new(n);
    for i in x..=len {
        v.arr[i - 1] = v.arr[i - 2] + a.arr[i - 1];
    }

    let mut ret = v.clone();
    for e in &mut ret.arr {
        *e = (*e + INVS[1]) * 6 * INVS[1].pow(2);
    }

    let mut v_2 = mult_i64(&v, &v);
    for i in x..=len {
        ret.arr[i - 1] += v_2.arr[i - 1] * 3 * INVS[1];
    }
    {
        //v_2 = mult_sparse(&v, &v_2);
        v.arr[rt..].fill(0);
        for i in x..=rt {
            let y = v.arr[i - 1] - v.arr[i - 2];
            if y != 0 {
                for j in 1..=rt / i {
                    v.arr[len - j] += y * v_2.arr[len - i * j];
                }
            }
        }
        v.arr[..rt].fill(0);
        core::mem::swap(&mut v_2, &mut v);
    }

    for i in 1..=len {
        ret.arr[i - 1] += v_2.arr[i - 1];
        ret.arr[i - 1] /= const { 6 * INVS[1].pow(3) };
    }
    let mut ret = DirichletFenwickI64::from(ret);
    for i in (2..x).rev() {
        let ai = a.arr[i - 1] / INVS[1];
        if ai == 0 {
            continue;
        }
        ret.sparse_mul_unlimited(i, ai);
    }
    ret.into()
}
 */
// faster by a log factor, but much more susceptible to overflow - multiplies input by a factor of 16 before reducing
#[must_use]
pub fn pseudo_euler_transform_fraction_i64(a: FIArrayI64) -> FIArrayI64 {
    const fn inv_odd(mut k: i64) -> i64 {
        let mut exp = (1u64 << 63) - 1;

        let mut r: i64 = 1;
        while exp > 1 {
            r = r.overflowing_mul(k).0;
            k = k.overflowing_mul(k).0;
            exp >>= 1;
        }
        r.overflowing_mul(k).0
    }

    const INVS: [i64; 4] = [0, 2, 1, inv_odd(3) << 1];

    let mut a = a;
    let len = a.arr.len();
    let n = a.x;
    let rt = a.isqrt;
    let x = iroot::<4>(n) + 1;

    for i in (1..len).rev() {
        a.arr[i] -= a.arr[i - 1];
        a.arr[i] *= INVS[1];
    }
    a.arr[0] *= INVS[1]; // kinda pointless tbh
    for i in (x..=rt).rev() {
        let v = a.arr[i - 1] / INVS[1];
        if v == 0 {
            continue;
        }
        let mut e = 1;
        let mut pi = i;
        let mut pv = v;
        while pi <= n / i {
            e += 1;
            pi *= i;
            pv *= v;
            a[pi] += pv * INVS[e];
        }
    }

    let mut v = FIArrayI64::new(n);
    for i in x..=len {
        v.arr[i - 1] = v.arr[i - 2] + a.arr[i - 1];
    }

    let mut ret = v.clone();
    for e in &mut ret.arr {
        *e = (*e + INVS[1]) * 2 * INVS[1].pow(2);
    }

    let mut v_2 = mult_i64(&v, &v);
    for i in x..=len {
        ret.arr[i - 1] += v_2.arr[i - 1] * INVS[1];
    }
    {
        //v_2 = mult_sparse(&v, &v_2);
        v.arr[rt..].fill(0);
        for i in x..=rt {
            let y = v.arr[i - 1] - v.arr[i - 2];
            if y != 0 {
                for j in 1..=rt / i {
                    v.arr[len - j] += y * v_2.arr[len - i * j];
                }
            }
        }
        v.arr[..rt].fill(0);
        core::mem::swap(&mut v_2, &mut v);
    }

    for i in 1..=len {
        ret.arr[i - 1] += v_2.arr[i - 1] * const { inv_odd(3) };
        ret.arr[i - 1] /= const { 2 * INVS[1].pow(3) };
    }
    /* let mut ret = DirichletFenwickI64::from(ret);
    for i in (2..x).rev() {
        let ai = a.arr[i - 1] / INVS[1];
        if ai == 0 {
            continue;
        }
        ret.sparse_mul_unlimited(i, ai);
    }
    ret.into() */
    ret.adjacent_difference();
    ret.arr[0] = 0;
    for i in 2..x {
        ret.arr[i - 1] = a.arr[i - 1] / INVS[1];
    }
    ret.partial_sum();
    let mut ret_fenwick = DynamicPrefixSumI64(ret.arr, len);
    let get_index = |v| -> usize {
        if v == 0 {
            unsafe { core::hint::unreachable_unchecked() };
        } else if v <= rt {
            v - 1
        } else {
            len - (n / v)
        }
    };
    let mut sp = ret_fenwick.sum(x - 2);
    for p in (2..x).rev() {
        let sp1 = ret_fenwick.sum(p - 2);
        if sp1 == sp {
            continue;
        }

        ret_fenwick.shrink_flattened_prefix(1 + get_index(p * p));
        let w = sp - sp1;

        let lim = n / p;
        let mut prev = sp1;
        let mut i = p;
        while i <= lim / i {
            let cur = ret_fenwick.sum(i - 1);
            if cur != prev {
                ret_fenwick.add(get_index(i * p), w * (cur - prev));
                prev = cur;
            }
            i += 1;
        }
        for j in (1..=lim / i).rev() {
            let cur = ret_fenwick.sum(get_index(lim / j));
            if cur != prev {
                ret_fenwick.add(len - j, w * (cur - prev));
                prev = cur;
            }
        }
        sp = sp1;
    }

    ret_fenwick.shrink_flattened_prefix(1);
    ret_fenwick.inc(0);

    ret.arr = ret_fenwick.flatten();
    ret
}

// faster by a log factor, but more susceptible to overflow - multiplies input by a factor of 4 before reducing
/* #[must_use]
pub fn inverse_pseudo_euler_transform_fraction_i64(a: FIArrayI64) -> FIArrayI64 {
    const fn inv_odd(mut k: i64) -> i64 {
        let mut exp = (1u64 << 63) - 1;

        let mut r: i64 = 1;
        while exp > 1 {
            r = r.overflowing_mul(k).0;
            k = k.overflowing_mul(k).0;
            exp >>= 1;
        }
        r.overflowing_mul(k).0
    }
    const INVS: [i64; 5] = [0, 4, 2, inv_odd(3) << 2, 1];
    let mut a = DirichletFenwickI64::from(a);
    let rt = a.isqrt;
    let n = a.x;
    let len = a.bit.0.len();

    let mut ret = FIArrayI64::new(n);

    let x = iroot::<5>(n) + 1;
    for i in 2..x {
        let vi = a.bit.sum(i - 1) - 1;
        if vi == 0 {
            continue;
        }
        ret.arr[i - 1] = vi;
        a.sparse_mul_at_most_one(i, vi);
    }
    a.bit.dec(0);
    let a = FIArrayI64::from(a);

    for i in x..=len {
        ret.arr[i - 1] = a.arr[i - 1] * INVS[1];
    }

    let a_2 = mult_i64(&a, &a);
    let mut pow_x = x * x;
    for i in ret.get_index(pow_x) + 1..=len {
        ret.arr[i - 1] -= a_2.arr[i - 1] * INVS[2];
    }
    let a_3 = mult_i64(&a, &a_2);
    pow_x *= x;
    for i in ret.get_index(pow_x) + 1..=len {
        ret.arr[i - 1] += a_3.arr[i - 1] * INVS[3];
    }
    let a_4 = mult_i64(&a_2, &a_2);
    pow_x *= x;
    for i in ret.get_index(pow_x) + 1..=len {
        ret.arr[i - 1] -= a_4.arr[i - 1] * INVS[4];
    }

    for i in (x..len).rev() {
        ret.arr[i] -= ret.arr[i - 1];
    }

    for x in x..=rt {
        let v = ret.arr[x - 1] / INVS[1];
        if v == 0 {
            continue;
        }
        let mut e = 1;
        let mut pv = v;
        let mut px = x;
        while px <= n / x {
            e += 1;
            px *= x;
            pv *= v;

            ret[px] -= pv * INVS[e];
        }
    }
    for i in 1..x - 1 {
        ret.arr[i] += ret.arr[i - 1];
    }
    for i in x - 1..len {
        ret.arr[i] /= INVS[1];
        ret.arr[i] += ret.arr[i - 1];
    }
    ret
}
 */
#[must_use]
pub fn inverse_pseudo_euler_transform_fraction_i64(a: FIArrayI64) -> FIArrayI64 {
    const fn inv_odd(mut k: i64) -> i64 {
        let mut exp = (1u64 << 63) - 1;

        let mut r: i64 = 1;
        while exp > 1 {
            r = r.overflowing_mul(k).0;
            k = k.overflowing_mul(k).0;
            exp >>= 1;
        }
        r.overflowing_mul(k).0
    }
    const INVS: [i64; 7] = [
        0,
        4,
        2,
        inv_odd(3) << 2,
        1,
        inv_odd(5) << 2,
        inv_odd(3) << 1,
    ];
    //let start = std::time::Instant::now();
    let mut a = a;
    let n = a.x;
    let rt = a.isqrt;
    let len = a.arr.len();

    let mut ret = FIArrayI64::new(n);

    let x = iroot::<7>(n) + 1;
    // remove contributions of small primes
    let keys = FIArrayI64::keys(n).collect_vec();
    for i in 2..x {
        let vi = a.arr[i - 1] - 1;
        if vi == 0 {
            continue;
        }
        ret.arr[i - 1] = vi;
        for (j, &v) in keys.iter().enumerate().rev() {
            if v < i {
                break;
            }
            a.arr[j] -= vi * a[v / i];
        }
    }
    for e in &mut a.arr {
        *e -= 1;
    }
    let a = a;
    //println!("Finished sieving: {:?}",start.elapsed());
    // zeta now equals zeta_t - 1
    // compute log(zeta_t) using log(x + 1) = x - x^2 / 2 + x^3 / 3 - x^4 / 4 + x^5 / 5 - x^6 / 6
    // in order to not have to deal with rational numbers, we compute 4 * log(zeta_t)
    // and adjust later

    for i in x..=len {
        ret.arr[i - 1] = a.arr[i - 1] * INVS[1];
    }
    //println!("started first mul: {:?}",start.elapsed());
    let a_2 = mult_i64(&a, &a);
    //println!("finished first mul: {:?}",start.elapsed());

    let mut x_pow = x * x;

    for i in ret.get_index(x_pow) + 1..=len {
        ret.arr[i - 1] -= a_2.arr[i - 1] * INVS[2];
    }
    //println!("started second mul: {:?}",start.elapsed());
    let mut a_3 = mult_i64(&a, &a_2);
    //println!("finished second mul: {:?}",start.elapsed());

    x_pow *= x;
    for i in ret.get_index(x_pow) + 1..=len {
        ret.arr[i - 1] += a_3.arr[i - 1] * INVS[3];
    }
    //println!("started third mul: {:?}",start.elapsed());
    let a_4 = mult_i64(&a_2, &a_2);
    //println!("finished third mul: {:?}",start.elapsed());
    x_pow *= x;
    for i in ret.get_index(x_pow) + 1..=len {
        ret.arr[i - 1] -= a_4.arr[i - 1] * INVS[4];
    }
    //println!("started fourth mul: {:?}",start.elapsed());
    let a_5 = mult_sparse_i64(&a, &a_4);
    //println!("finished fourth mul: {:?}",start.elapsed());
    x_pow *= x;
    for i in ret.get_index(x_pow) + 1..=len {
        ret.arr[i - 1] += a_5.arr[i - 1] * INVS[5];
    }
    //println!("started fifth mul: {:?}",start.elapsed());
    let a_6 = mult_sparse_i64(&a_2, &a_4);
    //println!("finished fifth mul: {:?}",start.elapsed());
    x_pow *= x;
    for i in ret.get_index(x_pow) + 1..=len {
        ret.arr[i - 1] -= a_6.arr[i - 1] * INVS[6];
    }
    // correction phase: get rid of contributions of prime powers
    for i in (x..len).rev() {
        ret.arr[i] -= ret.arr[i - 1];
    }
    for x in x..=rt {
        let v = ret.arr[x - 1] / INVS[1];
        if v == 0 {
            continue;
        }
        let mut e = 1;
        let mut pv = v;
        let mut px = x;
        while px <= n / x {
            e += 1;
            px *= x;
            pv *= v;

            ret[px] -= pv * INVS[e];
        }
    }
    for i in 1..x - 1 {
        ret.arr[i] += ret.arr[i - 1];
    }
    for i in x - 1..len {
        ret.arr[i] /= INVS[1];
        ret.arr[i] += ret.arr[i - 1];
    }
    ret
}

fn mult_simple_i64(a: &FIArrayI64, b: &FIArrayI64) -> FIArrayI64 {
    unsafe { core::hint::assert_unchecked(a.x == b.x) };
    let R2 = a.isqrt;
    let n = a.x;
    let mut res = FIArrayI64::new(n);

    let collect_nonzero = |ds: &FIArrayI64| {
        let mut vec = vec![];
        if ds.arr[0] != 0 {
            vec.push((1, ds.arr[0]));
        }
        for i in 1..R2 {
            if ds.arr[i] != ds.arr[i - 1] {
                vec.push((i + 1, ds.arr[i] - ds.arr[i - 1]));
            }
        }
        vec.push((R2 + 1, 0));
        vec
    };
    let len = res.arr.len();

    let pa = collect_nonzero(a);
    let va = &pa[..pa.len() - 1];
    let pb = collect_nonzero(b);
    let vb = &pb[..pb.len() - 1];
    res.arr[0] += a.arr[0] * b.arr[0];
    for i in 2..=R2 {
        res[i * i] += (a.arr[i - 1] - a.arr[i - 2]) * (b.arr[i - 1] - b.arr[i - 2]);
    }
    let mut r = vb.len();
    let mut l = 0;
    for &(x, fx) in va {
        while pb[l].0 <= x {
            l += 1;
        }

        let Nx = n / x;
        while r > l && Nx < pb[r - 1].0 * r - 1 {
            r -= 1;
        }
        if r < l {
            r = l;
        }
        for (y, fy) in &pb[l..r] {
            res[x * y] += fx * fy;
        }
        if r != 0 && pb[r].0 <= Nx {
            res[x * pb[r].0] -= fx * b.arr[pb[r - 1].0 - 1];
        }
    }
    r = va.len();
    l = 0;
    for &(y, fy) in vb {
        while pa[l].0 <= y {
            l += 1;
        }
        let Nx = n / y;
        while r > l && Nx / pa[r - 1].0 < r - 1 {
            r -= 1;
        }
        if r < l {
            r = l;
        }

        for (x, fx) in &pa[l..r] {
            res[y * x] += fy * fx;
        }
        if r != 0 && y * pa[r].0 <= n {
            res[y * pa[r].0] -= fy * a.arr[pa[r - 1].0 - 1];
        }
    }
    res.partial_sum();
    r = vb.len();
    l = 0;
    for &(x, fx) in va {
        while pb[l].0 <= x {
            l += 1;
        }
        let Nx = n / x;
        while r > l && Nx / pb[r - 1].0 < r - 1 {
            r -= 1;
        }
        if r < l {
            r = l;
        }
        for i in (1..=Nx / pb[r].0).rev() {
            res.arr[len - i] += fx * b[Nx / i];
        }
    }
    r = va.len();
    l = 0;
    for &(y, fy) in vb {
        while pa[l].0 <= y {
            l += 1;
        }
        let Nx = n / y;
        while r > l && Nx / pa[r - 1].0 < r - 1 {
            r -= 1;
        }
        if r < l {
            r = l;
        }

        for i in (1..=Nx / pa[r].0).rev() {
            res.arr[len - i] += fy * a[Nx / i];
        }
    }

    res
}

fn mult_simple(F: &FIArray, G: &FIArray) -> FIArray {
    unsafe { core::hint::assert_unchecked(F.x == G.x) };
    let R2 = F.isqrt;
    let n = F.x;
    let mut res = FIArray::new(n);

    let len = res.arr.len();

    let f = |i: usize| {
        if i == 1 {
            F.arr[0]
        } else {
            F.arr[i - 1] - F.arr[i - 2]
        }
    };
    let g = |i: usize| {
        if i == 1 {
            G.arr[0]
        } else {
            G.arr[i - 1] - G.arr[i - 2]
        }
    };

    for i in 1..=R2 {
        res[i * i] += f(i) * g(i);
    }
    let mut r = R2;
    for x in 1..=R2 {
        let fx = f(x);
        let gx = g(x);
        let Nx = n / x;

        while r > x && Nx / r < r - 1 {
            r -= 1;
        }
        if r < x {
            r = x;
        }
        for y in x + 1..=r {
            res[x * y] += fx * g(y) + gx * f(y);
        }
        if r + 1 <= Nx {
            res[x * (r + 1)] -= fx * G.arr[r - 1] + gx * F.arr[r - 1];
        }
    }

    res.partial_sum();
    r = R2;
    for x in 1..=R2 {
        let fx = f(x);
        let gx = g(x);

        let Nx = n / x;

        while r > x && Nx / r < r - 1 {
            r -= 1;
        }
        if r < x {
            r = x;
        }
        for y in 1..=Nx / (r + 1) {
            res.arr[len - y] += fx * G[Nx / y] + gx * F[Nx / y];
        }
    }

    res
}
