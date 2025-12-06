use crate::utils::{FIArray::FIArrayI64, multiplicative_function_summation::totient_sum};
use std::{collections::HashMap, time::Instant};

const MOD: i64 = 1e9 as i64;
const N: usize = 510_510;
const PHI_N: usize = 92160;

const M: i64 = 1e11 as i64;

const PRIMES: [i64; 7] = [2, 3, 5, 7, 11, 13, 17];
const fn lookup_create() -> [(i64, i64); 128] {
    let mut ret = [(1, 1); 128];
    let mut i = 1_usize;
    while i < 128 {
        let mut b = i;
        while b != 0 {
            let p = b.trailing_zeros() as usize;
            ret[i].0 *= PRIMES[p];
            ret[i].1 *= PRIMES[p] - 1;
            b &= b - 1;
        }
        i += 1;
    }
    ret
}
const LOOKUP: [(i64, i64); 128] = lookup_create();
const fn mobius(pmask: i8) -> i64 {
    1 - ((pmask.count_ones() as i64 & 1) << 1)
}

fn original_s(
    n: i8,
    m: i64,
    totient_sums: &FIArrayI64,
    cache: &mut HashMap<(i8, i64), i64>,
) -> i64 {
    let phi_n = LOOKUP[n as usize].1;
    if m == 1 {
        return phi_n;
    }

    if let Some(&v) = cache.get(&(n, m)) {
        return v;
    }
    assert_ne!(n, 0); // only times when n would be 0 are unrolled

    let mut sum = 0;
    let mut d_submask = n;
    while d_submask != 0 {
        let mut inner = 0;
        let complement = n & !d_submask;
        let mut d_prime_submask = complement;
        while d_prime_submask != 0 {
            let ndp = d_prime_submask | d_submask;
            let prod = LOOKUP[ndp as usize].0;
            if m >= prod {
                inner += mobius(d_prime_submask) * original_s(ndp, m / prod, totient_sums, cache);
                inner %= MOD;
                if inner < 0 {
                    inner += MOD;
                }
            }
            d_prime_submask = (d_prime_submask - 1) & complement;
        }

        let (d, phi_d) = LOOKUP[d_submask as usize];

        if m >= d {
            inner += original_s(d_submask, m / d, totient_sums, cache);
            inner %= MOD;
        }

        let coeff = (phi_n * d) / phi_d;
        sum += coeff * inner;
        sum %= MOD;
        d_submask = (d_submask - 1) & n;
    }
    {
        let mut inner = 0;
        let complement = n;
        let mut d_prime_submask = complement;
        while d_prime_submask != 0 {
            let ndp = d_prime_submask;
            let prod = LOOKUP[ndp as usize].0;
            if m >= prod {
                inner += mobius(d_prime_submask) * original_s(ndp, m / prod, totient_sums, cache);
                inner %= MOD;
                if inner < 0 {
                    inner += MOD;
                }
            }
            d_prime_submask = (d_prime_submask - 1) & complement;
        }

        inner += totient_sums[m];
        if inner >= MOD {
            inner -= MOD;
        }

        sum += phi_n * inner;
        sum %= MOD;
    }
    cache.insert((n, m), sum);
    sum
}

fn s(n: i8, m: i64, totient_sums: &FIArrayI64) -> i64 {
    let phi_n = LOOKUP[n as usize].1;
    if m == 0 {
        0
    } else if m == 1 {
        phi_n
    } else {
        assert_ne!(n, 0); // only times when n would be 0 are unrolled
        let mut sum = 0;
        let mut d_submask = n;
        loop {
            let (d, phi_d) = LOOKUP[d_submask as usize];
            sum += ((phi_n / phi_d) * s(d_submask, m / d, totient_sums)) % MOD;
            sum %= MOD;
            d_submask = (d_submask - 1) & n;
            if d_submask == 0 {
                break;
            }
        }
        sum += phi_n * totient_sums[m];
        sum % MOD
    }
}

fn s_prime_rec(n: i8, m: i64, totient_sums: &FIArrayI64) -> i64 {
    let phi_n = LOOKUP[n as usize].1;
    if m == 1 {
        phi_n
    } else if n == 0 {
        totient_sums[m]
    } else {
        let p = PRIMES[n.trailing_zeros() as usize];
        let mut sum = (p - 1) * s_prime_rec(n & (n - 1), m, totient_sums);
        if m >= p {
            sum += s_prime_rec(n, m / p, totient_sums);
        }
        sum % MOD
    }
}

pub fn main() {
    let start = Instant::now();
    let sums = totient_sum::<MOD>(M);
    let end = start.elapsed();
    println!("all totient summations precomputed in {end:?}");
    // FIArrayI64::keys(M) = all possible values of M div n, hence all possible totient queries are already in the cache
    // note, that only queries where n is 17-smooth might occur, but inserting only those values is a pain in the ass

    let start = Instant::now();
    let sum = s_prime_rec(0x7f, M, &sums);
    let end = start.elapsed();
    println!("Recursion over primes: {sum}, took {end:?}");

    let start = Instant::now();
    let sum = s(0x7f, M, &sums);
    let end = start.elapsed();
    println!("Recursion over divisors: {sum}, took {end:?}");

    let start = Instant::now();
    let mut cache = HashMap::new();
    let sum = original_s(0x7f, M, &sums, &mut cache);
    let end = start.elapsed();
    println!("Recursion over divisors (original): {sum}, took {end:?}");
    //dbg!(cache.len());
}
/*
fn common_factor_map() -> Vec<u8> {
    let mut s = vec![0u8; N + 2];
    for (i, &p) in PRIMES.iter().enumerate() {
        for m in (p..=N).step_by(p) {
            s[m] |= 1 << i;
        }
    }
    s
}
//res = 622077701, took 2963.2892923s
//res = 546874880, took 2900.8909665s
pub fn sieve_solve() {
    // can easily precompute gcd(n,i) and phi(gcd(n,i)) bc n is product of first 7 primes
    // can use segmented sieve for phi(i) for big i
    let start = Instant::now();
    let mut sum = 0;
    let s = common_factor_map();
    let mut totients = (0..=N).collect_vec();
    dbg!((M / N, M % N));
    let mut primes = vec![];
    for p in 2..=N {
        if totients[p] == p {
            totients[p] -= 1;
            let mut multiple = p << 1;
            while multiple <= N {
                totients[multiple] /= p;
                totients[multiple] *= p - 1;
                multiple += p;
            }
            if p <= M.isqrt() {
                primes.push(p);
            }
        }
    }
    for i in 1..=N {
        let (gcd, phi_gcd) = LOOKUP[s[i] as usize];
        assert_eq!(totients[i].checked_mul(gcd).unwrap() % phi_gcd, 0, "damn");
        sum += (totients[i].checked_mul(gcd).unwrap() / phi_gcd) % MOD;
        if sum >= MOD {
            sum -= MOD;
        }
        assert!(sum < MOD);
    }
    dbg!(sum);
    let mut l = N + 1;
    let mut large_prime_factor = vec![0; N + 1];
    let primes = primes;
    while l < M {
        let r = M.min(l + N);
        totients
            .iter_mut()
            .enumerate()
            .for_each(|(p, t)| *t = l + p);
        large_prime_factor
            .iter_mut()
            .enumerate()
            .for_each(|(p, t)| *t = l + p);
        for &p in &primes {
            if p * p > r {
                break;
            }
            let mut multiple = l.next_multiple_of(p);
            while multiple <= r {
                totients[multiple - l] /= p;
                totients[multiple - l] *= p - 1;
                let mut k = large_prime_factor[multiple - l];
                while k % p == 0 {
                    k /= p;
                }
                large_prime_factor[multiple - l] = k;
                multiple += p;
            }
        }
        for i in l..=r {
            if totients[i - l] == i {
                sum += (i - 1) % MOD;
            } else {
                let (gcd, phi_gcd) = LOOKUP[s[i - l + 1] as usize];
                let large_prime = large_prime_factor[i - l];
                if large_prime != 1 {
                    totients[i - l] /= large_prime;
                    totients[i - l] *= large_prime - 1;
                }
                assert_eq!(
                    (totients[i - l].checked_mul(gcd).unwrap()) % phi_gcd,
                    0,
                    "damn"
                );
                sum += (totients[i - l].checked_mul(gcd).unwrap() / phi_gcd) % MOD;
            }
            if sum >= MOD {
                sum -= MOD;
            }
            assert!(sum < MOD);
        }
        if ((l - 1) / N).trailing_zeros() >= 7 {
            dbg!(start.elapsed());
            dbg!((l, r));
            dbg!(r as f64 / M as f64);
            dbg!(M as f64 / r as f64);
        }
        if r == M {
            dbg!(start.elapsed());
            dbg!((l, r));
            dbg!((r - l, N));
            dbg!(r as f64 / M as f64);
            dbg!(M as f64 / r as f64);
        }
        l = r;
    }
    sum = sum.checked_mul(PHI_N).unwrap() % MOD;
    let end = start.elapsed();
    println!("res = {sum}, took {end:?}");
    //assert_eq!(45_480_596_821_125_120_usize % MOD, sum);
}
 */
