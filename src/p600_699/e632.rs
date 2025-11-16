use itertools::Itertools;

use crate::utils::prime_sieves::sift;
const N: u64 = 1e19 as _;
const SQRT: u64 = N.isqrt();
const CBRT: u64 = icbrt(N);
const MOD: u64 = 1e9 as u64 + 7;
const fn icbrt(x: u64) -> u64 {
    let mut rt = 0;
    let mut rt_squared = 0;
    let mut rt_cubed = 0;
    while rt_cubed <= x {
        rt += 1;
        rt_squared += 2 * rt - 1;
        rt_cubed += 3 * rt_squared - 3 * rt + 1;
    }
    rt - 1
}
pub fn main() {
    let start = std::time::Instant::now();
    let primes = sift(SQRT).into_boxed_slice();
    println!("Generated primes: {:?}", start.elapsed());
    let large_keys = (1..CBRT)
        .map(|d| N / (d * d))
        .chain((CBRT != N / (CBRT * CBRT)).then_some(N / (CBRT * CBRT)))
        .collect_vec()
        .into_boxed_slice();
    let mut large_sqf = vec![0; large_keys.len()].into_boxed_slice();
    let mut small_sqf = vec![0; CBRT as usize].into_boxed_slice();

    for v in 1..=CBRT {
        small_sqf[v as usize - 1] = v;
    }
    for (i, &v) in large_keys.iter().enumerate() {
        large_sqf[i] = v;
    }
    for &p in &primes {
        for (i, &v) in large_keys.iter().enumerate() {
            if v < p * p {
                break;
            }
            large_sqf[i] -= if v / (p * p) <= CBRT {
                small_sqf[(v / (p * p)) as usize - 1]
            } else {
                large_sqf[((i + 1) * (p as usize)) - 1]
            };
        }
        for v in (1..=CBRT).rev() {
            if v < p * p {
                break;
            }
            small_sqf[v as usize - 1] -= small_sqf[(v / (p * p)) as usize - 1];
        }
    }
    println!("Counted squarefree numbers: {:?}", start.elapsed());
    let mut sum = 0;
    let mut res = 1;
    for k in 0.. {
        let c_k = squarefree(N, 1, k, &primes, &small_sqf, &large_sqf);
        if c_k == 0 {
            break;
        }
        sum += c_k;
        println!("C{k} = {c_k}");
        res *= c_k % MOD;
        res %= MOD;
    }
    println!("res = {res}, took {:?}", start.elapsed());
    dbg!(N - sum);
}
fn squarefree(
    limit: u64,
    div: u64,
    k: u8,
    primes: &[u64],
    small_sqf: &[u64],
    large_sqf: &[u64],
) -> u64 {
    if k == 0 {
        return if limit <= CBRT {
            small_sqf[limit as usize - 1]
        } else {
            large_sqf[div as usize - 1]
        };
    }
    let mut res = 0;
    for (i, &p) in primes.iter().enumerate() {
        let pp = p * p;
        if pp > limit {
            break;
        }
        let mut new_lim = limit / pp;
        let mut new_div = div * p;
        loop {
            res += squarefree(
                new_lim,
                new_div,
                k - 1,
                &primes[i + 1..],
                small_sqf,
                large_sqf,
            );
            new_lim /= pp;
            new_div *= p;
            if new_lim == 0 {
                break;
            }
        }
    }
    res
}
