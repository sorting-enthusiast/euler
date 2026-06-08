use crate::utils::primes::wheel_sieve;

const MOD: u32 = 1e9 as u32 + 7;
const fn mulmod(a: u32, b: u32) -> u32 {
    ((a as u64 * b as u64) % MOD as u64) as u32
}
fn SubsetZetaTransform(v: &mut [u32]) {
    let n = v.len();
    let mut j = 1;
    while j < n {
        for i in 0..n {
            if i & j != 0 {
                v[i] += v[i ^ j];
                if v[i] >= MOD {
                    v[i] -= MOD;
                }
            }
        }
        j <<= 1
    }
}
fn SubsetMobiusTransform(v: &mut [u32]) {
    let n = v.len();
    let mut j = 1;
    while j < n {
        for i in 0..n {
            if i & j != 0 {
                v[i] += MOD - v[i ^ j];
                if v[i] >= MOD {
                    v[i] -= MOD;
                }
            }
        }
        j <<= 1
    }
}
pub fn main() {
    let start = std::time::Instant::now();
    let n = 1e6 as usize + 1;
    let mut exp = 999983 as usize;
    let primes = wheel_sieve(n as u64 - 1);
    let mut x = vec![0; n].into_boxed_slice();
    for &p in &primes {
        x[p as usize] = 1;
    }
    let mut r = vec![0; n].into_boxed_slice();
    r[0] = 1;

    SubsetZetaTransform(&mut r);
    SubsetZetaTransform(&mut x);
    while exp > 1 {
        if exp & 1 != 0 {
            for i in 0..n {
                r[i] = mulmod(r[i], x[i]);
            }
        }
        for i in 0..n {
            x[i] = mulmod(x[i], x[i]);
        }
        exp >>= 1;
    }
    for i in 0..n {
        r[i] = mulmod(r[i], x[i]);
    }
    SubsetMobiusTransform(&mut r);
    let res = primes
        .into_iter()
        .map(|p| r[p as usize] as u64)
        .sum::<u64>()
        % MOD as u64;
    println!("res = {res}, took {:?}", start.elapsed());
}
