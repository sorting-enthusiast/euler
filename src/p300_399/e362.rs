use itertools::Itertools;

use crate::utils::{FIArray::FIArrayU64, primes::wheel_sieve};
#[must_use]
fn sqf(x: u64) -> FIArrayU64 {
    let mut s = FIArrayU64::unit(x);
    let keys = FIArrayU64::keys(x).collect_vec().into_boxed_slice();
    let primes = wheel_sieve(x.isqrt() as _).into_boxed_slice();

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
// can optimize using fenwick trees, but the code is fast enough as is
pub fn main() {
    const N: u64 = 1e10 as _;
    const SQRT_N: u64 = N.isqrt();
    let start = std::time::Instant::now();
    let sqf = sqf(N);
    let mut fsf = FIArrayU64::eps(N);
    let keys = FIArrayU64::keys(N).collect_vec().into_boxed_slice();
    let len = keys.len();

    let mut res = 0;
    for q in 2..=SQRT_N {
        if sqf.arr[q as usize - 1] == sqf.arr[q as usize - 2] {
            continue;
        }
        for (i, &v) in keys[..=len - q as usize]
            .iter()
            .enumerate()
            .skip(q as usize - 1)
        {
            fsf.arr[i] += fsf[v / q];
        }
        res += fsf[N / q];
    }

    /*for q in SQRT_N + 1..=N {
        if squarefree[q as usize] == 0 {
            continue;
        }
        fsf[N] += fsf[N / q];
    }*/
    let mut q = SQRT_N + 1;
    while q <= N {
        let k = N / q;
        let q_max = N / k;

        let cnt = sqf[q_max] - sqf[q - 1];

        res += cnt * fsf[k];

        q = q_max + 1;
    }
    println!("res = {res}, took {:?}", start.elapsed());
}
fn test1() {
    const N: u64 = 1e2 as _;
    const SQRT_N: u64 = N.isqrt();

    let mut squarefree = vec![1; N as usize + 1];
    squarefree[0] = 0;
    for q in 2..=SQRT_N as usize {
        if squarefree[q] == 0 {
            continue;
        }
        for m in (q * q..=N as usize).step_by(q * q) {
            squarefree[m] = 0;
        }
    }
    let mut fsf = FIArrayU64::eps(N);
    let keys = FIArrayU64::keys(N).collect_vec().into_boxed_slice();
    for q in (2..=SQRT_N).rev() {
        if squarefree[q as usize] == 0 {
            continue;
        }
        for (i, &v) in keys.iter().enumerate() {
            if q > v {
                continue;
            }
            fsf.arr[i] += fsf[v / q];
        }
    }
    for q in SQRT_N + 1..=N {
        if squarefree[q as usize] == 0 {
            continue;
        }
        fsf[N] += fsf[N / q];
    }

    dbg!(fsf[N] - 1);
}
