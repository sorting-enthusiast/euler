use itertools::Itertools;

use crate::utils::FIArray::FIArray;
use crate::utils::primes::prime_sieves::sift;
// 501
// ez: count composites with either 3 distinct prime factors,
// or 2 distinct prime factors, with one of them being of multplicity 3 and the other of 1,
// or p^7

pub fn main() {
    const N: usize = 1e12 as _;
    let start = std::time::Instant::now();
    let primes = sift(N.isqrt() as u64);
    let mut pi = FIArray::new(N);
    let keys = FIArray::keys(N).collect_vec();

    unsafe { core::hint::assert_unchecked(pi.arr.len() == keys.len()) };
    for (i, v) in keys.iter().enumerate() {
        pi.arr[i] = (v + 1) >> 1;
    }
    pi.arr[0] = 0;
    pi.arr[2] = 2;
    for &p in &primes[1..] {
        let p = p as usize;
        let sp = pi.arr[p - 2];

        let pdiv = strength_reduce::StrengthReducedUsize::new(p);
        for (i, &v) in keys.iter().enumerate().rev() {
            if v < p * p {
                break;
            }
            pi.arr[i] -= pi[v / pdiv] - sp;
        }
    }
    let mut sum = primes.partition_point(|p| p.pow(7) <= N as u64);
    for (i, &p) in primes.iter().enumerate() {
        let p = p as usize;
        if p * p * p > N {
            break;
        }
        sum += pi[N / (p * p * p)] - usize::from(p <= N / (p * p * p));

        for &q in &primes[i + 1..] {
            let q = q as usize;

            if p * q * q >= N {
                break;
            }
            sum += pi[N / (p * q)] - pi.arr[q - 1];
        }
    }
    println!("res = {sum}, took {:?}", start.elapsed());
}
