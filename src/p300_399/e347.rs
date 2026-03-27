use crate::utils::primes::wheel_sieve;

const N: usize = 1e7 as _;
const SQRT_N: usize = N.isqrt();
pub fn main() {
    let primes = wheel_sieve(N as _);
    let mut res = 0;
    for (i, &p) in primes.iter().enumerate() {
        let p = p as usize;
        let np = N / p;
        if p >= np {
            break;
        }
        for &q in &primes[i + 1..] {
            let q = q as usize;
            if q > np {
                break;
            }
            res += find_largest(p, q);
        }
    }
    dbg!(res);
}
fn find_largest(p: usize, q: usize) -> usize {
    let mut max = 0;
    for i in 1..=N.ilog(p) {
        let pp = p.pow(i);
        for j in 1..=(N / pp).ilog(q) {
            let ppqq = pp * q.pow(j);
            max = max.max(ppqq);
        }
    }
    max
}
