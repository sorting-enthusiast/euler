use crate::utils::{fast_divisor_sums::icbrt, primes::wheel_sieve};

fn r2(t: i64) -> i64 {
    let mut res = 0;
    let tsqrt = t.isqrt();
    for k in 1..=tsqrt {
        let tk = t / k;
        res += [0, 1, 1, 0][(tk & 3) as usize];
        res += [0, 1, 0, -1][(k & 3) as usize] * tk;
    }
    res - [0, 1, 1, 0][(tsqrt & 3) as usize] * tsqrt
}
const N: i64 = (47 * 47) as _;
const SQRT_N: i64 = N.isqrt();
fn mobius_sieve(n: usize) -> Vec<i8> {
    let mut res = vec![0; n];
    if n < 2 {
        return res;
    }
    let mut composite = crate::utils::bit_array::BitArray::zeroed(n);
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

pub fn main() {
    dbg!(r2(8 * 8));
    let mob = mobius_sieve(SQRT_N as usize + 1);
    dbg!(&mob[1..]);
    let mut res = r2(N);
    for i in 2..=SQRT_N {
        if mob[i as usize] != 0 {
            res += i64::from(mob[i as usize]) * r2(N / (i * i));
        }
    }
    dbg!(res);
    solve();
    dbg!(
        wheel_sieve(SQRT_N as _)
            .into_iter()
            .filter(|&p| p & 3 == 3)
            .count()
    );
}
fn solve() {
    fn mobius_sieve(n: usize) -> Vec<i64> {
        let mut res = vec![0; n];
        if n < 2 {
            return res;
        }
        let mut composite = crate::utils::bit_array::BitArray::zeroed(n);
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

    const I: i64 = icbrt(N);
    const D: i64 = (N / I).isqrt();
    let start = std::time::Instant::now();
    let mut res = 0;
    let mut mertens_small = mobius_sieve(D as usize + 1);
    for d in 1..=D {
        if mertens_small[d as usize] != 0 {
            res += mertens_small[d as usize] * r2(N / (d * d));
        }
        mertens_small[d as usize] += mertens_small[d as usize - 1];
    }
    let mut small_diff = vec![0; I as usize - 1].into_boxed_slice();
    for i in 1..I {
        small_diff[i as usize - 1] = r2(i);
    }
    res -= mertens_small[D as usize] * small_diff[I as usize - 2];

    for i in (2..I).rev() {
        small_diff[i as usize - 1] -= small_diff[i as usize - 2];
    }

    let mut mertens_big = vec![0; I as usize - 1].into_boxed_slice(); // indexed by denominator
    for i in (1..I).rev() {
        let v = (N / i).isqrt() as usize;
        let vsqrt = v.isqrt();
        let mut m = 1 - v as i64 + vsqrt as i64 * mertens_small[vsqrt];
        for d in 2..=vsqrt {
            m -= if v / d <= D as usize {
                mertens_small[v / d]
            } else {
                mertens_big[i as usize * d * d - 1]
            };
            m -= (mertens_small[d] - mertens_small[d - 1]) * (v / d) as i64;
        }
        mertens_big[i as usize - 1] = m;
        res += small_diff[i as usize - 1] * m;
    }

    println!("res = {res}, took {:?}", start.elapsed());
}
