use crate::utils::multiplicative_function_summation::mobius_sieve;

const N: i64 = 1e12 as _;
const SQRT_N: i64 = N.isqrt();

const fn icbrt(x: i64) -> i64 {
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

// runs in O(n^2/3)
fn inner(n: i64) -> i64 {
    // inner2 is around 4 times faster than inner1
    fn inner1(n: i64) -> i64 {
        #[inline(always)]
        fn s(lim: i64, x1: i64, x2: i64) -> i64 {
            let mut ret = 0;
            let mut q = x1;
            while q <= lim.min(x2) {
                let k = lim / q;
                let q_next = 1 + (lim / k).min(x2);

                ret += (q_next - q) * k;

                q = q_next;
            }
            ret
        }
        let mut ret = 0;
        let nsqrt = n.isqrt();
        for i in 1..=nsqrt {
            let lim = n / i;
            ret += s(lim, i + 1, 2 * i - 1);
        }
        ret
    }

    fn inner2(n: i64) -> i64 {
        #[inline(always)]
        fn s(lim: i64, x1: i64, x2: i64) -> i64 {
            let mut ret = 0;
            let mut prev = x1 - 1;
            let z = lim / x1;
            let zL = lim / x2;

            for z in (zL + 1..=z).rev() {
                let curr = lim / z;
                ret += z * (curr - prev);
                prev = curr;
            }

            ret += zL * (x2 - prev);
            ret
        }

        let mut ret = 0;
        let nsqrt = n.isqrt();
        let ncbrt = icbrt(n);
        for i in 1..=ncbrt {
            let lim = n / i;
            for j in i + 1..2 * i {
                ret += lim / j;
            }
        }
        for i in ncbrt + 1..=nsqrt {
            let lim = n / i;
            ret += s(lim, i + 1, 2 * i - 1);
        }
        ret
    }
    inner2(n)
}

pub fn main() {
    solve();
    let start = std::time::Instant::now();
    let mob = mobius_sieve(SQRT_N as usize + 1);
    let mut res = inner(N);
    for i in 2..=SQRT_N {
        if mob[i as usize] != 0 {
            res += i64::from(mob[i as usize]) * inner(N / (i * i));
        }
    }
    println!("res = {res}, took {:?}", start.elapsed());
}

// O(n^2/3) time, O(n^1/3) space
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
            res += mertens_small[d as usize] * inner(N / (d * d));
        }
        mertens_small[d as usize] += mertens_small[d as usize - 1];
    }
    let mut small_diff = vec![0; I as usize - 1].into_boxed_slice();
    for i in 1..I {
        small_diff[i as usize - 1] = inner(i);
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
