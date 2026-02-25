use crate::{
    div_i64, mult_i64, mult_sparse_i64,
    utils::{
        FIArray::{DirichletFenwickI64, FIArrayI64},
        math::iroot,
        multiplicative_function_summation::mertens,
        primes::wheel_sieve,
    },
};
const fn sum_squares(x: i128) -> i128 {
    (x * (x + 1) * (2 * x + 1)) / 6
}
const fn sum_n(x: i128) -> i128 {
    if x & 1 == 0 {
        (x / 2) * (x + 1)
    } else {
        ((x + 1) / 2) * x
    }
}
const N: i128 = 1e12 as _;
const SQRT_N: i128 = N.isqrt();

// count coprime pairs * 3 for x,y,z planes, using the fact that the count is 2*totient sum - 1
// add contribution of all coprime triplets a,b,c, using generalized mobius inversion and dirichlet hyperbola method on
// n^3 = sum_{d<=n} C(n/d)
pub fn main() {
    let start = std::time::Instant::now();
    let m = mertens(N as _);

    let len = m.arr.len();

    let mut r = i128::from(m[N as _]) + N.pow(3) - i128::from(m[SQRT_N as _]) * SQRT_N.pow(3);
    let mut tot = i128::from(m[N as _]) + sum_n(N) - sum_n(SQRT_N) * i128::from(m[SQRT_N as _]);
    for i in 2..=SQRT_N as usize {
        r += i128::from(m.arr[i - 1] - m.arr[i - 2]) * (N / i as i128).pow(3);
        r += (3 * i * i - 3 * i + 1) as i128 * i128::from(m.arr[len - i]);

        tot += i as i128 * i128::from(m.arr[len - i]);
        tot += i128::from(m.arr[i - 1] - m.arr[i - 2]) * sum_n(N / i as i128);
    }
    let res = r + 3 + 3 * (2 * tot - 1);
    println!("res = {res}, took {:?}", start.elapsed());
    alt();
}
// alternate solution
fn alt() {
    let start = std::time::Instant::now();
    let m = mertens(N as _);

    let len = m.arr.len();

    let mut j2 =
        i128::from(m[N as _]) + sum_squares(N) - sum_squares(SQRT_N) * i128::from(m[SQRT_N as _]);
    let mut j1 = i128::from(m[N as _]) + sum_n(N) - sum_n(SQRT_N) * i128::from(m[SQRT_N as _]);
    for i in 2..=SQRT_N as usize {
        let mu_i = i128::from(m.arr[i - 1] - m.arr[i - 2]);

        j2 += (i * i) as i128 * i128::from(m.arr[len - i]);
        j2 += mu_i * sum_squares(N / i as i128);

        j1 += i as i128 * i128::from(m.arr[len - i]);
        j1 += mu_i * sum_n(N / i as i128);
    }
    let res = 1 + 3 * (j2 + j1);
    println!("res = {res}, took {:?}", start.elapsed());
}
// 1e15:
// 1e16: 242.1647739s
pub fn solve() {
    const N: usize = 1e16 as _;
    let start = std::time::Instant::now();
    let mob = {
        let mut zeta = DirichletFenwickI64::zeta(N);
        let lim = iroot::<7>(N) + 1;
        let mut primes = vec![];
        for p in 2..lim {
            if zeta.get_bucket_prefix(p - 1) == 1 {
                continue;
            }
            primes.push(p);
            zeta.sparse_mul_at_most_one(p, 1);
        }
        let zeta_lim = FIArrayI64::from(zeta);
        let mut mu = DirichletFenwickI64::from(div_i64(&FIArrayI64::eps(N), &zeta_lim));

        for &p in primes.iter().rev() {
            mu.sparse_mul_at_most_one(p, 1);
        }
        FIArrayI64::from(mu)
    };
    dbg!(start.elapsed());
    let rt_n = mob.isqrt;
    let len = mob.arr.len();

    let mut ret = (N as i128 + 1).pow(3) - 1 + 7 * mob.arr[len - 1] as i128
        - mob.arr[rt_n - 1] as i128 * ((rt_n as i128 + 1).pow(3) - 1);
    for i in 2..=rt_n {
        ret += (mob.arr[i - 1] - mob.arr[i - 2]) as i128 * (((N / i) as i128 + 1).pow(3) - 1)
            + (3 * i as i128 * (i as i128 + 1) + 1) * mob.arr[len - i] as i128;
    }
    println!("res = {ret}, took {:?}", start.elapsed());
}

pub fn solve_alt() {
    const N: usize = 1e16 as _;
    let start = std::time::Instant::now();
    let mob = {
        let mut zeta = DirichletFenwickI64::zeta(N);
        let lim = iroot::<8>(N) + 1;
        let primes = wheel_sieve(zeta.isqrt as u64);
        for &p in &primes {
            let p = p as usize;
            if p >= lim {
                break;
            }
            zeta.sparse_mul_at_most_one(p, 1);
        }
        let zeta_lim = FIArrayI64::from(zeta);

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
        let mut mu = DirichletFenwickI64::from(mu);
        let i = primes.partition_point(|&p| p < lim as u64);
        for &p in primes[..i].iter().rev() {
            mu.sparse_mul_at_most_one(p as _, 1);
        }
        FIArrayI64::from(mu)
    };
    dbg!(start.elapsed());
    let rt_n = mob.isqrt;
    let len = mob.arr.len();

    let mut ret = (N as i128 + 1).pow(3) - 1 + 7 * mob.arr[len - 1] as i128
        - mob.arr[rt_n - 1] as i128 * ((rt_n as i128 + 1).pow(3) - 1);
    for i in 2..=rt_n {
        ret += (mob.arr[i - 1] - mob.arr[i - 2]) as i128 * (((N / i) as i128 + 1).pow(3) - 1)
            + (3 * i as i128 * (i as i128 + 1) + 1) * mob.arr[len - i] as i128;
    }
    println!("res = {ret}, took {:?}", start.elapsed());
}
