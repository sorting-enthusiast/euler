use std::time::Instant;

use itertools::Itertools;

use crate::utils::powerful_numbers::PowerfulExtSkipZero;
const N: i64 = 1e18 as i64;
const MOD: i64 = 1e9 as i64 + 7;
const fn powmod(mut x: i64, mut exp: i64) -> i64 {
    if exp == 0 {
        return 1;
    }
    let mut r = 1;
    x %= MOD;
    while exp > 1 {
        if exp & 1 == 1 {
            r = (r * x) % MOD;
        }
        x = (x * x) % MOD;
        exp >>= 1;
    }
    (r * x) % MOD
}

// powerful number trick ftw: g(p) = 1(p) = 1
// g(p^e) = p^(2 * floor(e/2))
// h(p^e) = (f*g^-1)(p^e) = 0 if e odd, p^(e-2) * (p^2 - 1) if e even
pub fn main() {
    overkill();
    let start = Instant::now();

    let h = |p, e| {
        if e & 1 == 0 {
            (powmod(p, e - 2) * ((p * p - 1) % MOD)) % MOD
        } else {
            0
        }
    };
    let mut sum = 0;
    for (n, hn) in PowerfulExtSkipZero::<_, MOD>::new(N, h) {
        sum += hn * ((N / n) % MOD);
        sum %= MOD;
    }

    let end = start.elapsed();
    println!("res = {sum}, took {end:?}");
}
// O(n^4/9) time O(n^1/3) space
fn overkill() {
    const N: usize = self::N as _;
    const MOD: usize = self::MOD as _;
    const fn sum_squares(x: usize) -> usize {
        let x = (x % (6 * MOD)) as u128;
        (((x * (x + 1) * (2 * x + 1)) / 6) % MOD as u128) as usize
    }
    const fn icbrt(x: usize) -> usize {
        let mut rt = 1 << (1 + x.ilog2().div_ceil(3));
        let mut x_div_rt2 = (x / rt) / rt;
        while rt > x_div_rt2 {
            rt = ((rt << 1) + x_div_rt2) / 3;
            x_div_rt2 = (x / rt) / rt;
        }
        rt
    }
    fn sqf_sieve(n: usize) -> Vec<usize> {
        unsafe { core::hint::assert_unchecked(n >= 1) };
        let mut sqf = vec![1; n];
        sqf[0] = 0;
        if n < 2 {
            return sqf;
        }
        let sqrtn = n.isqrt();
        let mut d2 = 1;
        for d in 2..=sqrtn {
            d2 += (d << 1) - 1;
            if sqf[d2] == 0 {
                continue;
            }
            for m in (d2..n).step_by(d2) {
                sqf[m] = 0;
            }
        }
        sqf
    }
    const B: usize = icbrt(N);
    const A: usize = N / (B * B);
    const SQRT_N: usize = N.isqrt();
    let start = Instant::now();

    let sqrts = (2..=A)
        .map(|i| (N / i).isqrt())
        .collect_vec()
        .into_boxed_slice(); // precompute once, used very often
    let mut sqf_small = sqf_sieve(A + 1);
    let mut res = sum_squares(SQRT_N);
    for i in 2..=A {
        if sqf_small[i] == 1 {
            res += sum_squares(sqrts[i - 2]);
            if res >= MOD {
                res -= MOD;
            }
        }
        sqf_small[i] += sqf_small[i - 1];
    }
    let mut sqf_big = vec![0; B - 1].into_boxed_slice(); // indexed by denominator

    for d in (1..B).rev() {
        let v = N / (d * d);
        let b = icbrt(v);
        let a = v / (b * b);

        let mut sqf = v + sqf_small[a] * b - SQRT_N / d;
        for i in 2..=a {
            sqf -= (sqf_small[i] - sqf_small[i - 1]) * sqrts[i - 2] / d; //(v / i).isqrt();
        }
        for i in 2..=b {
            sqf -= if i * d < B {
                sqf_big[(i * d) - 1]
            } else {
                sqf_small[v / (i * i)]
            };
        }
        sqf_big[d - 1] = sqf;
        res += (d * d) % MOD * (sqf % MOD);
        res %= MOD;
    }
    res += MOD - (sqf_small[A] % MOD * sum_squares(B - 1)) % MOD;
    if res >= MOD {
        res -= MOD;
    }

    let end = start.elapsed();
    println!("res = {res}, took {end:?}");
}
