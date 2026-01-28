use rayon::iter::{IntoParallelIterator, ParallelIterator};

use crate::utils::{
    FIArray::FIArrayI64,
    multiplicative_function_summation::{sum_n_i64, totient_sieve},
};
const N: usize = 1e11 as _;
const MOD: i64 = 1e8 as i64;
const N1: i64 = N as i64 + 1;
const N2: i64 = (N1 % MOD * N1 % MOD) % MOD;
const fn powmod(mut x: i64, mut exp: i64) -> i64 {
    if exp == 0 {
        return 1;
    }
    let mut r: i64 = 1;
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
//https://oeis.org/A119437
pub fn main() {
    let start = std::time::Instant::now();

    let mut res = const { (powmod(powmod(2, N1), N1) - N2 - 1) % MOD }; // all sets of points minus empty set and singletons
    if res < 0 {
        res += MOD;
    }
    dbg!(res);
    res -= s(); // subtract all collinear sets
    if res < 0 {
        res += MOD;
    }
    println!("{res}, took {:?}", start.elapsed());
}
fn s() -> i64 {
    let f1 = sum_f(N1 as _);
    let f2 = sum_f2(N1 as _);
    let f3 = sum_f3(N1 as _);
    let f = |k| {
        let mut s = (N2 * (2 * f1[(N1 / k) as _] - 1)) % MOD
            + (powmod(k, 2) * f3[(N1 / k) as _]) % MOD
            - ((const { N1 % MOD } * (k % MOD)) % MOD * (3 * f2[(N1 / k) as _] - 1)) % MOD;
        s %= MOD;
        if s < 0 {
            s += MOD;
        }
        s
    };
    let S = |k_lo, k_hi| f(k_hi + 1) - f(k_hi) - (f(k_lo) - f(k_lo - 1));

    /* let helper = || {
        let mut sum = 0;
        let mut l = 3;
        while l <= N1 {
            let div = N1 / l;
            let r = N1 / div;
            //sum += MOD - (s[div] * (id[r] + MOD - id[l - 1]) % MOD) % MOD;
            if sum >= MOD {
                sum -= MOD;
            }
            l = r + 1;
        }
        sum
    }; */
    dbg!(S(3, 3));
    dbg!(S(4, 4));

    let coeff = |i| {
        let mut ret = powmod(2, i) - i % MOD - sum_n_i64::<MOD>(i as usize - 1) - 1;
        ret %= MOD;
        if ret < 0 {
            ret += MOD;
        }
        ret
    };
    // TODO: optimise sum using sqrt trick
    let mut s = (3..=N1)
        .into_par_iter()
        .map(|i| (coeff(i) * S(i, i)) % MOD)
        .sum::<i64>();
    s %= MOD;
    if s < 0 {
        s += MOD;
    }
    s <<= 1;
    //dbg!(helper() - S(3, N1));
    dbg!(s % MOD)
}

const fn sum_squares(x: i64) -> i64 {
    let x = (x % (6 * MOD)) as u128;
    (((x * (x + 1) * (2 * x + 1)) / 6) % MOD as u128) as i64
}
const fn sum_cubes(x: usize) -> i64 {
    let s = sum_n_i64::<MOD>(x);
    (s * s) % MOD
}

// (N phi) = (N (N * mu)) = (N^2) * (N mu) =>
// N * (N phi) = N^2
fn sum_f2(x: usize) -> FIArrayI64 {
    let y = if x > 15 {
        (1e8 as usize).min((x as f64).powf(2. / 3.) as usize >> 1)
    } else {
        x as usize
    };
    let mut small = totient_sieve(y + 1);
    for i in 2..=y {
        small[i] %= MOD;
        small[i] *= i as i64 % MOD;
        small[i] %= MOD;
        small[i] += small[i - 1];
        small[i] %= MOD;
    }

    let mut ret = FIArrayI64::new(x);

    for (i, v) in FIArrayI64::keys(x).enumerate() {
        if v as usize <= y {
            ret.arr[i] = small[v as usize];
            continue;
        }
        let mut f_v = (sum_squares(v as _) - sum_n_i64::<MOD>(v)) % MOD;
        let vsqrt = v.isqrt();
        for i in 2..=vsqrt {
            f_v -= ((i as i64 % MOD) * ret[v / i]) % MOD;
            f_v -= ((ret.arr[i as usize - 1] - ret.arr[i as usize - 2]) % MOD
                * sum_n_i64::<MOD>(v / i))
                % MOD;
            f_v %= MOD;
        }
        f_v += (ret[vsqrt] * sum_n_i64::<MOD>(vsqrt)) % MOD;
        f_v %= MOD;
        if f_v < 0 {
            f_v += MOD;
        }
        ret.arr[i] = f_v;
    }

    ret
}

// u * phi = N
fn sum_f(x: usize) -> FIArrayI64 {
    let y = if x > 15 {
        (1e8 as usize).min((x as f64).powf(2. / 3.) as usize >> 1)
    } else {
        x as usize
    };
    let mut small_phi = totient_sieve(y + 1);
    for i in 2..=y {
        small_phi[i] %= MOD;
        small_phi[i] += small_phi[i - 1];
        small_phi[i] %= MOD;
    }
    let mut Phi = FIArrayI64::new(x);

    for (ind, v) in FIArrayI64::keys(x).enumerate() {
        if v as usize <= y {
            Phi.arr[ind] = small_phi[v as usize];
            continue;
        }
        let vsqrt = v.isqrt();

        let mut phi_v = (sum_n_i64::<MOD>(v) + MOD - v as i64 % MOD) % MOD;
        for i in 2..=vsqrt {
            phi_v -= ((Phi.arr[i as usize - 1] - Phi.arr[i as usize - 2]) % MOD
                * ((v / i) as i64 % MOD))
                % MOD;
            phi_v -= Phi[v / i];
            phi_v %= MOD;
        }
        phi_v += (Phi[vsqrt] * vsqrt as i64) % MOD;
        phi_v %= MOD;
        if phi_v < 0 {
            phi_v += MOD;
        }
        Phi.arr[ind] = phi_v;
    }
    Phi
}

// (N^2 phi) = (N^2 (N * mu)) = (N^3) * (N^2 mu) =>
// (N^2) * (N^2 phi) = N^3
fn sum_f3(x: usize) -> FIArrayI64 {
    let y = if x > 15 {
        (1e8 as usize).min((x as f64).powf(2. / 3.) as usize >> 1)
    } else {
        x as usize
    };
    let mut small = totient_sieve(y + 1);
    for i in 2..=y {
        small[i] %= MOD;
        small[i] *= powmod(i as i64, 2);
        small[i] %= MOD;
        small[i] += small[i - 1];
        small[i] %= MOD;
    }

    let mut ret = FIArrayI64::new(x);

    for (i, v) in FIArrayI64::keys(x).enumerate() {
        if v as usize <= y {
            ret.arr[i] = small[v as usize];
            continue;
        }
        let mut f_v = (sum_cubes(v as _) - sum_squares(v as _)) % MOD;
        let vsqrt = v.isqrt();
        for i in 2..=vsqrt {
            f_v -= (powmod(i as _, 2) * ret[v / i]) % MOD;
            f_v -= ((ret.arr[i as usize - 1] - ret.arr[i as usize - 2]) % MOD
                * sum_squares((v / i) as _))
                % MOD;
            f_v %= MOD;
        }
        f_v += (ret[vsqrt] * sum_squares(vsqrt as _)) % MOD;
        f_v %= MOD;
        if f_v < 0 {
            f_v += MOD;
        }
        ret.arr[i] = f_v;
    }

    ret
}
