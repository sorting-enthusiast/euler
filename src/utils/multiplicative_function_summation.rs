use itertools::Itertools;

use crate::{
    mult_correction, mult_i64, mult_sparse_i64,
    p300_399::e362::mult,
    utils::{
        FIArray::{
            DirichletFenwick, DirichletFenwickI64, FIArray, FIArrayI32, FIArrayI64, FIArrayI128,
            FIArrayIsize, FIArrayU32, FIArrayU64, FIArrayU128, FIArrayUsize,
        },
        bit_array::BitArray,
        math::iroot,
        primes::wheel_sieve,
    },
};
#[must_use]
pub fn totient_sieve(n: usize) -> Vec<i64> {
    let mut res = vec![0; n];
    if n < 2 {
        return res;
    }
    let mut composite = BitArray::zeroed(n);
    let mut primes = vec![];
    res[1] = 1;
    for i in 2..n {
        if !composite.get(i) {
            if (i << 1) < n {
                primes.push(i);
            }
            res[i] = i as i64 - 1;
        }
        for &p in &primes {
            if i * p >= n {
                break;
            }
            composite.set(i * p);
            if i % p == 0 {
                res[i * p] = res[i] * p as i64;
                break;
            }
            res[i * p] = res[i] * (p as i64 - 1);
        }
    }
    res
}
#[must_use]
pub fn mobius_sieve(n: usize) -> Vec<i8> {
    let mut res = vec![0; n];
    if n < 2 {
        return res;
    }
    let mut composite = BitArray::zeroed(n);
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
#[must_use]
pub fn mobius_sieve_i16(n: usize) -> Vec<i16> {
    let mut res = vec![0; n];
    if n < 2 {
        return res;
    }
    let mut composite = BitArray::zeroed(n);
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
#[must_use]
pub fn divisor_sieve(n: usize) -> Vec<i64> {
    let mut res = vec![0; n];
    if n < 2 {
        return res;
    }
    let mut composite = BitArray::zeroed(n);
    let mut pow = vec![0; n];
    let mut primes = vec![];
    res[1] = 1;
    for i in 2..n {
        if !composite.get(i) {
            if i << 1 < n {
                primes.push(i);
            }
            res[i] = 2;
            pow[i] = i;
        }
        for &p in &primes {
            if i * p >= n {
                break;
            }
            composite.set(i * p);
            if i.is_multiple_of(p) {
                res[i * p] = res[i] + res[i / pow[i]];
                pow[i * p] = pow[i] * p;
                break;
            }
            res[i * p] = res[i] << 1;
            pow[i * p] = p;
        }
    }
    res
}

#[must_use]
pub fn totient_sum<const MOD: i64>(x: usize) -> FIArrayI64 {
    let y = if x > 1023 {
        (1e8 as usize).min((x as f64).powf(2. / 3.) as usize >> 1)
    } else {
        x as usize
    };
    let mut small_phi = totient_sieve(y + 1);
    for i in 2..=y {
        small_phi[i] += small_phi[i - 1];
        if MOD != 0 {
            small_phi[i] %= MOD;
        }
    }
    let mut Phi = FIArrayI64::new(x);

    for (ind, v) in FIArrayI64::keys(x).enumerate() {
        if v as usize <= y {
            Phi.arr[ind] = small_phi[v as usize];
            continue;
        }
        let vsqrt = v.isqrt();

        let mut phi_v = if MOD != 0 {
            (sum_n_i64::<MOD>(v) + MOD - v as i64 % MOD) % MOD
        } else {
            sum_n_i64::<MOD>(v) - v as i64
        };
        for i in 2..=vsqrt {
            let c = (small_phi[i as usize] - small_phi[i as usize - 1]) * (v / i) as i64;
            phi_v -= if MOD == 0 { c } else { c % MOD };
            phi_v -= Phi[v / i];
            if MOD != 0 {
                phi_v %= MOD;
            }
        }
        phi_v += Phi[vsqrt] * vsqrt as i64;
        if MOD != 0 {
            phi_v %= MOD;
            if phi_v < 0 {
                phi_v += MOD;
            }
        }
        Phi.arr[ind] = phi_v;
    }
    Phi
}
#[must_use]
pub fn mertens(x: usize) -> FIArrayI64 {
    let y = if x > 1023 {
        (1e9 as usize).min((x as f64).powf(2. / 3.) as usize >> 2)
    } else {
        x as usize
    };
    let mut small_m = mobius_sieve_i16(y + 1);
    for i in 2..=y {
        small_m[i] += small_m[i - 1];
    }
    let mut M = FIArrayI64::new(x);

    for (i, v) in FIArrayI64::keys(x).enumerate() {
        if v as usize <= y {
            M.arr[i] = i64::from(small_m[v as usize]);
            continue;
        }
        let vsqrt = v.isqrt();

        let mut mu_v = 1;
        mu_v -= v as i64;

        for i in 2..=vsqrt {
            mu_v -= (M.arr[i as usize - 1] - M.arr[i as usize - 2]) * (v / i) as i64;
            mu_v -= M[v / i];
        }
        mu_v += vsqrt as i64 * M[vsqrt];
        M.arr[i] = mu_v;
    }
    M
}
#[must_use]
pub fn divisor_summatory(x: usize) -> FIArrayI64 {
    /* let u = FIArrayI64::unit(x);
    dirichlet_mul_i64(&u, &u, x) */
    let mut zeta = DirichletFenwickI64::zeta(x);
    let lim = iroot::<8>(x) + 1;
    let mut primes = vec![];
    for p in 2..lim {
        if zeta.get_bucket_prefix(p - 1) == 1 {
            continue;
        }
        primes.push(p);
        zeta.sparse_mul_at_most_one(p, 1);
    }
    let zeta_lim = FIArrayI64::from(zeta);
    let mut zeta_2 = DirichletFenwickI64::from(mult_i64(&zeta_lim, &zeta_lim));
    for &p in primes.iter().rev() {
        zeta_2.sparse_mul_unlimited(p, 2);
    }
    let approx = FIArrayI64::from(zeta_2);
    mult_correction(&approx, &primes, |_, _, e| e as i64 + 1)
}

/// O(n^\frac12 \log n \log \log n) time, O(n^\frac12) space
#[must_use]
pub fn count_squarefree(x: usize) -> FIArray {
    let mut s = DirichletFenwick::zeta(x);
    for p in wheel_sieve(s.isqrt as u64) {
        let p = p as usize;
        s.sparse_mul_at_most_one(p * p, 1); // https://github.com/ecnerwala/cp-book/blob/master/src/dirichlet_series.hpp
    }
    s.into()
}
#[must_use]
pub fn sqf(x: usize) -> FIArrayI64 {
    let zeta = FIArrayI64::unit(x);
    let R2 = zeta.isqrt;
    let len = zeta.arr.len();

    let mob = mobius_sieve(R2 + 1);
    let mut mertens_sqrt = FIArrayI64::new(x);
    for i in 1..=R2 {
        mertens_sqrt[i * i] += mob[i] as i64;
    }
    for i in 1..len {
        mertens_sqrt.arr[i] += mertens_sqrt.arr[i - 1];
    }
    mult_sparse_i64(&zeta, &mertens_sqrt)
}
#[must_use]
pub fn sqf_icy(x: usize) -> FIArrayI64 {
    let mut Sqf = FIArrayI64::new(x);
    let mob = mobius_sieve(Sqf.isqrt + 1);
    for i in 1..=Sqf.isqrt {
        if mob[i] != 0 {
            Sqf.arr[i - 1] = 1;
        }
    }
    let s_N = Sqf.isqrt;
    let len = Sqf.arr.len();
    for d in 1..=s_N {
        if mob[d] == 0 {
            continue;
        }
        let M = x / (d * d);
        let mut l = 1;
        while l <= M && l <= s_N {
            let val = M / l;
            if val == 0 {
                break;
            }
            let r_bound = s_N.min(M / val);
            let term = mob[d] as i64 * val as i64;
            Sqf.arr[len - l] += term;
            if s_N < len - r_bound {
                Sqf.arr[len - r_bound - 1] -= term;
            } else {
                break;
            }

            l = r_bound + 1;
        }
    }
    for i in 1..s_N {
        Sqf.arr[i] += Sqf.arr[i - 1];
    }
    for i in 2..=len - s_N {
        Sqf.arr[len - i] += Sqf.arr[len - i + 1];
    }

    Sqf
}
/// Time complexity depends on sparsity:
/// On sparse inputs, i.e. only taking values at a \frac{1}{\log n} fraction of integers, takes O(n^2/3) time.
/// On dense inputs the function will take O(n^2/3 \log n) time.
#[must_use]
pub fn pseudo_euler_transform(a: FIArray) -> FIArray {
    let n = a.x;
    let rt = a.isqrt;
    let x = iroot::<3>(n) + 1;
    let mut r = a;
    let a = r.clone();
    for e in &mut r.arr[x - 1..] {
        *e -= a.arr[x - 2];
    }
    r.arr[..x - 1].fill(0);
    for i in x..=rt {
        let vi = a.arr[i - 1] - a.arr[i - 2];
        if vi == 0 {
            continue;
        }
        let n_i = n / i;
        for z in 1..=n_i / i {
            r[n / z] += vi * (a[n_i / z] - a.arr[i - 2]);
        }
    }
    let mut r = DirichletFenwick::from(r);
    r.bit.add(0, 1);
    for i in (2..x).rev() {
        let cur = a.arr[i - 1] - a.arr[i - 2];
        if cur == 0 {
            continue;
        }
        r.sparse_mul_unlimited(i, cur);
    }
    r.into()
}

#[must_use]
pub fn inverse_pseudo_euler_transform(a: FIArray) -> FIArray {
    let n = a.x;
    let len = a.arr.len();
    let mut r = FIArray::new(n);
    let mut a_bit = DirichletFenwick::from(a);
    let x = iroot::<3>(n) + 1;
    for i in 2..x {
        let cur = a_bit.bit.sum(i - 1) - 1;
        if cur == 0 {
            continue;
        }
        r.arr[i - 1] = cur;
        a_bit.sparse_mul_at_most_one(i, cur);
    }
    let a = FIArray::from(a_bit);
    for i in x..=len {
        r.arr[i - 1] = a.arr[i - 1] - a.arr[i - 2];
    }
    for i in x..=a.isqrt {
        let vi = r.arr[i - 1];
        if vi == 0 {
            continue;
        }
        let n_i = n / i;
        for z in 1..=n_i / i {
            let v = vi * (a[n_i / z] - a.arr[i - 2]);
            r.arr[len - z] -= v;
            if z > 1 {
                r.arr[len - z + 1] += v;
            }
        }
    }
    for i in 1..len {
        r.arr[i] += r.arr[i - 1];
    }
    r
}

// faster by a log factor, but much more susceptible to overflow - multiplies input by a factor of 1296 before reducing
#[must_use]
pub fn pseudo_euler_transform_fraction(a: FIArray) -> FIArray {
    const INVS: [usize; 4] = [0, 6, 3, 2];

    let mut a = a;
    let len = a.arr.len();
    let n = a.x;
    let rt = a.isqrt;
    let x = iroot::<4>(n) + 1;

    for i in (1..len).rev() {
        a.arr[i] -= a.arr[i - 1];
        a.arr[i] *= INVS[1];
    }
    a.arr[0] *= INVS[1]; // kinda pointless tbh
    for i in (x..=rt).rev() {
        let v = a.arr[i - 1] / INVS[1];
        if v == 0 {
            continue;
        }
        let mut e = 1;
        let mut pi = i;
        let mut pv = v;
        while pi <= n / i {
            e += 1;
            pi *= i;
            pv *= v;
            a[pi] += pv * INVS[e];
        }
    }

    let mut v = FIArray::new(n);
    for i in x..=len {
        v.arr[i - 1] = v.arr[i - 2] + a.arr[i - 1];
    }

    let mut ret = v.clone();
    for e in &mut ret.arr {
        *e = (*e + INVS[1]) * 6 * INVS[1].pow(2);
    }

    let mut v_2 = mult(&v, &v);
    for i in x..=len {
        ret.arr[i - 1] += v_2.arr[i - 1] * 3 * INVS[1];
    }
    {
        //v_2 = mult_sparse(&v, &v_2);
        v.arr[rt..].fill(0);
        for i in x..=rt {
            let y = v.arr[i - 1] - v.arr[i - 2];
            if y != 0 {
                for j in 1..=rt / i {
                    v.arr[len - j] += y * v_2.arr[len - i * j];
                }
            }
        }
        v.arr[..rt].fill(0);
        core::mem::swap(&mut v_2, &mut v);
    }

    for i in 1..=len {
        ret.arr[i - 1] += v_2.arr[i - 1];
        ret.arr[i - 1] /= const { 6 * INVS[1].pow(3) };
    }
    let mut ret = DirichletFenwick::from(ret);
    for i in (2..x).rev() {
        let ai = a.arr[i - 1] / INVS[1];
        if ai == 0 {
            continue;
        }
        ret.sparse_mul_unlimited(i, ai);
    }
    ret.into()
}
// faster by a log factor, but more susceptible to overflow - multiplies input by a factor of 6 before reducing
#[must_use]
pub fn inverse_pseudo_euler_transform_fraction(a: FIArray) -> FIArray {
    const INVS: [usize; 4] = [0, 6, 3, 2];
    let mut a = DirichletFenwick::from(a);
    let rt = a.isqrt;
    let n = a.x;
    let len = a.bit.0.len();

    let mut ret = FIArray::new(n);

    let x = iroot::<4>(n) + 1;
    for i in 2..x {
        let vi = a.bit.sum(i - 1) - 1;
        if vi == 0 {
            continue;
        }
        ret.arr[i - 1] = vi;
        a.sparse_mul_at_most_one(i, vi);
    }
    a.bit.dec(0);
    let mut a = FIArray::from(a);

    // a now equals a_t - 1
    // compute log(a_t) using log(x + 1) = x^3 / 3 - x^2 / 2 + x
    // in order to not have to deal with rational numbers, we compute 6 * log(a_t)
    // and adjust later

    for i in x..=len {
        ret.arr[i - 1] = a.arr[i - 1] * INVS[1];
    }

    let mut a_2 = crate::p300_399::e362::mult(&a, &a);
    /* let ind = pow_zeta.get_index(x.pow(2));
    assert!(pow_zeta.arr[..ind].iter().all(|&e| e == 0));
    dbg!(ind, x.pow(2), rt - 1, len, len - ind); */

    for i in rt + 1..=len {
        ret.arr[i - 1] -= a_2.arr[i - 1] * INVS[2];
    }
    {
        //pow_zeta = mult_sparse(&zeta, &pow_zeta);

        a.arr[rt..].fill(0);
        for i in x..=rt {
            let y = a.arr[i - 1] - a.arr[i - 2];
            if y != 0 {
                for j in 1..=rt / i {
                    a.arr[len - j] += y * a_2.arr[len - i * j];
                }
            }
        }
        //zeta.arr[..rt].fill(0);
        core::mem::swap(&mut a_2, &mut a);
    }
    let ind = a_2.get_index(x.pow(3));
    //dbg!(ind, len, len - ind);
    //assert!(pow_zeta.arr[..ind].iter().all(|&e| e == 0)); */
    for i in ind + 1..=len {
        ret.arr[i - 1] += a_2.arr[i - 1] * INVS[3];
    }

    // correction phase: get rid of contributions of prime powers
    for i in (x..len).rev() {
        ret.arr[i] -= ret.arr[i - 1];
    }

    for x in x..=rt {
        let v = ret.arr[x - 1] / 6;
        if v == 0 {
            continue;
        }
        let mut e = 1;
        let mut pv = v;
        let mut px = x;
        while px <= n / x {
            e += 1;
            px *= x;
            pv *= v;

            ret[px] -= pv * INVS[e];
        }
    }
    for i in 1..x - 1 {
        ret.arr[i] += ret.arr[i - 1];
    }
    for i in x - 1..len {
        ret.arr[i] /= 6;
        ret.arr[i] += ret.arr[i - 1];
    }
    ret
}

/* #[must_use]
pub fn sqf(x: usize) -> FIArray {
    let mut Sqf = FIArray::eps(x);
    let xsqrt = Sqf.isqrt;
    let mut d2 = 1;
    for d in 2..=xsqrt.isqrt() {
        d2 += (d << 1) - 1;
        if Sqf.arr[d2 - 1] == 0 {
            continue;
        }
        for m in (d2..=xsqrt).step_by(d2) {
            Sqf.arr[m - 1] = 0;
        }
    }
    let mut sqrts = FIArray::unit(x);
    for v in &mut sqrts.arr {
        *v = v.isqrt();
    }
    for (i, v) in FIArray::keys(x).enumerate().skip(1) {
        if v <= xsqrt {
            Sqf.arr[i] += Sqf.arr[i - 1];
            continue;
        }
        let b = iroot::<3>(v);
        let a = v / (b * b);
        let mut sqf = v + Sqf.arr[a - 1] * b - sqrts[v]; // v.isqrt();
        for i in 2..=a {
            if Sqf.arr[i - 1] != Sqf.arr[i - 2] {
                sqf -= sqrts[v / i]; //(v / i).isqrt();
            }
        }
        for i in 2..=b {
            sqf -= Sqf[v / (i * i)];
        }
        Sqf.arr[i] = sqf;
    }
    Sqf
} */
#[must_use]
pub fn totient_sum_single<const MOD: i64>(x: usize) -> i64 {
    let M = mertens(x);
    let mut sum = M[x] + sum_n_i64::<MOD>(x);
    if MOD != 0 {
        sum %= MOD;
    }
    let isqrt = x.isqrt();
    for i in 2..=isqrt {
        sum += i as i64 * M[x / i];
        sum += (M.arr[i as usize - 1] - M.arr[i as usize - 2]) * sum_n_i64::<MOD>(x / i);
        if MOD != 0 {
            sum %= MOD;
        }
    }
    sum -= sum_n_i64::<MOD>(isqrt) * M[isqrt];
    if MOD != 0 {
        sum %= MOD;
        if sum < 0 {
            sum += MOD;
        }
    }
    sum
}
#[must_use]
pub fn mertens_slow(x: usize) -> FIArrayI64 {
    let y = if x > 1023 {
        x.isqrt() as usize
    } else {
        x as usize
    };
    let mobius = mobius_sieve(y + 1);
    let mut M = FIArrayI64::new(x);

    M.arr[0] = 1;
    for (i, v) in FIArrayI64::keys(x).enumerate().skip(1) {
        if v as usize <= y {
            M.arr[i] = i64::from(mobius[i + 1]) + M.arr[i - 1];
            continue;
        }
        let vsqrt = v.isqrt();

        let mut mu_v = 1;
        mu_v -= v as i64;

        for i in 2..=vsqrt {
            mu_v -= i64::from(mobius[i as usize]) * (v / i) as i64;
            mu_v -= M[v / i];
        }
        mu_v += vsqrt as i64 * M[vsqrt];
        M.arr[i] = mu_v;
    }
    M
}

pub fn sum_over_primes<const MOD: i64>(
    x: usize,
    mut f: impl FnMut(usize) -> i64,
    mut F: impl FnMut(usize) -> i64,
) -> FIArrayI64 {
    let primes = wheel_sieve(x.isqrt() as u64 + 1);
    let mut s = FIArrayI64::new(x);
    let keys = FIArrayI64::keys(x).collect_vec();

    unsafe { core::hint::assert_unchecked(s.arr.len() == keys.len()) };
    for (i, &v) in keys.iter().enumerate() {
        s.arr[i] = F(v);
    }
    unsafe { core::hint::assert_unchecked(s.arr.len() > x.isqrt()) };
    for p in primes {
        let p = p as usize;
        let sp = s.arr[p - 2];
        let fp = f(p);
        for (i, &v) in keys.iter().enumerate().rev() {
            if v < p * p {
                break;
            }
            s.arr[i] -= fp * (s[v / p] - sp);
            if MOD != 0 {
                s.arr[i] %= MOD;
                if s.arr[i] < 0 {
                    s.arr[i] += MOD;
                }
            }
        }
    }
    s
}

pub fn trick<const MOD: i64>(
    n: i64,
    mut g: impl FnMut(i64) -> i64,
    mut F: impl FnMut(i64) -> i64,
) -> i64 {
    let mut res = 0;
    let mut m = n;
    while m > 0 {
        res += ((F(n / m) + MOD - F(n / (m + 1))) % MOD * g(m)) % MOD;
        if res >= MOD {
            res -= MOD;
        }
        m = n / (n / m + 1);
    }
    res
}

use paste::paste;
macro_rules! min25_sieve_impl_for {
    ($($type:ty),+) => { $(
        paste!{
            #[must_use] pub const fn [<divisor_summatory_ $type>](x: usize) -> $type {
                let mut sum = 0;
                let sqrt = x.isqrt();
                let mut n = 1;
                while n <= sqrt {
                    sum += (x / n) as $type;
                    n += 1;
                }
                sum <<= 1;
                sum - (sqrt * sqrt) as $type
            }
            pub const fn [<sum_n_ $type>]<const MOD: $type>(x: usize) -> $type {
                let x = x as $type;
                if MOD == 0 {
                    if x & 1 == 0 {
                        (x / 2) * (x + 1)
                    } else {
                        ((x + 1) / 2) * x
                    }
                } else {
                    let x = x % (MOD << 1);
                    (if x & 1 == 0 {
                        (x / 2) * (x + 1)
                    } else {
                        ((x + 1) / 2) * x
                    }) % MOD
                }
            }
            // note: does not require the functions f and g to be multiplicative
            pub fn [<dirichlet_mul_ $type>](F: &[<FIArray $type:camel>], G: &[<FIArray $type:camel>], n: usize) -> [<FIArray $type:camel>] {
                let mut H = [<FIArray $type:camel>]::new(n as _);
                let len = H.arr.len();

                let rt_n = n.isqrt();

                let to_ord = |x| {
                    if x <= rt_n { x } else { len + 1 - (n / x) }
                };
                let mut propogate = |(x0, x1), (y0, y1), (z0, z1)| {
                    let f_x1 = F.arr[x1 - 1];
                    let g_y1 = G.arr[y1 - 1];
                    let f_x0_1 = F.arr.get(x0 - 2).copied().unwrap_or_default();
                    let g_y0_1 = G.arr.get(y0 - 2).copied().unwrap_or_default();

                    let t = (f_x1 - f_x0_1) * (g_y1 - g_y0_1);
                    H.arr[z0 - 1] += t;
                    if let Some(v) = H.arr.get_mut(z1) {
                        *v -= t;
                    }
                };
                propogate((1, 1), (1, 1), (1, len));
                for k in 2..=len {
                    if k > rt_n {
                        let z = len + 1 - k;
                        for x in 2.. {
                            let y_lo_ord = 1 + to_ord(x).max(to_ord(z));
                            let y_hi_ord = to_ord(n / (x * z));
                            if y_hi_ord < y_lo_ord {
                                break;
                            }
                            propogate((x, x), (y_lo_ord, y_hi_ord), (k, k));
                            propogate((y_lo_ord, y_hi_ord), (x, x), (k, k));
                        }
                    }
                    propogate((1, 1), (k, k), (k, len));
                    propogate((k, k), (1, 1), (k, len));
                    let x = k;
                    if x <= rt_n {
                        for y in 2..k {
                            let z_lo_ord = to_ord(x * y);
                            let z_hi_ord = to_ord(n / x);
                            if z_hi_ord < z_lo_ord {
                                break;
                            }
                            propogate((x, x), (y, y), (z_lo_ord, z_hi_ord));
                            propogate((y, y), (x, x), (z_lo_ord, z_hi_ord));
                        }

                        propogate((x, x), (x, x), (to_ord(x * x), len));
                    }
                }

                for i in 1..len {
                    H.arr[i] += H.arr[i - 1];
                }
                H
            }

            pub fn [<dirichlet_mul_single_ $type>](F: &[<FIArray $type:camel>], G: &[<FIArray $type:camel>]) -> $type {
                let rt_n = F.isqrt;
                let len = F.arr.len();

                let mut ret = F.arr[0] * G.arr[len - 1] + G.arr[0] * F.arr[len - 1] - F.arr[rt_n - 1] * G.arr[rt_n - 1];
                for i in 2..=rt_n {
                    ret += (F.arr[i - 1] - F.arr[i - 2]) * G.arr[len - i]
                        + (G.arr[i - 1] - G.arr[i - 2]) * F.arr[len - i];
                }
                ret
            }

            pub fn [<dirichlet_div_ $type>](H: &[<FIArray $type:camel>], G: &[<FIArray $type:camel>], n: usize) -> [<FIArray $type:camel>] {
                let mut F = [<FIArray $type:camel>]::new(n as _);
                let len = F.arr.len();
                let mut modified = BitArray::zeroed(len + 2);

                let rt_n = n.isqrt();

                let mut H = H.clone();
                for i in (1..len).rev() {
                    H.arr[i] -= H.arr[i - 1];
                }
                let to_ord = |x| {
                    if x <= rt_n { x } else { len + 1 - (n / x) }
                };
                let mut propogate = |(x0, x1), (y0, y1), (z0, z1)| {
                    let g_y1 = G.arr[y1 - 1];
                    let f_x0_1 = F.arr.get(x0 - 2).copied().unwrap_or_default();
                    let g_y0_1 = G.arr.get(y0 - 2).copied().unwrap_or_default();
                    if !modified.get(x1) {
                        F.arr[x1 - 1] = f_x0_1 + H.arr[z0-1] / (g_y1 - g_y0_1);
                        modified.set(x1);

                    }
                    let f_x1 = F.arr[x1 - 1];
                    let t = (f_x1 - f_x0_1) * (g_y1 - g_y0_1);
                    H.arr[z0 - 1] -= t;
                    if let Some(v) = H.arr.get_mut(z1) {
                        *v += t;
                    }
                };
                propogate((1, 1), (1, 1), (1, len));
                for k in 2..=len {
                    let z = len + 1 - k;
                    for x in 2.. {
                        let y_lo_ord = 1 + to_ord(x).max(to_ord(z));
                        let y_hi_ord = to_ord(n / (x * z));
                        if y_hi_ord < y_lo_ord {
                            break;
                        }
                        propogate((x, x), (y_lo_ord, y_hi_ord), (k, k));
                        propogate((y_lo_ord, y_hi_ord), (x, x), (k, k));
                    }
                    propogate((1, 1), (k, k), (k, len));
                    propogate((k, k), (1, 1), (k, len));
                    let x = k;
                    for y in 2..k {
                        let z_lo_ord = to_ord(x * y);
                        let z_hi_ord = to_ord(n / x);
                        if z_hi_ord < z_lo_ord {
                            break;
                        }
                        propogate((x, x), (y, y), (z_lo_ord, z_hi_ord));
                        propogate((y, y), (x, x), (z_lo_ord, z_hi_ord));
                    }

                    if x <= rt_n {
                        propogate((x, x), (x, x), (to_ord(x * x), len));
                    }
                }

                F
            }

            pub fn [<dirichlet_inv_ $type>]( G: &[<FIArray $type:camel>], n: usize) -> [<FIArray $type:camel>] {
                let mut F = [<FIArray $type:camel>]::new(n as _);
                let len = F.arr.len();
                let mut modified = BitArray::zeroed(len + 2);

                let rt_n = n.isqrt();

                let mut H = [<FIArray $type:camel>]::eps(n as _);
                for i in (1..len).rev() {
                    H.arr[i] -= H.arr[i - 1];
                }
                let to_ord = |x| {
                    if x <= rt_n { x } else { len + 1 - (n / x) }
                };
                let mut propogate = |(x0, x1), (y0, y1), (z0, z1)| {
                    let g_y1 = G.arr[y1 - 1];
                    let f_x0_1 = F.arr.get(x0 - 2).copied().unwrap_or_default();
                    let g_y0_1 = G.arr.get(y0 - 2).copied().unwrap_or_default();
                    if !modified.get(x1) {
                        F.arr[x1 - 1] = f_x0_1 + H.arr[z0-1] / (g_y1 - g_y0_1);
                        modified.set(x1);

                    }
                    let f_x1 = F.arr[x1 - 1];
                    let t = (f_x1 - f_x0_1) * (g_y1 - g_y0_1);
                    H.arr[z0 - 1] -= t;
                    if let Some(v) = H.arr.get_mut(z1) {
                        *v += t;
                    }
                };
                propogate((1, 1), (1, 1), (1, len));
                for k in 2..=len {
                    let z = len + 1 - k;
                    for x in 2.. {
                        let y_lo_ord = 1 + to_ord(x).max(to_ord(z));
                        let y_hi_ord = to_ord(n / (x * z));
                        if y_hi_ord < y_lo_ord {
                            break;
                        }
                        propogate((x, x), (y_lo_ord, y_hi_ord), (k, k));
                        propogate((y_lo_ord, y_hi_ord), (x, x), (k, k));
                    }
                    propogate((1, 1), (k, k), (k, len));
                    propogate((k, k), (1, 1), (k, len));
                    let x = k;
                    for y in 2..k {
                        let z_lo_ord = to_ord(x * y);
                        let z_hi_ord = to_ord(n / x);
                        if z_hi_ord < z_lo_ord {
                            break;
                        }
                        propogate((x, x), (y, y), (z_lo_ord, z_hi_ord));
                        propogate((y, y), (x, x), (z_lo_ord, z_hi_ord));
                    }

                    if x <= rt_n {
                        propogate((x, x), (x, x), (to_ord(x * x), len));
                    }
                }

                F
            }


            pub fn [<dirichlet_mul_with_buffer_ $type>](F: &[<FIArray $type:camel>], G: &[<FIArray $type:camel>], n: usize, H: &mut [<FIArray $type:camel>]) {
                H.arr.fill(0);
                let len = H.arr.len();

                let rt_n = n.isqrt();

                let to_ord = |x| {
                    if x <= rt_n { x } else { len + 1 - (n / x) }
                };
                let mut propogate = |(x0, x1), (y0, y1), (z0, z1)| {
                    let f_x1 = F.arr[x1 - 1];
                    let g_y1 = G.arr[y1 - 1];
                    let f_x0_1 = F.arr.get(x0 - 2).copied().unwrap_or_default();
                    let g_y0_1 = G.arr.get(y0 - 2).copied().unwrap_or_default();

                    let t = (f_x1 - f_x0_1) * (g_y1 - g_y0_1);
                    H.arr[z0 - 1] += t;
                    if let Some(v) = H.arr.get_mut(z1) {
                        *v -= t;
                    }
                };
                propogate((1, 1), (1, 1), (1, len));
                for k in 2..=len {
                    let z = len + 1 - k;
                    for x in 2.. {
                        let y_lo_ord = 1 + to_ord(x).max(to_ord(z));
                        let y_hi_ord = to_ord(n / (x * z));
                        if y_hi_ord < y_lo_ord {
                            break;
                        }
                        propogate((x, x), (y_lo_ord, y_hi_ord), (k, k));
                        propogate((y_lo_ord, y_hi_ord), (x, x), (k, k));
                    }
                    propogate((1, 1), (k, k), (k, len));
                    propogate((k, k), (1, 1), (k, len));
                    let x = k;
                    for y in 2..k {
                        let z_lo_ord = to_ord(x * y);
                        let z_hi_ord = to_ord(n / x);
                        if z_hi_ord < z_lo_ord {
                            break;
                        }
                        propogate((x, x), (y, y), (z_lo_ord, z_hi_ord));
                        propogate((y, y), (x, x), (z_lo_ord, z_hi_ord));
                    }

                    if x <= rt_n {
                        propogate((x, x), (x, x), (to_ord(x * x), len));
                    }
                }

                for i in 1..len {
                    H.arr[i] += H.arr[i - 1];
                }
            }

            // note: does not require the functions f and g to be multiplicative
            pub fn [<dirichlet_mulmod_ $type>]<const MOD: $type>(F: &[<FIArray $type:camel>], G: &[<FIArray $type:camel>], n: usize) -> [<FIArray $type:camel>] {
                let mut H = [<FIArray $type:camel>]::new(n as _);
                let len = H.arr.len();

                let rt_n = n.isqrt();

                let to_ord = |x| {
                    if x <= rt_n { x } else { len + 1 - (n / x) }
                };
                let mut propogate = |(x0, x1), (y0, y1), (z0, z1)| {
                    let f_x1 = F.arr[x1 - 1];
                    let g_y1 = G.arr[y1 - 1];
                    let f_x0_1 = F.arr.get(x0 - 2).copied().unwrap_or_default();
                    let g_y0_1 = G.arr.get(y0 - 2).copied().unwrap_or_default();

                    let t = ((f_x1 - f_x0_1) % MOD * (g_y1 - g_y0_1) % MOD) % MOD;
                    H.arr[z0 - 1] += t;
                    //H.arr[z0 - 1] %= MOD;
                    if let Some(v) = H.arr.get_mut(z1) {
                        *v -= t;
                        //*v %= MOD;
                    }
                };
                propogate((1, 1), (1, 1), (1, len));
                for k in 2..=len {
                    let z = len + 1 - k;
                    for x in 2.. {
                        let y_lo_ord = 1 + to_ord(x).max(to_ord(z));
                        let y_hi_ord = to_ord(n / (x * z));
                        if y_hi_ord < y_lo_ord {
                            break;
                        }
                        propogate((x, x), (y_lo_ord, y_hi_ord), (k, k));
                        propogate((y_lo_ord, y_hi_ord), (x, x), (k, k));
                    }
                    propogate((1, 1), (k, k), (k, len));
                    propogate((k, k), (1, 1), (k, len));
                    let x = k;
                    for y in 2..k {
                        let z_lo_ord = to_ord(x * y);
                        let z_hi_ord = to_ord(n / x);
                        if z_hi_ord < z_lo_ord {
                            break;
                        }
                        propogate((x, x), (y, y), (z_lo_ord, z_hi_ord));
                        propogate((y, y), (x, x), (z_lo_ord, z_hi_ord));
                    }

                    if x <= rt_n {
                        propogate((x, x), (x, x), (to_ord(x * x), len));
                    }
                }

                for i in 1..len {
                    H.arr[i] %= MOD;
                    H.arr[i] += H.arr[i - 1];
                    H.arr[i] %= MOD;
                    if H.arr[i] < 0 {
                        H.arr[i] += MOD;
                    }
                }
                H
            }

            pub fn [<dirichlet_mulmod_with_buffer_ $type>]<const MOD: $type>(
                F: &[<FIArray $type:camel>],
                G: &[<FIArray $type:camel>],
                n: usize,
                H: &mut [<FIArray $type:camel>],
            ) {
                H.arr.fill(0);
                let len = H.arr.len();

                let rt_n = n.isqrt();

                let to_ord = |x| {
                    if x <= rt_n { x } else { len + 1 - (n / x) }
                };
                let mut propogate = |(x0, x1), (y0, y1), (z0, z1)| {
                    let f_x1 = F.arr[x1 - 1];
                    let g_y1 = G.arr[y1 - 1];
                    let f_x0_1 = F.arr.get(x0 - 2).copied().unwrap_or_default();
                    let g_y0_1 = G.arr.get(y0 - 2).copied().unwrap_or_default();

                    let t = ((f_x1 - f_x0_1) % MOD * (g_y1 - g_y0_1) % MOD) % MOD;
                    H.arr[z0 - 1] += t;
                    //H.arr[z0 - 1] %= MOD;
                    if let Some(v) = H.arr.get_mut(z1) {
                        *v -= t;
                        //*v %= MOD;
                    }
                };
                propogate((1, 1), (1, 1), (1, len));
                for k in 2..=len {
                    let z = len + 1 - k;
                    for x in 2.. {
                        let y_lo_ord = 1 + to_ord(x).max(to_ord(z));
                        let y_hi_ord = to_ord(n / (x * z));
                        if y_hi_ord < y_lo_ord {
                            break;
                        }
                        propogate((x, x), (y_lo_ord, y_hi_ord), (k, k));
                        propogate((y_lo_ord, y_hi_ord), (x, x), (k, k));
                    }
                    propogate((1, 1), (k, k), (k, len));
                    propogate((k, k), (1, 1), (k, len));
                    let x = k;
                    for y in 2..k {
                        let z_lo_ord = to_ord(x * y);
                        let z_hi_ord = to_ord(n / x);
                        if z_hi_ord < z_lo_ord {
                            break;
                        }
                        propogate((x, x), (y, y), (z_lo_ord, z_hi_ord));
                        propogate((y, y), (x, x), (z_lo_ord, z_hi_ord));
                    }

                    if x <= rt_n {
                        propogate((x, x), (x, x), (to_ord(x * x), len));
                    }
                }

                for i in 1..len {
                    H.arr[i] %= MOD;
                    H.arr[i] += H.arr[i - 1];
                    H.arr[i] %= MOD;
                    if H.arr[i] < 0 {
                        H.arr[i] += MOD;
                    }
                }
            }


        }
    )+ };
}
min25_sieve_impl_for!(u32, i32, u64, i64, usize, isize, u128, i128);

// O(log(k)x^(2/3)) time, O(x^(1/2)) space, specifically ~6x^(1/2) i64's
#[must_use]
pub fn general_divisor_summatory<const MOD: i64>(x: usize, mut k: usize) -> FIArrayI64 {
    assert!(k > 1);
    let mut buffer = FIArrayI64::new(x);
    let mut u = FIArrayI64::unit(x);
    let mut r = if k & 1 == 1 {
        k >>= 1;
        dirichlet_mulmod_with_buffer_i64::<MOD>(&u, &u, x as _, &mut buffer);
        core::mem::swap(&mut u, &mut buffer);
        FIArrayI64::unit(x)
    } else {
        FIArrayI64::eps(x)
    };
    while k > 1 {
        if k & 1 == 1 {
            dirichlet_mulmod_with_buffer_i64::<MOD>(&r, &u, x as _, &mut buffer);
            core::mem::swap(&mut r, &mut buffer);
            //r = r * u;
        }
        //u = u * u;
        dirichlet_mulmod_with_buffer_i64::<MOD>(&u, &u, x as _, &mut buffer);
        core::mem::swap(&mut u, &mut buffer);
        k >>= 1;
    }
    //(r * u)
    dirichlet_mulmod_with_buffer_i64::<MOD>(&r, &u, x as _, &mut buffer);
    buffer
}

#[must_use]
pub fn general_divisor_summatory_alt<const MOD: i64>(x: usize, mut k: usize) -> FIArrayI64 {
    assert!(k > 1);
    let mut buffer = FIArrayI64::new(x);
    let mut u = FIArrayI64::unit(x);
    let mut r = {
        while k > 1 && k & 1 == 0 {
            dirichlet_mulmod_with_buffer_i64::<MOD>(&u, &u, x as _, &mut buffer);
            core::mem::swap(&mut u, &mut buffer);
            k >>= 1;
        }
        if k == 1 {
            return u;
        }
        let ret = u.clone();
        dirichlet_mulmod_with_buffer_i64::<MOD>(&u, &u, x as _, &mut buffer);
        core::mem::swap(&mut u, &mut buffer);
        k >>= 1;
        ret
    };
    while k > 1 {
        if k & 1 == 1 {
            dirichlet_mulmod_with_buffer_i64::<MOD>(&r, &u, x as _, &mut buffer);
            core::mem::swap(&mut r, &mut buffer);
            //r = r * u;
        }
        //u = u * u;
        dirichlet_mulmod_with_buffer_i64::<MOD>(&u, &u, x as _, &mut buffer);
        core::mem::swap(&mut u, &mut buffer);
        k >>= 1;
    }
    //(r * u)
    dirichlet_mulmod_with_buffer_i64::<MOD>(&r, &u, x as _, &mut buffer);
    buffer
}
