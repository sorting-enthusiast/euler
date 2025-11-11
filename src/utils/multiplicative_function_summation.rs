use itertools::Itertools;

use crate::utils::{
    FIArray::{
        FIArrayI32, FIArrayI64, FIArrayI128, FIArrayIsize, FIArrayU32, FIArrayU64, FIArrayU128,
        FIArrayUsize,
    },
    bit_array::BitArray,
    prime_sieves::sift,
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
pub fn totient_sum<const MOD: i64>(x: i64) -> FIArrayI64 {
    let y = if x > 1023 {
        (1e8 as usize).min((x as f64).powf(2. / 3.) as usize >> 1)
    } else {
        x as usize
    };
    let mut small_phi = totient_sieve(y + 1);
    for i in 2..=y {
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

        let mut phi_v = (sum_n_i64::<MOD>(v) + MOD - v % MOD) % MOD;
        for i in 2..=vsqrt {
            phi_v -= ((small_phi[i as usize] - small_phi[i as usize - 1]) * (v / i)) % MOD;
            phi_v -= Phi[v / i];
            phi_v %= MOD;
        }
        phi_v += Phi[vsqrt] * vsqrt;
        phi_v %= MOD;
        if phi_v < 0 {
            phi_v += MOD;
        }
        Phi.arr[ind] = phi_v;
    }
    Phi
}
#[must_use]
pub fn mertens(x: i64) -> FIArrayI64 {
    let y = if x > 1023 {
        (1e8 as usize).min((x as f64).powf(2. / 3.) as usize >> 2)
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
        mu_v -= v;

        for i in 2..=vsqrt {
            mu_v -= i64::from(small_m[i as usize] - small_m[i as usize - 1]) * (v / i);
            mu_v -= M[v / i];
        }
        mu_v += vsqrt * M[vsqrt];
        M.arr[i] = mu_v;
    }
    M
}
#[must_use]
pub fn divisor_summatory(x: i64) -> FIArrayI64 {
    let u = FIArrayI64::unit(x);
    dirichlet_mul_i64(&u, &u, x as _)
}

#[must_use]
pub fn totient_sum_single<const MOD: i64>(x: i64) -> i64 {
    let M = mertens(x);
    let mut sum = M[x] + sum_n_i64::<MOD>(x);
    sum %= MOD;
    let isqrt = x.isqrt();
    for i in 2..=isqrt {
        sum += i * M[x / i];
        sum += (M.arr[i as usize - 1] - M.arr[i as usize - 2]) * sum_n_i64::<MOD>(x / i);
        sum %= MOD;
    }
    sum -= sum_n_i64::<MOD>(isqrt) * M[isqrt];
    sum %= MOD;
    if sum < 0 {
        sum += MOD;
    }
    sum
}
#[must_use]
pub fn mertens_slow(x: i64) -> FIArrayI64 {
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
        mu_v -= v;

        for i in 2..=vsqrt {
            mu_v -= i64::from(mobius[i as usize]) * (v / i);
            mu_v -= M[v / i];
        }
        mu_v += vsqrt * M[vsqrt];
        M.arr[i] = mu_v;
    }
    M
}

pub fn sum_over_primes<const MOD: i64>(
    x: i64,
    mut f: impl FnMut(i64) -> i64,
    mut F: impl FnMut(i64) -> i64,
) -> FIArrayI64 {
    let primes = sift(x.isqrt() as u64 + 1);
    let mut s = FIArrayI64::new(x);
    let keys = FIArrayI64::keys(x).collect_vec();

    unsafe { core::hint::assert_unchecked(s.arr.len() == keys.len()) };
    for (i, &v) in keys.iter().enumerate() {
        s.arr[i] = F(v);
    }
    unsafe { core::hint::assert_unchecked(s.arr.len() as i64 > x.isqrt()) };
    for p in primes {
        let p = p as i64;
        let sp = s.arr[p as usize - 2];
        for (i, &v) in keys.iter().enumerate().rev() {
            if v < p * p {
                break;
            }
            s.arr[i] -= f(p) * (s[v / p] - sp);
            s.arr[i] %= MOD;
            if s.arr[i] < 0 {
                s.arr[i] += MOD;
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
            #[must_use] pub const fn [<divisor_summatory_ $type>](x: $type) -> $type {
                let mut sum = 0;
                let sqrt = x.isqrt();
                let mut n = 1;
                while n <= sqrt {
                    sum += x / n;
                    n += 1;
                }
                sum <<= 1;
                sum - sqrt * sqrt
            }
            pub const fn [<sum_n_ $type>]<const MOD: $type>(x: $type) -> $type {
                let x = x % (MOD << 1);
                (if x & 1 == 0 {
                    (x / 2) * (x + 1)
                } else {
                    ((x + 1) / 2) * x
                }) % MOD
            }
            pub fn [<min25_sieve_ $type>]<const MOD: $type>(
                x: $type,
                mut g: impl FnMut($type) -> $type,
                mut G: impl FnMut($type) -> $type,
                mut f: impl FnMut($type, $type) -> $type,
            ) -> [<FIArray $type:camel>] {
                let primes = super::prime_sieves::sift(x.isqrt() as u64 + 1);
                let mut s = [<FIArray $type:camel>]::new(x);
                let keys = [<FIArray $type:camel>]::keys(x).collect_vec();

                unsafe { core::hint::assert_unchecked(s.arr.len() == keys.len()) };
                for (i, &v) in keys.iter().enumerate() {
                    s.arr[i] = G(v);
                }
                for &p in &primes {
                    let sp = s.arr[p as usize - 2];
                    let p = p as $type;
                    for (i, &v) in keys.iter().enumerate().rev() {
                        if v < p * p {
                            break;
                        }
                        s.arr[i] -= g(p) * (s[v / p] - sp);
                        s.arr[i] %= MOD;
                        if s.arr[i] < 0 {
                            s.arr[i] += MOD;
                        }
                    }
                }

                for &p in primes.iter().rev() {
                    let sp = s.arr[p as usize - 1];
                    let p = p as $type;
                    for (i, &v) in keys.iter().enumerate().rev() {
                        if v < p * p {
                            break;
                        }
                        let mut e = 1;
                        let mut u = v / p;
                        while u >= p {
                            s.arr[i] += f(p, e) * (s[u] - sp) + f(p, e + 1);
                            s.arr[i] %= MOD;
                            if s.arr[i] < 0 {
                                s.arr[i] += MOD;
                            }
                            e += 1;
                            u /= p;
                        }
                    }
                }
                s
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
                H
            }

            pub fn [<dirichlet_mul_single_ $type>](F: &[<FIArray $type:camel>], G: &[<FIArray $type:camel>], n: usize) -> $type {
                let rt_n = n.isqrt();
                let mut ret = F.arr[0] * G[n as $type] + G.arr[0] * F[n as $type] - F[rt_n as $type] * G[rt_n as $type];
                for i in 2..=rt_n {
                    let ni = n / i;
                    ret += (F.arr[i - 1] - F.arr[i - 2]) * G[ni as $type]
                        + (G.arr[i - 1] - G.arr[i - 2]) * F[ni as $type];
                }
                ret
            }

            pub fn [<dirichlet_mul_single_zero_prefix_ $type>](F: &[<FIArray $type:camel>], G: &[<FIArray $type:camel>], n: usize, prefix_f: usize, prefix_g: usize) -> $type {
                let rt_n = n.isqrt();
                let mut ret = F.arr[0] * G[n as $type] + G.arr[0] * F[n as $type] - F[rt_n as $type] * G[rt_n as $type];
                for i in prefix_f..prefix_g {
                    let ni = n / i;
                    ret += (F.arr[i - 1] - F.arr[i - 2]) * G[ni as $type];
                }
                for i in prefix_g..=rt_n {
                    let ni = n / i;
                    ret += (F.arr[i - 1] - F.arr[i - 2]) * G[ni as $type]
                        + (G.arr[i - 1] - G.arr[i - 2]) * F[ni as $type];
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
pub fn general_divisor_summatory<const MOD: i64>(x: i64, mut k: u8) -> FIArrayI64 {
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
pub fn general_divisor_summatory_alt<const MOD: i64>(x: i64, mut k: u8) -> FIArrayI64 {
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
