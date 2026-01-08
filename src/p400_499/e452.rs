use crate::utils::{
    FIArray::FIArrayI64, fenwick::FenwickTree,
    multiplicative_function_summation::general_divisor_summatory_alt,
};

const N: i64 = 1e9 as _;
const MOD: i64 = 1_234_567_891;
// bruh
// sum of coefficients up to n in \zeta(s)^n
pub fn main() {
    //dense_pseudo_euler_transform_based();
    let start = std::time::Instant::now();
    let res = general_divisor_summatory_alt::<MOD>(N, N as _)[N];
    println!("res = {}, took {:?}", res, start.elapsed());
}
// TODO: fix
fn dense_pseudo_euler_transform_based() {
    const SQRT_N: i64 = N.isqrt();
    const x: i64 = 1 + SQRT_N.isqrt();
    const INVS: [i64; 4] = [0, 1, modinv(2), modinv(3)];
    const INV_FACTS: [i64; 4] = [0, 1, modinv(2), modinv(6)];
    const { assert!(x.pow(4) > N) };

    let start = std::time::Instant::now();
    let mut log_zeta = log_zeta(N);
    println!("Finished computing log zeta: {:?}", start.elapsed());
    let len = log_zeta.arr.len();

    for i in (1..len).rev() {
        log_zeta.arr[i] += MOD - log_zeta.arr[i - 1];
        if log_zeta.arr[i] >= MOD {
            log_zeta.arr[i] -= MOD;
        }
        log_zeta.arr[i] *= N % MOD;
        log_zeta.arr[i] %= MOD;
    }
    /* for i in (x..=SQRT_N).rev() {
        let v = log_zeta.arr[(i - 1) as usize];
        if v == 0 {
            continue;
        }
        let mut pv = v;
        let mut e = 1;
        let mut pi = i;
        while pi <= N / i {
            e += 1;
            pi *= i;
            pv *= v;
            pv %= MOD;
            log_zeta[pi] += pv * INVS[e];
            log_zeta[pi] %= MOD;
        }
    } */
    println!(
        "Finished computing the contribution of powers: {:?}",
        start.elapsed()
    );

    let mut v = FIArrayI64::new(N);
    for i in x as usize..=len {
        v.arr[i - 1] += v.arr[i - 2] + log_zeta.arr[i - 1];
        v.arr[i - 1] %= MOD;
    }
    let v = v;

    let mut zeta_n = v.clone();
    for e in &mut zeta_n.arr {
        *e += 1;
    }

    let mut r = dirichlet_mul_zero_prefix(&v, &v, N as _, x as _, x as _);
    for i in x as usize..=len {
        zeta_n.arr[i - 1] += r.arr[i - 1] * INV_FACTS[2];
        zeta_n.arr[i - 1] %= MOD;
    }

    r = dirichlet_mul_zero_prefix(&v, &r, N as _, x as _, SQRT_N as _);

    for i in 1..=len {
        zeta_n.arr[i - 1] += r.arr[i - 1] * INV_FACTS[3];
        zeta_n.arr[i - 1] %= MOD;
    }
    println!(
        "Finished computing exp((n log(zeta))_trunc): {:?}",
        start.elapsed()
    );

    for i in (1..len).rev() {
        zeta_n.arr[i] += MOD - zeta_n.arr[i - 1];
        if zeta_n.arr[i] >= MOD {
            zeta_n.arr[i] -= MOD;
        }
    }
    let mut f = FenwickTree::new(0, 0);
    core::mem::swap(&mut zeta_n.arr, &mut f.0);
    f.construct();
    for q in (2..x).rev() {
        let v = log_zeta.arr[(q - 1) as usize];
        if v == 0 {
            continue;
        }
        let lim = N / q;
        // sparse_mul_unlimited(q,v)
        let mut prev = 0;
        for i in FIArrayI64::keys(lim) {
            let cur = f.sum(log_zeta.get_index(i));
            if cur != prev {
                f.add(log_zeta.get_index(i * q), v * (cur - prev));
                prev = cur;
            }
        }
        /* let mut i = 1;
        while i <= lim / i {
            let cur = f.sum(i - 1);
            if cur != prev {
                f.add(log_zeta.get_index(i * q), cur - prev);
                prev = cur;
            }
            i += 1;
        }
        for j in (1..=lim / i).rev() {
            let cur = f.sum(log_zeta.get_index(lim / j));
            if cur != prev {
                f.add(len - j, cur - prev);
                prev = cur;
            }
        } */
    }
    zeta_n.arr = f.flatten();
    for e in &mut zeta_n.arr {
        *e %= MOD;
        if *e < 0 {
            *e += MOD;
        }
    }
    let res = zeta_n[N];
    println!("res = {res}, took {:?}", start.elapsed());
}

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
const fn modinv(x: i64) -> i64 {
    powmod(x, MOD - 2)
}
const fn icbrt(x: i64) -> i64 {
    let mut rt = 1 << (1 + x.ilog2().div_ceil(3));
    let mut x_div_rt2 = (x / rt) / rt;
    while rt > x_div_rt2 {
        rt = ((rt << 1) + x_div_rt2) / 3;
        x_div_rt2 = (x / rt) / rt;
    }
    rt
}
fn log_zeta(n: i64) -> FIArrayI64 {
    const INVS: [i64; 6] = [0, 1, modinv(2), modinv(3), modinv(4), modinv(5)];
    let rt = n.isqrt();
    let mut zeta = FIArrayI64::unit(n as _);
    let len = zeta.arr.len();

    let mut buffer = zeta.clone();

    let mut ret = FIArrayI64::new(n);
    let x = ((icbrt(rt) + 1) * (n as f64).ln() as i64) as usize; // since primes are sparse, can afford to increase x by logarithmic factor without hurting complexity
    // remove contributions of small primes (first ~n^1/6 of them)
    for p in 2..x {
        let val = zeta.arr[p - 1] - 1;
        if val == 0 {
            //not prime
            continue;
        }
        ret.arr[p - 1] = 1;
        for (i, nk) in buffer.arr.iter().enumerate().rev() {
            if i < p {
                break;
            }
            zeta.arr[i] -= zeta[nk / p as i64];
        }
        zeta.arr[p - 1] = 1;
    }
    //let prime_count = zeta.arr[..x - 1].iter().filter(|&&e| e == 1).count();
    zeta.arr[..x - 1].fill(0);

    for i in x..=len {
        zeta.arr[i - 1] -= 1;
    }

    // zeta now equals zeta_t - 1
    // compute log(zeta_t) using log(x + 1) = x^5 / 5 - x^4 / 4 + x^3 / 3 - x^2 / 2 + x
    // x is zeta_t - 1.
    // in order to not have to deal with rational numbers, we compute 60 * log(zeta_t)
    // and adjust later
    // the contributions of x^4 and x^5 are 0 for essentially all reasonable n

    for i in x..=len {
        ret.arr[i - 1] = zeta.arr[i - 1];
    }
    //let start = std::time::Instant::now();
    let pow_zeta = dirichlet_mul_zero_prefix(&zeta, &zeta, n as _, x - 1, x - 1);
    //dbg!(start.elapsed());
    let z2_pref = pow_zeta.arr.iter().take_while(|&&e| e == 0).count();

    for i in z2_pref..=len {
        ret.arr[i - 1] += MOD - (pow_zeta.arr[i - 1] * INVS[2]) % MOD;
        if ret.arr[i - 1] >= MOD {
            ret.arr[i - 1] -= MOD;
        }
    }

    if z2_pref * z2_pref < n as _ {
        dirichlet_mul_zero_prefix_with_buffer(
            &pow_zeta,
            &pow_zeta,
            n as _,
            &mut buffer,
            z2_pref,
            z2_pref,
        );

        for i in z2_pref..=len {
            ret.arr[i - 1] += MOD - (buffer.arr[i - 1] * INVS[4]) % MOD;
            if ret.arr[i - 1] >= MOD {
                ret.arr[i - 1] -= MOD;
            }
        }
    }
    if (x - 1) * z2_pref < n as _ {
        dirichlet_mul_zero_prefix_with_buffer(
            &zeta,
            &pow_zeta,
            n as _,
            &mut buffer,
            x - 1,
            z2_pref,
        );

        for i in z2_pref..=len {
            ret.arr[i - 1] += (buffer.arr[i - 1] * INVS[3]) % MOD;
            if ret.arr[i - 1] >= MOD {
                ret.arr[i - 1] -= MOD;
            }
        }
        let z3_pref = buffer.arr.iter().take_while(|&&e| e == 0).count();

        if z2_pref * z3_pref < n as _ {
            dirichlet_mul_zero_prefix_with_buffer(
                &pow_zeta, &buffer, n as _, &mut zeta, z2_pref, z3_pref,
            );
            for i in z3_pref..=len {
                ret.arr[i - 1] += (zeta.arr[i - 1] * INVS[5]) % MOD;
                if ret.arr[i - 1] >= MOD {
                    ret.arr[i - 1] -= MOD;
                }
            }
        }
    }
    //dbg!(prime_count + ret.arr[len - 1] / 60); // approximate final result
    // correction phase: get rid of contributions of prime powers
    for i in (x + 1..=len).rev() {
        ret.arr[i - 1] -= ret.arr[i - 2];
    }

    for x in 2..=x {
        let v = ret.arr[x - 1];
        if v == 0 {
            continue;
        }
        let mut e = 1;
        //let mut pv = v;
        let mut px = x;
        while px <= n as usize / x {
            e += 1;
            px *= x;
            //pv *= v;
            let ret_px = &mut ret[px as _];
            *ret_px += /* pv * */ modinv(e);
            if *ret_px >= MOD {
                *ret_px -= MOD;
            }
        }
    }
    for i in 1..len {
        ret.arr[i] %= MOD;
        ret.arr[i] += ret.arr[i - 1];
        ret.arr[i] %= MOD;
        if ret.arr[i] < 0 {
            ret.arr[i] += MOD;
        }
    }
    ret
}
#[must_use]
pub fn dirichlet_mul_zero_prefix(
    F: &FIArrayI64,
    G: &FIArrayI64,
    n: usize,
    prefix_f: usize,
    prefix_g: usize,
) -> FIArrayI64 {
    assert!(prefix_f > 0);
    assert!(prefix_g > 0);
    assert!(prefix_f <= prefix_g);
    let mut H = FIArrayI64::new(n as _);
    let len = H.arr.len();
    let rt_n = H.isqrt;

    let real_pref_f = if prefix_f <= rt_n as _ {
        prefix_f
    } else {
        n / (len - prefix_f) - 1
    };
    let real_pref_g = if prefix_g <= rt_n as _ {
        prefix_g
    } else {
        n / (len - prefix_g) - 1
    };

    if real_pref_f * real_pref_g >= n {
        return H;
    }

    let to_ord = |x| {
        if x <= rt_n as _ { x } else { len + 1 - (n / x) }
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

    let prefix_h = to_ord(real_pref_f * real_pref_g) - 1;

    for k in prefix_f..=len {
        let z = len + 1 - k;
        if k >= prefix_h {
            for x in prefix_f.. {
                let y_lo_ord = 1 + to_ord(x).max(to_ord(z));
                let y_hi_ord = to_ord(n / (x * z));
                if y_hi_ord < y_lo_ord {
                    break;
                }
                propogate((x, x), (y_lo_ord, y_hi_ord), (k, k));
                propogate((y_lo_ord, y_hi_ord), (x, x), (k, k));
            }
        }

        let x = k;
        for y in prefix_f..k {
            let z_lo_ord = to_ord(x * y);
            let z_hi_ord = to_ord(n / x);
            if z_hi_ord < z_lo_ord {
                break;
            }
            propogate((x, x), (y, y), (z_lo_ord, z_hi_ord));
            propogate((y, y), (x, x), (z_lo_ord, z_hi_ord));
        }

        if prefix_g <= x && x <= rt_n as _ {
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

pub fn dirichlet_mul_zero_prefix_with_buffer(
    F: &FIArrayI64,
    G: &FIArrayI64,
    n: usize,
    H: &mut FIArrayI64,
    prefix_f: usize,
    prefix_g: usize,
) {
    assert!(prefix_f > 0);
    assert!(prefix_g > 0);
    assert!(prefix_f <= prefix_g);

    H.arr.fill(0);
    let len = H.arr.len();
    if prefix_f == len || prefix_g == len {
        return;
    }
    let rt_n = n.isqrt();

    let real_pref_f = if prefix_f <= rt_n {
        prefix_f
    } else {
        n / (len - prefix_f) - 1
    };
    let real_pref_g = if prefix_g <= rt_n {
        prefix_g
    } else {
        n / (len - prefix_g) - 1
    };

    if real_pref_f * real_pref_g >= n {
        return;
    }

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
    let prefix_h = to_ord(real_pref_f * real_pref_g) - 1;

    for k in prefix_f..=len {
        let z = len + 1 - k;
        if k >= prefix_h {
            for x in prefix_f.. {
                let y_lo_ord = 1 + to_ord(x).max(to_ord(z));
                let y_hi_ord = to_ord(n / (x * z));
                if y_hi_ord < y_lo_ord {
                    break;
                }
                propogate((x, x), (y_lo_ord, y_hi_ord), (k, k));
                propogate((y_lo_ord, y_hi_ord), (x, x), (k, k));
            }
        }
        let x = k;
        for y in prefix_f..k {
            let z_lo_ord = to_ord(x * y);
            let z_hi_ord = to_ord(n / x);
            if z_hi_ord < z_lo_ord {
                break;
            }
            propogate((x, x), (y, y), (z_lo_ord, z_hi_ord));
            propogate((y, y), (x, x), (z_lo_ord, z_hi_ord));
        }

        if prefix_g <= x && x <= rt_n {
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
