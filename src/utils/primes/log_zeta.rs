#![allow(non_snake_case)]
use crate::{
    p300_399::e362::{mult, mult_sparse},
    utils::{
        FIArray::{DirichletFenwick, FIArray},
        math::iroot,
    },
};
// TODO: try implementing the original variant
// based on https://codeforces.com/blog/entry/91632?#comment-802482, https://codeforces.com/blog/entry/117783
// O(n^\frac23  \log^{-(1+c)} n) time, O(n^(1/2)) space. (not entirely sure what the value of c is)
// equivalent to inverse_pseudo_euler_transform_fraction(FIArray::unit(n))
// 1e17: res = 2623557157654233, took 896.8224632s
// 1e16: res = 279238341033925, took 190.0312267s
// 1e15: res = 29844570422669, took 43.0123352s
// 1e14: res = 3204941750802, took 10.3329261s
// 1e13: res = 346065536839, took 2.3627085s
// 1e12: res = 37607912018, took 493.5248ms
// 1e11: res = 4118054813, took 120.5901ms
// 1e10: res = 455052511, took 23.8138ms
#[must_use]
pub fn log_zeta(n: usize) -> FIArray {
    const INVS: [usize; 4] = [0, 6, 3, 2];
    let mut zeta = DirichletFenwick::zeta(n);
    let rt = zeta.isqrt;
    let len = zeta.bit.0.len();

    let mut ret = FIArray::new(n);

    let x = iroot::<4>(n) + 1;
    // remove contributions of small primes
    for p in 2..x {
        if zeta.bit.sum(p - 1) == 1 {
            //not prime
            continue;
        }
        ret.arr[p - 1] = 1;
        zeta.sparse_mul_at_most_one(p, 1);
    }
    zeta.bit.dec(0);
    let mut zeta = FIArray::from(zeta);

    // zeta now equals zeta_t - 1
    // compute log(zeta_t) using log(x + 1) = x^3 / 3 - x^2 / 2 + x
    // in order to not have to deal with rational numbers, we compute 6 * log(zeta_t)
    // and adjust later

    for i in x..=len {
        ret.arr[i - 1] = zeta.arr[i - 1] * INVS[1];
    }

    let pa = {
        let mut vec = vec![];
        for i in x..=rt {
            if zeta.arr[i - 1] != zeta.arr[i - 2] {
                vec.push(i);
            }
        }
        vec.push(rt + 1);
        vec
    };
    let va = &pa[..pa.len() - 1];
    let mut pow_zeta = //mult(&zeta, &zeta);
    {
        let mut res = FIArray::new(n);
        let mut r = va.len();
        let mut l = 0;
        for &x in va {
            res[x * x] += 1;

            l += 1;
            let Nx = n / x;
            while r > l && Nx / pa[r - 1] < r - 1 {
                r -= 1;
            }
            if r < l {
                r = l;
            }
            let X = rt / x;
            let Nx = n / x;

            let mut i = l;
            /* while i != r {
                res[x * pa[i]] += 2;
                i += 1;
            } */
            while i != r && pa[i] <= X {
                res.arr[x * pa[i] - 1] += 2 ;
                i += 1;
            }
            while i != r {
                res.arr[len - Nx / pa[i]] += 2 ;
                i += 1;
            }
            
            if r != 0 && pa[r] <= Nx {
                res[x * pa[r]] -= zeta.arr[pa[r - 1] - 1] * 2;
            }
        }
        for i in 1..len {
            res.arr[i] += res.arr[i - 1];
        }
        r = va.len();
        l = 0;
        for &x in va {
            l += 1;
            let Nx = n / x;
            while r > l && Nx / pa[r - 1] < r - 1 {
                r -= 1;
            }
            if r < l {
                r = l;
            }
            /* for i in (1..=Nx / pa[r]).rev() {
                res[n / i] += 2 * zeta[Nx / i];
            } */
            let mut i = Nx / pa[r];
            let X = rt / x;
            while i > X {
                res.arr[len - i] += 2 * zeta.arr[Nx / i - 1];
                i -= 1;
            }
            while i > 0 {
                res.arr[len - i] += 2 * zeta.arr[len - x * i];
                i -= 1;
            }
        }
        res
    };
    /* let ind = pow_zeta.get_index(x.pow(2));
    assert!(pow_zeta.arr[..ind].iter().all(|&e| e == 0));
    dbg!(ind, x.pow(2), rt - 1, len, len - ind); */

    for i in rt + 1..=len {
        ret.arr[i - 1] -= pow_zeta.arr[i - 1] * INVS[2];
    }
    {
        //pow_zeta = mult_sparse(&zeta, &pow_zeta);

        zeta.arr[rt..].fill(0);
        for &i in va {
            for j in 1..=rt / i {
                //zeta[n / j] += pow_zeta[n / (i * j)];
                zeta.arr[len - j] += pow_zeta.arr[len - i * j];
            }
        }
        //zeta.arr[..rt].fill(0);
        core::mem::swap(&mut pow_zeta, &mut zeta);
    }
    let ind = pow_zeta.get_index(x.pow(3));
    //dbg!(ind, len, len - ind);
    //assert!(pow_zeta.arr[..ind].iter().all(|&e| e == 0)); */
    for i in ind + 1..=len {
        ret.arr[i - 1] += pow_zeta.arr[i - 1] * INVS[3];
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
        let mut px = x;
        while px <= n / x {
            e += 1;
            px *= x;

            ret[px] -= INVS[e];
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

#[must_use]
pub fn log_zeta_2(n: usize) -> FIArray {
    const INVS: [usize; 5] = [0, 12, 6, 4, 3];
    let mut zeta = DirichletFenwick::zeta(n);
    let rt = zeta.isqrt;
    let len = zeta.bit.0.len();

    let mut ret = FIArray::new(n);

    let x = iroot::<5>(n) + 1;
    // remove contributions of small primes
    for p in 2..x {
        if zeta.bit.sum(p - 1) == 1 {
            //not prime
            continue;
        }
        ret.arr[p - 1] = 1;
        zeta.sparse_mul_at_most_one(p, 1);
    }
    zeta.bit.dec(0);
    let zeta = FIArray::from(zeta);

    // zeta now equals zeta_t - 1
    // compute log(zeta_t) using log(x + 1) = x - x^2 / 2 + x^3 / 3 - x^4 / 4
    // in order to not have to deal with rational numbers, we compute 12 * log(zeta_t)
    // and adjust later

    for i in x..=len {
        ret.arr[i - 1] = zeta.arr[i - 1] * INVS[1];
    }
    let zeta_2 = mult(&zeta, &zeta);
    let mut x_pow = x * x;

    for i in ret.get_index(x_pow)..=len {
        ret.arr[i - 1] -= zeta_2.arr[i - 1] * INVS[2];
    }    
    let zeta_3 = mult(&zeta, &zeta_2);
    x_pow *= x;
    for i in ret.get_index(x_pow)..=len {
        ret.arr[i - 1] += zeta_3.arr[i - 1] * INVS[3];
    }
    let zeta_4 = mult(&zeta_2, &zeta_2);
    x_pow *= x;
    for i in ret.get_index(x_pow)..=len {
        ret.arr[i - 1] -= zeta_4.arr[i - 1] * INVS[4];
    }
    // correction phase: get rid of contributions of prime powers
    for i in (x..len).rev() {
        ret.arr[i] -= ret.arr[i - 1];
    }

    for x in x..=rt {
        let v = ret.arr[x - 1] / INVS[1];
        if v == 0 {
            continue;
        }
        let mut e = 1;
        let mut px = x;
        while px <= n / x {
            e += 1;
            px *= x;

            ret[px] -= INVS[e];
        }
    }
    for i in 1..x - 1 {
        ret.arr[i] += ret.arr[i - 1];
    }
    for i in x - 1..len {
        ret.arr[i] /= INVS[1];
        ret.arr[i] += ret.arr[i - 1];
    }
    ret
}

// TODO: finish implementing ecnerwala's approach: sieve up to n^1/4, flatten, and compute P2 and P3: https://codeforces.com/blog/entry/117783
#[must_use]
pub fn log_zeta_single(n: usize) -> usize {
    let mut f = DirichletFenwick::zeta(n);
    let rt = f.isqrt;
    let len = f.bit.0.len();

    let mut ret = 0;

    let x = iroot::<4>(n) + 1;
    // remove contributions of small primes
    for p in 2..x {
        if f.bit.sum(p - 1) == 1 {
            //not prime
            continue;
        }
        ret += 1;
        f.sparse_mul_at_most_one(p, 1);
    }
    let f = FIArray::from(f);
    ret += f[n];

    ret
}

#[must_use]
pub fn dirichlet_mul_zero_prefix(
    F: &FIArray,
    G: &FIArray,
    n: usize,
    prefix_f: usize,
    prefix_g: usize,
) -> FIArray {
    assert!(prefix_f > 0);
    assert!(prefix_g > 0);
    assert!(prefix_f <= prefix_g);
    let mut H = FIArray::new(n as _);
    let len = H.arr.len();
    let rt_n = H.isqrt;

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
        return H;
    }

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
        H.arr[i] += H.arr[i - 1];
    }
    H
}
/*
pub fn dirichlet_mul_zero_prefix_with_buffer(
    F: &FIArray,
    G: &FIArray,
    n: usize,
    H: &mut FIArray,
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

        let t = (f_x1 - f_x0_1) * (g_y1 - g_y0_1);
        H.arr[z0 - 1] += t;
        if let Some(v) = H.arr.get_mut(z1) {
            *v -= t;
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
        H.arr[i] += H.arr[i - 1];
    }
}
 */
