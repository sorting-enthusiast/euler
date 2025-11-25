use crate::utils::FIArray::FIArray;
// TODO: optimize using fenwick trees as described in a comment in this post https://codeforces.com/blog/entry/117783

// based on https://codeforces.com/blog/entry/91632?#comment-802482, https://codeforces.com/blog/entry/117783
// O(n^(2/3)) time, O(n^(1/2)) space. Pretty slow, despite noticeably superior time complexity.
// Likely due to repeated calls to dirichlet_mul, which is not particularly fast.
// Moreover, this function needs 2 times more memory than the O(n^0.75/log(n)) lucy_hedgehog based functions.
// O(n^0.75/log(n)) with low constant factors is better than O(n^(2/3)) with medium constant factors out to quite large n
// uses the fact that the coefficients of the logarithm of the DGF of u(n) = 1 (i.e. the zeta function)
// are exactly 1/k for p^k for some prime p, and 0 otherwise.
// Note: similarly to lucy_hedgehog, this code can be adapted to calculate the sum of totally multiplicative functions
// over the primes, though tbh you should probably just use lucy's algorithm for that.
// TODO: try to speed up the convolution steps more somehow, as they are the main bottleneck
const fn icbrt(x: usize) -> usize {
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
#[must_use]
pub fn log_zeta(n: usize) -> FIArray {
    const INVS: [usize; 6] = [0, 60, 30, 20, 15, 12];
    let rt = n.isqrt();
    let mut zeta = FIArray::unit(n as _);
    let len = zeta.arr.len();

    let mut buffer = zeta.clone();

    let mut ret = FIArray::new(n);
    let x = (icbrt(rt) + 1) * (n as f64).ln() as usize;
    // remove contributions of small primes
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
            zeta.arr[i] -= zeta[nk / p];
        }
        zeta.arr[p - 1] = 1;
    }
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
        ret.arr[i - 1] = zeta.arr[i - 1] * INVS[1];
    }

    let mut pow_zeta = dirichlet_mul_zero_prefix(&zeta, &zeta, n, x - 1, x - 1);
    for i in x..=len {
        ret.arr[i - 1] -= pow_zeta.arr[i - 1] * INVS[2];
    }

    dirichlet_mul_zero_prefix_with_buffer(
        &zeta,
        &pow_zeta,
        n,
        &mut buffer,
        x - 1,
        pow_zeta.arr.iter().take_while(|&&e| e == 0).count(),
    );
    core::mem::swap(&mut pow_zeta.arr, &mut buffer.arr);

    for i in x..=len {
        ret.arr[i - 1] += pow_zeta.arr[i - 1] * INVS[3];
    }

    dirichlet_mul_zero_prefix_with_buffer(
        &zeta,
        &pow_zeta,
        n,
        &mut buffer,
        x - 1,
        pow_zeta.arr.iter().take_while(|&&e| e == 0).count(),
    );
    core::mem::swap(&mut pow_zeta.arr, &mut buffer.arr);

    for i in x..=len {
        ret.arr[i - 1] -= pow_zeta.arr[i - 1] * INVS[4];
    }

    dirichlet_mul_zero_prefix_with_buffer(
        &zeta,
        &pow_zeta,
        n,
        &mut buffer,
        x - 1,
        pow_zeta.arr.iter().take_while(|&&e| e == 0).count(),
    );
    core::mem::swap(&mut pow_zeta.arr, &mut buffer.arr);

    for i in x..=len {
        ret.arr[i - 1] += pow_zeta.arr[i - 1] * INVS[5];
    }

    // correction phase: get rid of contributions of prime powers
    for i in (x + 1..=len).rev() {
        ret.arr[i - 1] -= ret.arr[i - 2];
    }

    for x in x..=rt {
        let v = ret.arr[x - 1] / 60;
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
    for i in 1..len {
        if i >= x - 1 {
            ret.arr[i] /= 60;
        }
        ret.arr[i] += ret.arr[i - 1];
    }
    ret
}

// identical to log_zeta, but convolutions are reordered in order to maximise shared 0 prefix, minor speedup for large n
// 1e17: 14411.9862525s
// 1e16: 2035.5664288s
// 1e15: 362.4408137s
// 1e14: 65.0882843s
// 1e13: 12.9792081s
// 1e12: 2.5241903s
// 1e11: 477.0335ms
// 1e10: 90.5864ms
// 1e9: 17.5497ms
// 1e8: 3.5963ms
// can try to write version which only computes final result:
// TODO: better exploit shared zero prefix
#[must_use]
pub fn log_zeta_reordered(n: usize) -> FIArray {
    const INVS: [usize; 6] = [0, 60, 30, 20, 15, 12];
    let rt = n.isqrt();
    let mut zeta = FIArray::unit(n as _);
    let len = zeta.arr.len();

    let mut buffer = zeta.clone();

    let mut ret = FIArray::new(n);
    let x = ((icbrt(rt) + 1) * (n as f64).ln() as usize); // since primes are sparse, can afford to increase x by logarithmic factor without hurting complexity
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
            zeta.arr[i] -= zeta[nk / p];
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
        ret.arr[i - 1] = zeta.arr[i - 1] * INVS[1];
    }
    //let start = std::time::Instant::now();
    let pow_zeta = dirichlet_mul_zero_prefix(&zeta, &zeta, n, x - 1, x - 1);
    //dbg!(start.elapsed());
    let z2_pref = pow_zeta.arr.iter().take_while(|&&e| e == 0).count();

    for i in z2_pref..=len {
        ret.arr[i - 1] -= pow_zeta.arr[i - 1] * INVS[2];
    }

    if z2_pref * z2_pref < n {
        dirichlet_mul_zero_prefix_with_buffer(
            &pow_zeta,
            &pow_zeta,
            n,
            &mut buffer,
            z2_pref,
            z2_pref,
        );

        for i in z2_pref..=len {
            ret.arr[i - 1] -= buffer.arr[i - 1] * INVS[4];
        }
    }
    if (x - 1) * z2_pref < n {
        dirichlet_mul_zero_prefix_with_buffer(&zeta, &pow_zeta, n, &mut buffer, x - 1, z2_pref);

        for i in z2_pref..=len {
            ret.arr[i - 1] += buffer.arr[i - 1] * INVS[3];
        }
        let z3_pref = buffer.arr.iter().take_while(|&&e| e == 0).count();

        if z2_pref * z3_pref < n {
            dirichlet_mul_zero_prefix_with_buffer(
                &pow_zeta, &buffer, n, &mut zeta, z2_pref, z3_pref,
            );
            for i in z3_pref..=len {
                ret.arr[i - 1] += zeta.arr[i - 1] * INVS[5];
            }
        }
    }
    //dbg!(prime_count + ret.arr[len - 1] / 60); // approximate final result
    // correction phase: get rid of contributions of prime powers
    for i in (x + 1..=len).rev() {
        ret.arr[i - 1] -= ret.arr[i - 2];
    }

    for x in x..=rt {
        let v = ret.arr[x - 1];
        if v == 0 {
            continue;
        }
        let mut e = 1;
        //let mut pv = v;
        let mut px = x;
        while px <= n / x {
            e += 1;
            px *= x;
            //pv *= v;

            ret[px] -= /* pv * */ INVS[e];
        }
    }
    for i in 1..len {
        if i >= x - 1 {
            ret.arr[i] /= 60;
        }
        ret.arr[i] += ret.arr[i - 1];
    }
    ret
}

// TODO: fix
/* pub fn log_zeta_reordered_single(n: usize) -> usize {
    const INVS: [usize; 6] = [0, 60, 30, 20, 15, 12];
    let rt = n.isqrt();
    let mut zeta = FIArray::unit(n as _);
    let len = zeta.arr.len();

    let mut buffer = zeta.clone();

    let mut ret = 0;
    let x = {
        let mut x = 2;
        let mut x_cubed = 8;
        while x_cubed <= rt {
            x += 1;
            x_cubed += 3 * x * (x - 1) + 1;
        }
        x
    } * (n as f64).ln() as usize; // since primes are sparse, can afford to increase x by logarithmic factor without hurting complexity
    // remove contributions of small primes (first ~n^1/6 of them)
    for p in 2..x {
        let val = zeta.arr[p - 1] - 1;
        if val == 0 {
            //not prime
            continue;
        }
        ret += 1;
        for (i, nk) in buffer.arr.iter().enumerate().rev() {
            if i < p {
                break;
            }
            zeta.arr[i] -= zeta[nk / p];
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

    ret += zeta.arr[len - 1] * INVS[1];

    //let start = std::time::Instant::now();
    let pow_zeta = dirichlet_mul_zero_prefix(&zeta, &zeta, n, x - 1, x - 1);
    //dbg!(start.elapsed());
    let z2_pref = pow_zeta.arr.iter().take_while(|&&e| e == 0).count();

    ret -= pow_zeta.arr[len - 1] * INVS[2];

    if z2_pref * z2_pref < n {
        ret -= INVS[4]
            * dirichlet_mul_single_zero_prefix_usize(&pow_zeta, &pow_zeta, n, z2_pref, z2_pref);
    }
    if (x - 1) * z2_pref < n {
        //if (x - 1) * z2_pref * z2_pref < n {
        println!("hello 1");
        dirichlet_mul_zero_prefix_with_buffer(&zeta, &pow_zeta, n, &mut buffer, x - 1, z2_pref);

        ret += buffer.arr[len - 1] * INVS[3];

        let z3_pref = buffer.arr.iter().take_while(|&&e| e == 0).count();

        if z2_pref * z3_pref < n {
            ret += dirichlet_mul_single_zero_prefix_usize(&pow_zeta, &buffer, n, z2_pref, z3_pref)
                * INVS[5];
        }
        /* } else {
        println!("hello 2");

        ret +=
            dirichlet_mul_single_zero_prefix_usize(&zeta, &pow_zeta, n, x - 1, z2_pref) * INVS[3];
        } */
    }
    //dbg!(prime_count + ret.arr[len - 1] / 60); // approximate final result
    // correction phase: get rid of contributions of prime powers
    let primes = sift(rt as _);

    for x in primes
        .into_iter()
        .filter_map(|p| (p as usize >= x).then_some(p as usize))
    {
        let mut e = 1;
        let mut px = x;
        while px <= n / x {
            e += 1;
            px *= x;
            ret -= INVS[e];
        }
    }

    ret / 60
}
 */
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
