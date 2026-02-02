use itertools::Itertools;

use crate::utils::{
    FIArray::{DirichletFenwick, DirichletFenwickU128, FIArray, FIArrayI64, FIArrayU128},
    math::iroot,
    multiplicative_function_summation::sum_n_i64,
    primes::{log_zeta::log_zeta, prime_sieves::sift},
};
// TODO: optimize to O(n^2/3)
const N: usize = 2018_2018_2018;
const SQRT_N: usize = N.isqrt();
const MOD: usize = 1e9 as _;
pub fn main() {
    let start = std::time::Instant::now();
    let primes = sift(SQRT_N as _);
    println!("Sieved primes: {:?}", start.elapsed());
    let mut pi = FIArray::new(N);
    let mut pi_sums = FIArrayI64::new(N);
    let keys = FIArray::keys(N).collect_vec();

    for (i, &v) in keys.iter().enumerate() {
        pi.arr[i] = v - 1;
        pi_sums.arr[i] = (sum_n_i64::<{ MOD as _ }>(v as _) + MOD as i64 - 1) % MOD as i64;
    }

    for &p in &primes {
        let p = p as usize;
        let sp = pi.arr[p - 2];
        let sp2 = pi_sums.arr[p - 2];

        for (i, &v) in keys.iter().enumerate().rev() {
            if v < p * p {
                break;
            }
            pi.arr[i] -= pi[v / p] - sp;
            pi_sums.arr[i] -= p as i64 * (pi_sums[v / p] - sp2);
            pi_sums.arr[i] %= MOD as i64;
            if pi_sums.arr[i] < 0 {
                pi_sums.arr[i] += MOD as i64;
            }
        }
    }
    let mut sum = ((1..SQRT_N).fold(0, |acc, k| acc + pi_sums[N / k] as usize)
        - (SQRT_N - 1) * pi_sums[const { N / SQRT_N }] as usize)
        % MOD;
    if sum >= MOD {
        sum -= MOD;
    }
    println!("Initialized sum: {:?}", start.elapsed());

    let mut smooth = DirichletFenwick::eps(N);
    let split = primes.partition_point(|p| p * p * p <= N as u64);
    for &p in &primes[..split] {
        let p = p as usize;
        //smooth.sparse_mul_unlimited(p, 1);
        {
            let x = p;
            let w = 1;
            let lim = smooth.x / x;
            let mut prev = 0;
            let mut i = 1;
            while i <= lim / i {
                let cur = smooth.bit.sum(i - 1);
                if cur != prev {
                    smooth.bit.add(smooth.get_index(i * x), w * (cur - prev));
                    prev = cur;
                }
                i += 1;
            }
            for j in (p..=lim / i).rev() {
                let cur = smooth.get_prefix(lim / j);
                if cur != prev {
                    smooth.bit.add(smooth.bit.0.len() - j, w * (cur - prev));
                    prev = cur;
                }
            }
        }
        sum += p * (smooth.get_prefix(N / p) % MOD);
        sum %= MOD;
    }
    for &p in &primes[split..] {
        // optimization, not required for 1-minute rule
        let p = p as usize;
        sum +=
            p * (N / p - (1..=(N / p) / p).fold(0, |acc, k| acc + pi[(N / p) / k] - pi.arr[p - 1]));
        sum %= MOD;
    }
    println!("res = {sum}, took {:?}", start.elapsed());
    solve_ext();
}

fn solve_ext() {
    // 2^52: ~613s
    // 2^53: res = 1872598915682836389066755424708, took 1029.6792178s
    // 2^54: res = 7347331436455433609288167879201, took 1548.9446494s
    const N: usize = 2018_2018_2018;
    const SQRT_N: usize = N.isqrt();
    const CBRT_N: usize = iroot::<3>(N);
    let start = std::time::Instant::now();
    let pi_sums = inverse_pseudo_euler_transform_fraction(FIArrayU128::id::<0>(N));
    println!("Summed primes: {:?}", start.elapsed());

    let mut sum = (1..=SQRT_N).map(|k| pi_sums[N / k] as u128).sum::<u128>()
        - SQRT_N as u128 * pi_sums.arr[SQRT_N - 1];
    println!("Initialized sum: {:?}", start.elapsed());

    let pi_counts = log_zeta(N);
    println!("Counted primes: {:?}", start.elapsed());

    let mut smooth = DirichletFenwick::eps(N);
    for p in 2..=CBRT_N {
        if pi_counts.arr[p - 1] == pi_counts.arr[p - 2] {
            continue;
        }
        //smooth.sparse_mul_unlimited(p, 1);
        {
            let x = p;
            let w = 1;
            let lim = smooth.x / x;
            let mut prev = 0;
            let mut i = 1;
            while i <= lim / i {
                let cur = smooth.bit.sum(i - 1);
                if cur != prev {
                    smooth.bit.add(smooth.get_index(i * x), w * (cur - prev));
                    prev = cur;
                }
                i += 1;
            }
            // j ends at p instead of 1 since we never use sums up to v >= N / p
            for j in (p..=lim / i).rev() {
                let cur = smooth.get_prefix(lim / j);
                if cur != prev {
                    smooth.bit.add(smooth.bit.0.len() - j, w * (cur - prev));
                    prev = cur;
                }
            }
        }
        sum += p as u128 * smooth.get_prefix(N / p) as u128;
    }
    println!("Added small prime contributions: {:?}", start.elapsed());

    for p in CBRT_N + 1..=SQRT_N {
        if pi_counts.arr[p - 1] == pi_counts.arr[p - 2] {
            continue;
        }
        let Np = N / p;
        let term = (Np + (Np / p) * pi_counts.arr[p - 1]) as u128
            - (1..=Np / p).map(|k| pi_counts[Np / k]).sum::<usize>() as u128;
        sum += p as u128 * term;
    }

    println!("res = {sum}, took {:?}", start.elapsed());
}
fn inverse_pseudo_euler_transform_fraction(a: FIArrayU128) -> FIArrayU128 {
    const INVS: [u128; 4] = [0, 6, 3, 2];
    let mut a = DirichletFenwickU128::from(a);
    let rt = a.isqrt;
    let n = a.x;
    let len = a.bit.0.len();

    let mut ret = FIArrayU128::new(n);

    let x = crate::utils::math::iroot::<4>(n) + 1;
    for i in 2..x {
        let vi = a.bit.sum(i - 1) - 1;
        if vi == 0 {
            continue;
        }
        ret.arr[i - 1] = vi;
        a.sparse_mul_at_most_one(i, vi);
    }
    a.bit.dec(0);
    let mut a = FIArrayU128::from(a);

    // a now equals a_t - 1
    // compute log(a_t) using log(x + 1) = x^3 / 3 - x^2 / 2 + x
    // in order to not have to deal with rational numbers, we compute 6 * log(a_t)
    // and adjust later

    for i in x..=len {
        ret.arr[i - 1] = a.arr[i - 1] * INVS[1];
    }

    let mut a_2 = mult_u128(&a, &a);
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
fn mult_u128(a: &FIArrayU128, b: &FIArrayU128) -> FIArrayU128 {
    unsafe { core::hint::assert_unchecked(a.x == b.x) };
    let R2 = a.isqrt;
    let n = a.x;
    let mut res = FIArrayU128::new(n);

    let s1 = |ds: &FIArrayU128| {
        let mut vec = vec![];
        if ds.arr[0] != 0 {
            vec.push((1, ds.arr[0]));
        }
        for i in 1..R2 {
            if ds.arr[i] != ds.arr[i - 1] {
                vec.push((i + 1, ds.arr[i] - ds.arr[i - 1]));
            }
        }
        vec.push((R2 + 1, 0));
        vec
    };
    let len = res.arr.len();
    if a == b {
        let pa = s1(a);
        let va = &pa[..pa.len() - 1];
        let mut r = va.len();
        let mut l = 0;
        for &(x, fx) in va {
            res[x * x] += fx * fx;

            l += 1;
            let Nx = n / x;
            while r > l && Nx / pa[r - 1].0 < r - 1 {
                r -= 1;
            }
            if r < l {
                r = l;
            }
            let X = R2 / x;
            let Nx = n / x;

            let mut i = l;

            while i != r && pa[i].0 <= X {
                res.arr[x * pa[i].0 - 1] += 2 * fx * pa[i].1;
                i += 1;
            }
            while i != r {
                res.arr[len - Nx / pa[i].0] += 2 * fx * pa[i].1;
                i += 1;
            }
            if r != 0 && pa[r].0 <= Nx {
                res[x * pa[r].0] -= fx * a.arr[pa[r - 1].0 - 1] * 2;
            }
        }
        for i in 1..len {
            res.arr[i] += res.arr[i - 1];
        }
        r = va.len();
        l = 0;
        for &(x, fx) in va {
            l += 1;
            let Nx = n / x;
            while r > l && Nx / pa[r - 1].0 < r - 1 {
                r -= 1;
            }
            if r < l {
                r = l;
            }

            let mut i = Nx / pa[r].0;
            let X = R2 / x;
            while i > X {
                res.arr[len - i] += 2 * fx * a.arr[Nx / i - 1];
                i -= 1;
            }
            while i > 0 {
                res.arr[len - i] += 2 * fx * a.arr[len - x * i];
                i -= 1;
            }
        }
    } else {
        let pa = s1(a);
        let va = &pa[..pa.len() - 1];
        let pb = s1(b);
        let vb = &pb[..pb.len() - 1];
        res.arr[0] += a.arr[0] * b.arr[0];
        for i in 2..=R2 {
            res[(i * i) as _] += (a.arr[i - 1] - a.arr[i - 2]) * (b.arr[i - 1] - b.arr[i - 2]);
        }
        let mut r = vb.len();
        let mut l = 0;
        for &(x, y) in va {
            while pb[l].0 <= x {
                l += 1;
            }
            let Nx = n / x;
            while r > l && Nx / pb[r - 1].0 < r - 1 {
                r -= 1;
            }
            if r < l {
                r = l;
            }
            let X = R2 / x;
            let Nx = n / x;

            let mut i = l;

            while i != r && pb[i].0 <= X {
                res.arr[x * pb[i].0 - 1] += y * pb[i].1;
                i += 1;
            }
            while i != r {
                res.arr[len - Nx / pb[i].0] += y * pb[i].1;
                i += 1;
            }
            if r != 0 && pb[r].0 <= Nx {
                res[(x * pb[r].0) as _] -= y * b.arr[pb[r - 1].0 - 1];
            }
        }
        r = va.len();
        l = 0;
        for &(x, y) in vb {
            while pa[l].0 <= x {
                l += 1;
            }
            let Nx = n / x;
            while r > l && Nx / pa[r - 1].0 < r - 1 {
                r -= 1;
            }
            if r < l {
                r = l;
            }
            let X = R2 / x;
            let Nx = n / x;

            let mut i = l;

            while i != r && pa[i].0 <= X {
                res.arr[x * pa[i].0 - 1] += y * pa[i].1;
                i += 1;
            }
            while i != r {
                res.arr[len - Nx / pa[i].0] += y * pa[i].1;
                i += 1;
            }
            if r != 0 && pa[r].0 <= Nx {
                res[(x * pa[r].0) as _] -= y * a.arr[pa[r - 1].0 - 1];
            }
        }
        for i in 1..len {
            res.arr[i] += res.arr[i - 1];
        }
        r = vb.len();
        l = 0;
        for &(x, y) in va {
            while pb[l].0 <= x {
                l += 1;
            }
            let Nx = n / x;
            while r > l && Nx / pb[r - 1].0 < r - 1 {
                r -= 1;
            }
            if r < l {
                r = l;
            }

            let mut i = Nx / pb[r].0;
            let X = R2 / x;
            while i > X {
                res.arr[len - i] += y * b.arr[Nx / i - 1];
                i -= 1;
            }
            while i > 0 {
                res.arr[len - i] += y * b.arr[len - x * i];
                i -= 1;
            }
        }
        r = va.len();
        l = 0;
        for &(x, y) in vb {
            while pa[l].0 <= x {
                l += 1;
            }
            let Nx = n / x;
            while r > l && Nx / pa[r - 1].0 < r - 1 {
                r -= 1;
            }
            if r < l {
                r = l;
            }

            let mut i = Nx / pa[r].0;
            let X = R2 / x;
            while i > X {
                res.arr[len - i] += y * a.arr[Nx / i - 1];
                i -= 1;
            }
            while i > 0 {
                res.arr[len - i] += y * a.arr[len - x * i];
                i -= 1;
            }
        }
    }
    res
}
