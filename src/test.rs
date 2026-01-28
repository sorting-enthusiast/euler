//https://loj.ac/s/2344195, https://loj.ac/s/1214183
use itertools::Itertools;

use crate::utils::{
    FIArray::{DirichletFenwickU32Mod, FIArrayU32},
    math::iroot,
    multiplicative_function_summation::sum_n_u64,
};

const N: usize = 9913488072488; //1e12 as _;
const MOD: u64 = 1e9 as u64 + 7;
const N_MOD: u64 = N as u64 % MOD;
const Reciprocal: [u64; N.ilog2() as usize + 1] = const {
    let mut ret = [0; N.ilog2() as usize + 1];
    ret[1] = 1;
    let mut i = 2;
    while i < ret.len() {
        ret[i] = ((MOD - MOD / i as u64) * ret[MOD as usize % i]) % MOD;
        i += 1;
    }
    ret
};
pub fn main() {
    let start = std::time::Instant::now();
    let mut zeta_mod = FIArrayU32::new(N);
    for (i, v) in FIArrayU32::keys(N).enumerate() {
        zeta_mod.arr[i] = (v % MOD as usize) as u32;
    }
    let mut id_mod = FIArrayU32::new(N);
    for (i, v) in FIArrayU32::keys(N).enumerate() {
        id_mod.arr[i] = sum_n_u64::<MOD>(v) as u32;
    }
    dbg!(start.elapsed());
    let pi = inverse_pseudo_euler_transform_mod(zeta_mod);
    dbg!(start.elapsed());

    let mut primes = vec![];
    for p in 2..=pi.isqrt {
        if pi.arr[p - 1] != pi.arr[p - 2] {
            primes.push(p);
        }
    }

    dbg!(start.elapsed());
    let mut pisums = inverse_pseudo_euler_transform_mod(id_mod);
    dbg!(start.elapsed());

    for i in 0..pi.arr.len() {
        pisums.arr[i] += MOD as u32 - pi.arr[i];
        if pisums.arr[i] >= MOD as u32 {
            pisums.arr[i] -= MOD as u32;
        }
    }
    for i in 1..pi.arr.len() {
        pisums.arr[i] += 2;
        if pisums.arr[i] >= MOD as u32 {
            pisums.arr[i] -= MOD as u32;
        }
    }

    dbg!(start.elapsed());
    let approx = pseudo_euler_transform_mod(pisums);
    dbg!(start.elapsed());
    let accurate = mult_correction(&approx, &primes, |_, p, e| p as u32 ^ e as u32);
    dbg!(start.elapsed());
    let res = accurate
        .arr
        .into_iter()
        .sorted_unstable()
        .dedup()
        .fold(0, |acc, v| acc ^ v);
    println!("res = {res}, took {:?}", start.elapsed());
}
fn mult_correction(
    d: &FIArrayU32,
    primes: &[usize],
    f: impl Fn(usize, usize, u8) -> u32,
) -> FIArrayU32 {
    struct Correction(FIArrayU32);
    impl Correction {
        fn fill(
            &mut self,
            primes: &[usize],
            lim: usize,
            x: usize,
            y: u64,
            f: &impl Fn(usize, usize, u8) -> u32,
        ) {
            let entry = &mut self.0[x];
            *entry += y as u32;
            if *entry >= MOD as u32 {
                *entry -= MOD as u32;
            }
            for (i, &p) in primes.iter().enumerate() {
                if p > lim / p {
                    break;
                }
                let fp = u64::from(f(p, p, 1));
                let mut prev = fp;
                let mut pp = p * p;
                let mut new_lim = lim / pp;
                for e in 2.. {
                    let cur = u64::from(f(pp, p, e));
                    let mut hp = cur + MOD - (fp * prev) % MOD;
                    if hp >= MOD {
                        hp -= MOD;
                    }
                    if hp != 0 {
                        self.fill(&primes[i + 1..], new_lim, x * pp, (y * hp) % MOD, f);
                    }
                    prev = cur;
                    if p > new_lim {
                        break;
                    }
                    pp *= p;
                    new_lim /= p;
                }
            }
        }
    }
    let mut correction = Correction(FIArrayU32::new(d.x));
    correction.fill(primes, d.x, 1, 1, &f);
    for i in 1..correction.0.arr.len() {
        correction.0.arr[i] += correction.0.arr[i - 1];
        if correction.0.arr[i] >= MOD as u32 {
            correction.0.arr[i] -= MOD as u32;
        }
    }
    mult_sparse_mod(d, &correction.0)
}

// taken from https://loj.ac/s/1214183
fn mult_sparse_mod(a: &FIArrayU32, b: &FIArrayU32) -> FIArrayU32 {
    unsafe { core::hint::assert_unchecked(a.x == b.x) };
    let mut res = a.clone();

    res.arr.fill(0);
    let R2 = a.isqrt;
    let n = a.x;
    let add = |v: &mut u32, a: u32, b: u32| {
        *v += ((u64::from(a) * u64::from(b)) % MOD) as u32;
        if *v >= MOD as u32 {
            *v -= MOD as u32;
        }
    };
    let sub = |v: &mut u32, a: u32, b: u32| {
        *v += MOD as u32 - ((u64::from(a) * u64::from(b)) % MOD) as u32;
        if *v >= MOD as u32 {
            *v -= MOD as u32;
        }
    };
    let s1 = |ds: &FIArrayU32| {
        let mut vec = vec![];
        if ds.arr[0] != 0 {
            vec.push((1, ds.arr[0]));
        }
        for i in 1..R2 {
            if ds.arr[i] != ds.arr[i - 1] {
                let mut v = ds.arr[i];
                sub(&mut v, 1, ds.arr[i - 1]);
                vec.push((i + 1, v));
            }
        }
        vec.push((R2 + 1, 0));
        vec
    };
    let pa = s1(a);
    let pb = s1(b);
    let va = &pa[..pa.len() - 1];
    let vb = &pb[..pb.len() - 1];
    let len = res.arr.len();

    let mut r = va.len();
    for &(x, y) in vb {
        let Nx = n / x;
        while r != 0 && Nx / pa[r - 1].0 < r - 1 {
            r -= 1;
        }
        let mut i = 0;
        let X = R2 / x;
        let Nx = n / x;
        while i != r && pa[i].0 <= X {
            add(&mut res.arr[x * pa[i].0 - 1], y, pa[i].1);
            i += 1;
        }
        while i != r {
            add(&mut res.arr[len - Nx / pa[i].0], y, pa[i].1);
            i += 1;
        }
        if r != 0 && pa[r].0 <= Nx {
            sub(&mut res[x * pa[r].0], y, a.arr[pa[r - 1].0 - 1]);
        }
    }
    for &(x, y) in va {
        sub(&mut res.arr[len - R2 / x], y, b.arr[R2 - 1]);
    }
    for i in 1..len {
        res.arr[i] += res.arr[i - 1];
        if res.arr[i] >= MOD as u32 {
            res.arr[i] -= MOD as u32;
        }
    }
    let mut r = va.len();
    for &(x, y) in vb {
        let Nx = n / x;
        while r != 0 && Nx / pa[r - 1].0 < r - 1 {
            r -= 1;
        }
        let mut i = Nx / pa[r].0;
        let X = R2 / x;
        while i > X {
            add(&mut res.arr[len - i], y, a.arr[Nx / i - 1]);
            i -= 1;
        }
        while i > 0 {
            add(&mut res.arr[len - i], y, a.arr[len - x * i]);
            i -= 1;
        }
    }
    for &(x, y) in va {
        for j in 1..=R2 / x {
            add(&mut res.arr[len - j], y, b.arr[len - x * j]);
        }
    }
    res
}
fn mult_mod(a: &FIArrayU32, b: &FIArrayU32) -> FIArrayU32 {
    unsafe { core::hint::assert_unchecked(a.x == b.x) };
    let R2 = a.isqrt;
    let n = a.x;
    let mut res = FIArrayU32::new(n);
    let add = |v: &mut u32, a: u32, b: u32| {
        *v += ((u64::from(a) * u64::from(b)) % MOD) as u32;
        if *v >= MOD as u32 {
            *v -= MOD as u32;
        }
    };
    let sub = |v: &mut u32, a: u32, b: u32| {
        *v += MOD as u32 - ((u64::from(a) * u64::from(b)) % MOD) as u32;
        if *v >= MOD as u32 {
            *v -= MOD as u32;
        }
    };

    let s1 = |ds: &FIArrayU32| {
        let mut vec = vec![];
        if ds.arr[0] != 0 {
            vec.push((1, ds.arr[0]));
        }
        for i in 1..R2 {
            if ds.arr[i] != ds.arr[i - 1] {
                let mut v = ds.arr[i];
                sub(&mut v, 1, ds.arr[i - 1]);
                vec.push((i + 1, v));
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
        for &(x, y) in va {
            add(&mut res[x * x], y, y);
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
                add(&mut res.arr[x * pa[i].0 - 1], y << 1, pa[i].1);
                i += 1;
            }
            while i != r {
                add(&mut res.arr[len - Nx / pa[i].0], y << 1, pa[i].1);
                i += 1;
            }
            if r != 0 && pa[r].0 <= Nx {
                sub(&mut res[x * pa[r].0], y << 1, a.arr[pa[r - 1].0 - 1]);
            }
        }
        for i in 1..len {
            res.arr[i] += res.arr[i - 1];
            if res.arr[i] >= MOD as u32 {
                res.arr[i] -= MOD as u32;
            }
        }
        r = va.len();
        l = 0;
        for &(x, y) in va {
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
                add(&mut res.arr[len - i], y << 1, a.arr[Nx / i - 1]);
                i -= 1;
            }
            while i > 0 {
                add(&mut res.arr[len - i], y << 1, a.arr[len - x * i]);
                i -= 1;
            }
        }
    } else {
        let pa = s1(a);
        let va = &pa[..pa.len() - 1];
        let pb = s1(b);
        let vb = &pb[..pb.len() - 1];
        add(&mut res.arr[0], a.arr[0], b.arr[0]);
        for i in 2..=R2 {
            let mut a_ = a.arr[i - 1];
            sub(&mut a_, 1, a.arr[i - 2]);
            let mut b_ = b.arr[i - 1];
            sub(&mut b_, 1, b.arr[i - 2]);
            add(&mut res[i * i], a_, b_);
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
                add(&mut res.arr[x * pb[i].0 - 1], y, pb[i].1);
                i += 1;
            }
            while i != r {
                add(&mut res.arr[len - Nx / pb[i].0], y, pb[i].1);
                i += 1;
            }
            if r != 0 && pb[r].0 <= Nx {
                sub(&mut res[x * pb[r].0], y, b.arr[pb[r - 1].0 - 1]);
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
                add(&mut res.arr[x * pa[i].0 - 1], y, pa[i].1);
                i += 1;
            }
            while i != r {
                add(&mut res.arr[len - Nx / pa[i].0], y, pa[i].1);
                i += 1;
            }
            if r != 0 && pa[r].0 <= Nx {
                sub(&mut res[x * pa[r].0], y, a.arr[pa[r - 1].0 - 1]);
            }
        }
        for i in 1..len {
            res.arr[i] += res.arr[i - 1];
            if res.arr[i] >= MOD as u32 {
                res.arr[i] -= MOD as u32;
            }
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
                add(&mut res.arr[len - i], y, b.arr[Nx / i - 1]);
                i -= 1;
            }
            while i > 0 {
                add(&mut res.arr[len - i], y, b.arr[len - x * i]);
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
                add(&mut res.arr[len - i], y, a.arr[Nx / i - 1]);
                i -= 1;
            }
            while i > 0 {
                add(&mut res.arr[len - i], y, a.arr[len - x * i]);
                i -= 1;
            }
        }
    }
    res
}

fn pseudo_euler_transform_mod(a: FIArrayU32) -> FIArrayU32 {
    let mut a = a;
    let len = a.arr.len();
    let n = a.x;
    let rt = a.isqrt;
    let x = iroot::<4>(n) + 1;

    for i in (1..len).rev() {
        a.arr[i] += MOD as u32 - a.arr[i - 1];
        if a.arr[i] >= MOD as u32 {
            a.arr[i] -= MOD as u32;
        }
    }
    for i in (x..=rt).rev() {
        let v = u64::from(a.arr[i - 1]);
        if v == 0 {
            continue;
        }
        let mut e = 1;
        let mut pi = i;
        let mut pv = v;
        while pi <= n / i {
            e += 1;
            pi *= i;
            pv = (pv * v) % MOD;
            let entry = &mut a[pi];
            *entry += ((pv * Reciprocal[e]) % MOD) as u32;
            if *entry >= MOD as u32 {
                *entry -= MOD as u32;
            }
        }
    }

    let mut v = FIArrayU32::new(n);
    for i in x..=len {
        v.arr[i - 1] = v.arr[i - 2] + a.arr[i - 1];
        if v.arr[i - 1] >= MOD as u32 {
            v.arr[i - 1] -= MOD as u32;
        }
    }

    let mut ret = v.clone();
    for e in &mut ret.arr {
        *e += 1;
        if *e >= MOD as u32 {
            *e -= MOD as u32;
        }
    }

    let mut v_2 = mult_mod(&v, &v);
    for i in rt + 1..=len {
        ret.arr[i - 1] += ((u64::from(v_2.arr[i - 1]) * Reciprocal[2]) % MOD) as u32;
        if ret.arr[i - 1] >= MOD as u32 {
            ret.arr[i - 1] -= MOD as u32;
        }
    }
    {
        //v_2 = mult_sparse(&v, &v_2);
        v.arr[rt..].fill(0);
        for i in x..=rt {
            if v.arr[i - 1] != v.arr[i - 2] {
                let mut y = u64::from(v.arr[i - 1]) + MOD - u64::from(v.arr[i - 2]);
                if y >= MOD {
                    y -= MOD;
                }
                for j in 1..=rt / i {
                    v.arr[len - j] += ((y * u64::from(v_2.arr[len - i * j])) % MOD) as u32;
                    if v.arr[len - j] >= MOD as u32 {
                        v.arr[len - j] -= MOD as u32;
                    }
                }
            }
        }
        v.arr[..rt].fill(0);
        core::mem::swap(&mut v_2, &mut v);
    }

    for i in rt + 1..=len {
        ret.arr[i - 1] += ((u64::from(v_2.arr[i - 1]) * Reciprocal[6]) % MOD) as u32;
        if ret.arr[i - 1] >= MOD as u32 {
            ret.arr[i - 1] -= MOD as u32;
        }
    }
    let mut ret = DirichletFenwickU32Mod::<{ MOD as u32 }>::from(ret);
    for i in (2..x).rev() {
        let ai = a.arr[i - 1];
        if ai == 0 {
            continue;
        }
        ret.sparse_mul_unlimited(i, ai);
    }
    ret.into()
}
fn inverse_pseudo_euler_transform_mod(a: FIArrayU32) -> FIArrayU32 {
    let mut a = DirichletFenwickU32Mod::<{ MOD as u32 }>::from(a);
    let rt = a.isqrt;
    let n = a.x;
    let len = a.bit.0.len();

    let mut ret = FIArrayU32::new(n);

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
    let mut a = FIArrayU32::from(a);

    // a now equals a_t - 1
    // compute log(a_t) using log(x + 1) = x^3 / 3 - x^2 / 2 + x

    ret.arr[x - 1..].copy_from_slice(&a.arr[x - 1..]);

    let mut a_2 = mult_mod(&a, &a);
    /* let ind = pow_zeta.get_index(x.pow(2));
    assert!(pow_zeta.arr[..ind].iter().all(|&e| e == 0));
    dbg!(ind, x.pow(2), rt - 1, len, len - ind); */

    for i in rt + 1..=len {
        ret.arr[i - 1] += MOD as u32 - ((u64::from(a_2.arr[i - 1]) * Reciprocal[2]) % MOD) as u32;
        if ret.arr[i - 1] >= MOD as u32 {
            ret.arr[i - 1] -= MOD as u32;
        }
    }
    {
        //pow_zeta = mult_sparse(&zeta, &pow_zeta);

        a.arr[rt..].fill(0);
        for i in x..=rt {
            if a.arr[i - 1] != a.arr[i - 2] {
                let mut y = u64::from(a.arr[i - 1]) + MOD - u64::from(a.arr[i - 2]);
                if y >= MOD {
                    y -= MOD;
                }
                for j in 1..=rt / i {
                    a.arr[len - j] += ((y * u64::from(a_2.arr[len - i * j])) % MOD) as u32;
                    if a.arr[len - j] >= MOD as u32 {
                        a.arr[len - j] -= MOD as u32;
                    }
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
        ret.arr[i - 1] += ((u64::from(a_2.arr[i - 1]) * Reciprocal[3]) % MOD) as u32;
        if ret.arr[i - 1] >= MOD as u32 {
            ret.arr[i - 1] -= MOD as u32;
        }
    }

    // correction phase: get rid of contributions of prime powers
    for i in (x..len).rev() {
        ret.arr[i] += MOD as u32 - ret.arr[i - 1];
        if ret.arr[i] >= MOD as u32 {
            ret.arr[i] -= MOD as u32;
        }
    }

    for x in x..=rt {
        let v = u64::from(ret.arr[x - 1]);
        if v == 0 {
            continue;
        }
        let mut e = 1;
        let mut pv = v;
        let mut px = x;
        while px <= n / x {
            e += 1;
            px *= x;
            pv = (pv * v) % MOD;
            let entry = &mut ret[px];
            *entry += MOD as u32 - ((pv * Reciprocal[e]) % MOD) as u32;
            if *entry >= MOD as u32 {
                *entry -= MOD as u32;
            }
        }
    }
    for i in 1..len {
        ret.arr[i] += ret.arr[i - 1];
        if ret.arr[i] >= MOD as u32 {
            ret.arr[i] -= MOD as u32;
        }
    }
    ret
}
