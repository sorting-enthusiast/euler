use crate::{
    incremental_flattening::DynamicPrefixSumU128,
    utils::{
        FIArray::{DirichletFenwickU64Mod, DirichletFenwickU128, FIArrayU64, FIArrayU128},
        math::iroot,
        multiplicative_function_summation::{dirichlet_mul_single_u64, dirichlet_mul_single_u128},
    },
};

const N: usize = 1e15 as _;
const MOD: u64 = 1e9 as u64 + 7;
const N_MOD: u64 = N as u64 % MOD;

pub fn main() {
    let start = std::time::Instant::now();
    let mut u = FIArrayU64::unit(N);
    for e in &mut u.arr {
        *e %= MOD;
    }
    let f1s = pseudo_euler_transform_mod(u);
    let mut log_deriv = FIArrayU64::unit(N);
    log_deriv.adjacent_difference();
    log_deriv.arr[0] = 0;
    for n in 2..=N.isqrt() {
        let mut nn = n;
        while nn <= N / n {
            nn *= n;
            log_deriv[nn] += 1;
        }
    }
    for i in 1..log_deriv.arr.len() {
        log_deriv.arr[i] %= MOD;
        log_deriv.arr[i] += log_deriv.arr[i - 1];
        if log_deriv.arr[i] >= MOD {
            log_deriv.arr[i] -= MOD;
        }
    }

    let deriv = dirichlet_mul_single_u64(&f1s, &log_deriv) % MOD;

    let res = (((N_MOD + 1) * f1s[N] as u64) % MOD + MOD - (1 + deriv) % MOD) % MOD;
    println!("res = {res}, took {:?}", start.elapsed());

    let start = std::time::Instant::now();
    let mut u = FIArrayU64::unit(N);
    for e in &mut u.arr {
        *e %= MOD;
    }
    let f1s = pseudo_euler_transform_mod_alt(u);
    let mut log_deriv = FIArrayU64::unit(N);
    log_deriv.adjacent_difference();
    log_deriv.arr[0] = 0;
    for n in 2..=N.isqrt() {
        let mut nn = n;
        while nn <= N / n {
            nn *= n;
            log_deriv[nn] += 1;
        }
    }
    for i in 1..log_deriv.arr.len() {
        log_deriv.arr[i] %= MOD;
        log_deriv.arr[i] += log_deriv.arr[i - 1];
        if log_deriv.arr[i] >= MOD {
            log_deriv.arr[i] -= MOD;
        }
    }

    let deriv = dirichlet_mul_single_u64(&f1s, &log_deriv) % MOD;

    let res = (((N_MOD + 1) * f1s[N] as u64) % MOD + MOD - (1 + deriv) % MOD) % MOD;
    println!("res = {res}, took {:?}", start.elapsed());

    let start = std::time::Instant::now();
    let mut u = FIArrayU64::unit(N);
    for e in &mut u.arr {
        *e %= MOD;
    }
    let f1s = pseudo_euler_transform_mod_alt2(u);
    let mut log_deriv = FIArrayU64::unit(N);
    log_deriv.adjacent_difference();
    log_deriv.arr[0] = 0;
    for n in 2..=N.isqrt() {
        let mut nn = n;
        while nn <= N / n {
            nn *= n;
            log_deriv[nn] += 1;
        }
    }
    for i in 1..log_deriv.arr.len() {
        log_deriv.arr[i] %= MOD;
        log_deriv.arr[i] += log_deriv.arr[i - 1];
        if log_deriv.arr[i] >= MOD {
            log_deriv.arr[i] -= MOD;
        }
    }

    let deriv = dirichlet_mul_single_u64(&f1s, &log_deriv) % MOD;

    let res = (((N_MOD + 1) * f1s[N] as u64) % MOD + MOD - (1 + deriv) % MOD) % MOD;
    println!("res = {res}, took {:?}", start.elapsed());
}
fn mult_mod(a: &FIArrayU64, b: &FIArrayU64) -> FIArrayU64 {
    unsafe { core::hint::assert_unchecked(a.x == b.x) };
    let R2 = a.isqrt;
    let n = a.x;
    let mut res = FIArrayU64::new(n);
    let add = |v: &mut u64, a: u64, b: u64| {
        *v += (a * b) % MOD;
        if *v >= MOD {
            *v -= MOD;
        }
    };
    let sub = |v: &mut u64, a: u64, b: u64| {
        *v += MOD - (a * b) % MOD;
        if *v >= MOD {
            *v -= MOD;
        }
    };

    let s1 = |ds: &FIArrayU64| {
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
            if res.arr[i] >= MOD {
                res.arr[i] -= MOD;
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
            if res.arr[i] >= MOD {
                res.arr[i] -= MOD;
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
const Reciprocal: [u64; 128] = const {
    let mut ret = [0; 128];
    ret[1] = 1;
    let mut i = 2;
    while i < ret.len() {
        ret[i] = ((MOD - MOD / i as u64) * ret[MOD as usize % i]) % MOD;
        i += 1;
    }
    ret
};
fn pseudo_euler_transform_mod(a: FIArrayU64) -> FIArrayU64 {
    let mut a = a;
    let len = a.arr.len();
    let n = a.x;
    let rt = a.isqrt;
    let x = iroot::<4>(n) + 1;

    for i in (1..len).rev() {
        a.arr[i] += MOD - a.arr[i - 1];
        if a.arr[i] >= MOD {
            a.arr[i] -= MOD;
        }
    }
    for i in (x..=rt).rev() {
        let v = a.arr[i - 1];
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
            *entry += (pv * Reciprocal[e]) % MOD;
            if *entry >= MOD {
                *entry -= MOD;
            }
        }
    }

    let mut v = FIArrayU64::new(n);
    for i in x..=len {
        v.arr[i - 1] = v.arr[i - 2] + a.arr[i - 1];
        if v.arr[i - 1] >= MOD {
            v.arr[i - 1] -= MOD;
        }
    }

    let mut ret = v.clone();
    for e in &mut ret.arr {
        *e += 1;
        if *e >= MOD {
            *e -= MOD;
        }
    }

    let mut v_2 = mult_mod(&v, &v);
    for i in rt + 1..=len {
        ret.arr[i - 1] += (v_2.arr[i - 1] * Reciprocal[2]) % MOD;
        if ret.arr[i - 1] >= MOD {
            ret.arr[i - 1] -= MOD;
        }
    }
    {
        //v_2 = mult_sparse(&v, &v_2);
        v.arr[rt..].fill(0);
        for i in x..=rt {
            if v.arr[i - 1] != v.arr[i - 2] {
                let mut y = v.arr[i - 1] + MOD - v.arr[i - 2];
                if y >= MOD {
                    y -= MOD;
                }
                for j in 1..=rt / i {
                    v.arr[len - j] += (y * v_2.arr[len - i * j]) % MOD;
                    if v.arr[len - j] >= MOD {
                        v.arr[len - j] -= MOD;
                    }
                }
            }
        }
        v.arr[..rt].fill(0);
        core::mem::swap(&mut v_2, &mut v);
    }

    for i in rt + 1..=len {
        ret.arr[i - 1] += (v_2.arr[i - 1] * Reciprocal[6]) % MOD;
        if ret.arr[i - 1] >= MOD {
            ret.arr[i - 1] -= MOD;
        }
    }
    let mut ret = DirichletFenwickU64Mod::<MOD>::from(ret);
    for i in (2..x).rev() {
        let ai = a.arr[i - 1];
        if ai == 0 {
            continue;
        }
        ret.sparse_mul_unlimited(i, ai);
    }
    ret.into()
}

fn pseudo_euler_transform_mod_alt(a: FIArrayU64) -> FIArrayU64 {
    let mut a = a;
    let len = a.arr.len();
    let n = a.x;
    let rt = a.isqrt;
    let x = iroot::<5>(n) + 1;

    for i in (1..len).rev() {
        a.arr[i] += MOD - a.arr[i - 1];
        if a.arr[i] >= MOD {
            a.arr[i] -= MOD;
        }
    }
    for i in (x..=rt).rev() {
        let v = a.arr[i - 1];
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
            *entry += (pv * Reciprocal[e]) % MOD;
            if *entry >= MOD {
                *entry -= MOD;
            }
        }
    }

    let mut v = FIArrayU64::new(n);
    for i in x..=len {
        v.arr[i - 1] = v.arr[i - 2] + a.arr[i - 1];
        if v.arr[i - 1] >= MOD {
            v.arr[i - 1] -= MOD;
        }
    }

    let mut ret = v.clone();
    for e in &mut ret.arr {
        *e += 1;
        if *e >= MOD {
            *e -= MOD;
        }
    }

    let v_2 = mult_mod(&v, &v);
    for i in x..=len {
        ret.arr[i - 1] += (v_2.arr[i - 1] * Reciprocal[2]) % MOD;
        if ret.arr[i - 1] >= MOD {
            ret.arr[i - 1] -= MOD;
        }
    }
    let v_3 = mult_mod(&v, &v_2);
    for i in x..=len {
        ret.arr[i - 1] += (v_3.arr[i - 1] * Reciprocal[6]) % MOD;
        if ret.arr[i - 1] >= MOD {
            ret.arr[i - 1] -= MOD;
        }
    }
    let v_4 = mult_mod(&v_2, &v_2);
    for i in x..=len {
        ret.arr[i - 1] += (v_4.arr[i - 1] * Reciprocal[24]) % MOD;
        if ret.arr[i - 1] >= MOD {
            ret.arr[i - 1] -= MOD;
        }
    }
    let mut ret = DirichletFenwickU64Mod::<MOD>::from(ret);
    for i in (2..x).rev() {
        let ai = a.arr[i - 1];
        if ai == 0 {
            continue;
        }
        ret.sparse_mul_unlimited(i, ai);
    }
    ret.into()
}

fn pseudo_euler_transform_mod_alt2(a: FIArrayU64) -> FIArrayU64 {
    let mut a = a;
    let len = a.arr.len();
    let n = a.x;
    let rt = a.isqrt;
    let x = iroot::<6>(n) + 1;

    for i in (1..len).rev() {
        a.arr[i] += MOD - a.arr[i - 1];
        if a.arr[i] >= MOD {
            a.arr[i] -= MOD;
        }
    }
    for i in (x..=rt).rev() {
        let v = a.arr[i - 1];
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
            *entry += (pv * Reciprocal[e]) % MOD;
            if *entry >= MOD {
                *entry -= MOD;
            }
        }
    }

    let mut v = FIArrayU64::new(n);
    for i in x..=len {
        v.arr[i - 1] = v.arr[i - 2] + a.arr[i - 1];
        if v.arr[i - 1] >= MOD {
            v.arr[i - 1] -= MOD;
        }
    }

    let mut ret = v.clone();
    for e in &mut ret.arr {
        *e += 1;
        if *e >= MOD {
            *e -= MOD;
        }
    }

    let v_2 = mult_mod(&v, &v);
    for i in x..=len {
        ret.arr[i - 1] += (v_2.arr[i - 1] * Reciprocal[2]) % MOD;
        if ret.arr[i - 1] >= MOD {
            ret.arr[i - 1] -= MOD;
        }
    }
    let v_3 = mult_mod(&v, &v_2);
    for i in x..=len {
        ret.arr[i - 1] += (v_3.arr[i - 1] * Reciprocal[6]) % MOD;
        if ret.arr[i - 1] >= MOD {
            ret.arr[i - 1] -= MOD;
        }
    }
    let v_4 = mult_mod(&v_2, &v_2);
    for i in x..=len {
        ret.arr[i - 1] += (v_4.arr[i - 1] * Reciprocal[24]) % MOD;
        if ret.arr[i - 1] >= MOD {
            ret.arr[i - 1] -= MOD;
        }
    }
    let v_5 = mult_mod(&v_2, &v_3);
    for i in x..=len {
        ret.arr[i - 1] += (v_5.arr[i - 1] * Reciprocal[120]) % MOD;
        if ret.arr[i - 1] >= MOD {
            ret.arr[i - 1] -= MOD;
        }
    }
    let mut ret = DirichletFenwickU64Mod::<MOD>::from(ret);
    for i in (2..x).rev() {
        let ai = a.arr[i - 1];
        if ai == 0 {
            continue;
        }
        ret.sparse_mul_unlimited(i, ai);
    }
    ret.into()
}
pub fn solve_ext() {
    // 10^10: res = 27523066955715610432408, took 183.7176ms
    // 10^11: res = 4146155072313131503875642, took 1.0545004s
    // 10^12: res = 614653534334750666012110236, took 5.5429549s
    // 10^13: res = 89833349988159904198204902249, took 30.8048768s
    // 10^14: res = 12963259897153673037053993566833, took 146.3552248s
    // 10^15: res = 1849270757707706829637849882798748, took 742.6764167s

    // 2^47: res = 27085375706684935393819763517531, took 200.979673s
    // 2^50: res = 2386818655053411238089407138215904, took 791.7089829s
    // 2^53: res = 208557046825781087145648840805722969, took 3347.7084533s
    // 2^54: res = 923832323918027518930451439968932656, took 5598.1159512s
    const N: usize = 1 << 45;
    let start = std::time::Instant::now();
    let f1s = pseudo_euler_transform_fraction_u128(FIArrayU128::unit(N));
    let mut log_deriv = FIArrayU128::unit(N);
    log_deriv.adjacent_difference();
    log_deriv.arr[0] = 0;
    for n in 2..=log_deriv.isqrt {
        let mut nn = n;
        while nn <= N / n {
            nn *= n;
            log_deriv[nn] += 1;
        }
    }
    log_deriv.partial_sum();

    let deriv = dirichlet_mul_single_u128(&f1s, &log_deriv);

    let res = (N + 1) as u128 * f1s[N] - deriv - 1;
    println!("res = {res}, took {:?}", start.elapsed());
}
pub fn solve_ext_lucy() {
    const N: usize = 1 << 45;
    let start = std::time::Instant::now();
    let f1s = pseudo_euler_transform_lucy_dense_u128(FIArrayU128::unit(N));
    let mut log_deriv = FIArrayU128::unit(N);
    log_deriv.adjacent_difference();
    log_deriv.arr[0] = 0;
    for n in 2..=log_deriv.isqrt {
        let mut nn = n;
        while nn <= N / n {
            nn *= n;
            log_deriv[nn] += 1;
        }
    }
    log_deriv.partial_sum();

    let deriv = dirichlet_mul_single_u128(&f1s, &log_deriv);

    let res = (N + 1) as u128 * f1s[N] - deriv - 1;
    println!("res = {res}, took {:?}", start.elapsed());
}

fn pseudo_euler_transform_fraction_u128(a: FIArrayU128) -> FIArrayU128 {
    const fn inv_odd(mut k: u128) -> u128 {
        let mut exp = (1u128 << 127) - 1;

        let mut r = 1u128;
        while exp > 1 {
            r = r.overflowing_mul(k).0;
            k = k.overflowing_mul(k).0;
            exp >>= 1;
        }
        r.overflowing_mul(k).0
    }

    const INVS: [u128; 4] = [0, 2, 1, inv_odd(3) << 1];
    let start = std::time::Instant::now();
    let mut a = a;
    let len = a.arr.len();
    let n = a.x;
    let rt = a.isqrt;
    let x = iroot::<4>(n) + 1;

    for i in (1..len).rev() {
        a.arr[i] -= a.arr[i - 1];
        a.arr[i] *= INVS[1];
    }
    //a.arr[0] *= INVS[1]; // kinda pointless tbh
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
    println!("Added power contributions: {:?}", start.elapsed());
    let mut v = FIArrayU128::new(n);
    for i in x..=len {
        v.arr[i - 1] = v.arr[i - 2] + a.arr[i - 1];
    }

    let mut ret = v.clone();
    for e in &mut ret.arr {
        *e = (*e + INVS[1]) * 2 * INVS[1].pow(2);
    }
    println!("Started convolution: {:?}", start.elapsed());
    let mut v_2 = mult_u128(&v, &v);
    println!("Finished convolution: {:?}", start.elapsed());
    for i in x..=len {
        ret.arr[i - 1] += v_2.arr[i - 1] * INVS[1];
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
        ret.arr[i - 1] += v_2.arr[i - 1] * const { inv_odd(3) };
        ret.arr[i - 1] /= const { 2 * INVS[1].pow(3) };
    }
    println!("Started appending small values: {:?}", start.elapsed());
    let mut ret = DirichletFenwickU128::from(ret);
    for i in (2..x).rev() {
        let ai = a.arr[i - 1] / INVS[1];
        if ai == 0 {
            continue;
        }
        ret.sparse_mul_unlimited(i, ai);
    }
    println!("Finished appending small values: {:?}", start.elapsed());

    ret.into()
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
        for i in 1..ds.isqrt {
            if ds.arr[i] != ds.arr[i - 1] {
                vec.push((i + 1, ds.arr[i] - ds.arr[i - 1]));
            }
        }
        vec.push((ds.isqrt + 1, 0));
        vec
    };
    let len = res.arr.len();
    if a == b {
        let pa = s1(a);
        let va = &pa[..pa.len() - 1];
        let mut r = va.len();
        let mut l = 0;
        for &(x, y) in va {
            res[x * x] += y * y;

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
                res.arr[x * pa[i].0 - 1] += 2 * y * pa[i].1;
                i += 1;
            }
            while i != r {
                res.arr[len - Nx / pa[i].0] += 2 * y * pa[i].1;
                i += 1;
            }
            if r != 0 && pa[r].0 <= Nx {
                res[x * pa[r].0] -= y * a.arr[pa[r - 1].0 - 1] * 2;
            }
        }
        for i in 1..len {
            res.arr[i] += res.arr[i - 1];
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
                res.arr[len - i] += 2 * y * a.arr[Nx / i - 1];
                i -= 1;
            }
            while i > 0 {
                res.arr[len - i] += 2 * y * a.arr[len - x * i];
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
            res[i * i] += (a.arr[i - 1] - a.arr[i - 2]) * (b.arr[i - 1] - b.arr[i - 2]);
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
                res[x * pb[r].0] -= y * b.arr[pb[r - 1].0 - 1];
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
                res[x * pa[r].0] -= y * a.arr[pa[r - 1].0 - 1];
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
#[must_use]
pub fn pseudo_euler_transform_lucy_dense_u128(mut a: FIArrayU128) -> FIArrayU128 {
    let x = a.x;
    let xsqrt = a.isqrt;
    let len = a.arr.len();
    let cutoff = xsqrt
        .isqrt()
        .max(2 * iroot::<3>((xsqrt / x.ilog2() as usize).pow(2)));
    let mut sp = a.arr[xsqrt - 1];

    for p in (cutoff + 1..=xsqrt).rev() {
        let sp1 = a.arr[p - 2];
        if sp1 == sp {
            continue;
        }
        let w = sp - sp1;
        //let mut ip = 0;
        for i in (1..=(x / p) / p).rev() {
            a.arr[len - i] += w
                * (if i * p <= xsqrt {
                    a.arr[len - i * p]
                } else {
                    a.arr[(x / (i * p)) - 1]
                } - sp1);
        }
        sp = sp1;
    }

    let mut a_fenwick = DynamicPrefixSumU128(a.arr, len);
    let get_index = |v| -> usize {
        if v == 0 {
            unsafe { core::hint::unreachable_unchecked() };
        } else if v <= xsqrt {
            v - 1
        } else {
            len - (x / v)
        }
    };
    for p in (2..=cutoff).rev() {
        let sp1 = a_fenwick.sum(p - 2);
        if sp1 == sp {
            continue;
        }

        a_fenwick.shrink_flattened_prefix(1 + get_index(p * p));

        let w = sp - sp1;
        let lim = x / p;
        let mut prev = sp1;
        let mut i = p;
        while i <= lim / i {
            let cur = a_fenwick.sum(i - 1);
            if cur != prev {
                a_fenwick.add(get_index(i * p), w * (cur - prev));
                prev = cur;
            }
            i += 1;
        }
        for j in (1..=lim / i).rev() {
            let cur = a_fenwick.sum(get_index(lim / j));
            if cur != prev {
                a_fenwick.add(len - j, w * (cur - prev));
                prev = cur;
            }
        }
        sp = sp1;
    }
    a_fenwick.shrink_flattened_prefix(1);
    if a_fenwick.sum(0) == 0 {
        a_fenwick.inc(0);
    }
    a.arr = a_fenwick.flatten();

    a
}
