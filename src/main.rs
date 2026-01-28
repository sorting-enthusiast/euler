#![warn(clippy::pedantic)]
#![allow(long_running_const_eval)]
#![allow(dead_code)]
#![allow(clippy::too_many_lines)]
#![allow(clippy::cast_possible_truncation)]
#![allow(clippy::inline_always)]
#![allow(clippy::cast_sign_loss)]
#![allow(clippy::cast_precision_loss)]
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]
#![allow(clippy::large_stack_arrays)]
use chrono::Local;

use crate::{
    p300_399::e362::mult,
    utils::{
        FIArray::{DirichletFenwickI64, FIArray, FIArrayI64},
        fast_divisor_sums,
        math::iroot,
        multiplicative_function_summation::{count_squarefree, divisor_summatory, mertens, sqf},
        primes::{
            log_zeta::log_zeta,
            primecount::{lucy_fenwick, mertens_min25},
        },
    },
};
pub mod p0_99;
pub mod p100_199;
pub mod p200_299;
pub mod p300_399;
pub mod p400_499;
pub mod p500_599;
pub mod p600_699;
pub mod p700_799;
pub mod p800_899;
pub mod p900_999;
pub mod test;
pub mod test2;
pub mod utils;

// digital root of n is just n mod 9 if n mod 9 != 0, otherwise 9
const fn is_target_little_endian() -> bool {
    u16::from_ne_bytes([1, 0]) == 1
}
// TODO: understand convex hull based lattice point counting, optimize dirichlet mul
pub fn main() {
    const { assert!(is_target_little_endian()) }; // some code relies on this
    println!("Started running at: {} ", Local::now().time());
    //p500_599::e580::main();
    //p800_899::e890::main();
    //test2::main();
    p400_499::e452::main();
    p300_399::e379::main();
    const N: usize = 1e11 as _;
    let start = std::time::Instant::now();
    let mut pi = inverse_pseudo_euler_transform_i64(FIArrayI64::unit(N));
    let mut primes = vec![];
    for p in 2..=pi.isqrt {
        if pi.arr[p - 1] != pi.arr[p - 2] {
            primes.push(p);
        }
    }
    for i in 0..pi.arr.len() {
        pi.arr[i] *= 3;
    }

    let s2 = mult_correction(&pseudo_euler_transform_i64(pi), &primes, |pp, p, e| {
        2 * e as i64 + 1
    });
    let end = start.elapsed();
    dbg!(end, (s2[N] + N as i64) >> 1); //2278024873139110749

    let start = std::time::Instant::now();
    let mut pi = inverse_pseudo_euler_transform_fraction_i64(FIArrayI64::unit(N));
    let mut primes = vec![];
    for p in 2..=pi.isqrt {
        if pi.arr[p - 1] != pi.arr[p - 2] {
            primes.push(p);
        }
    }
    for i in 0..pi.arr.len() {
        pi.arr[i] *= 3;
    }

    let s2 = mult_correction(
        &pseudo_euler_transform_fraction_i64(pi),
        &primes,
        |pp, p, e| 2 * e as i64 + 1,
    );
    let end = start.elapsed();
    dbg!(end, (s2[N] + N as i64) >> 1);
    assert_eq!(lucy_fenwick(N as _), log_zeta(N as _));
    println!("summing mobius");
    let start = std::time::Instant::now();
    let s1 = mertens(N);
    let end = start.elapsed();
    dbg!(end, s1[N]);
    let start = std::time::Instant::now();
    let s1 = mertens_min25(N);
    let end = start.elapsed();
    dbg!(end, s1[N]);
    let start = std::time::Instant::now();
    let mut pi = inverse_pseudo_euler_transform_fraction_i64(FIArrayI64::unit(N));
    let mut primes = vec![];
    for p in 2..=pi.isqrt {
        if pi.arr[p - 1] != pi.arr[p - 2] {
            primes.push(p);
        }
    }
    for e in &mut pi.arr {
        *e *= -1;
    }

    let s1 = mult_correction(
        &pseudo_euler_transform_fraction_i64(pi),
        &primes,
        |pp, p, e| -i64::from(e < 2),
    );
    let end = start.elapsed();
    dbg!(end, s1[N]);
    let start = std::time::Instant::now();
    let mut pi = inverse_pseudo_euler_transform_i64(FIArrayI64::unit(N));
    let mut primes = vec![];
    for p in 2..=pi.isqrt {
        if pi.arr[p - 1] != pi.arr[p - 2] {
            primes.push(p);
        }
    }
    for e in &mut pi.arr {
        *e *= -1;
    }

    let s1 = mult_correction(
        &pseudo_euler_transform_i64(pi),
        &primes,
        |pp, p, e| match e {
            0 => 1,
            1 => -1,
            _ => 0,
        },
    );
    let end = start.elapsed();
    dbg!(end, s1[N]);
    println!("hello and goodbye");
    println!("summing divisor counts");
    dbg!(fast_divisor_sums::divisor_summatory(N));
    let start = std::time::Instant::now();
    let s1 = divisor_summatory(N);
    let end = start.elapsed();
    dbg!(end, s1[N]);
    let start = std::time::Instant::now();
    let u = FIArray::unit(N);
    let s1 = mult(&u, &u);
    let end = start.elapsed();
    dbg!(end, s1[N]);
    let start = std::time::Instant::now();
    let mut pi = inverse_pseudo_euler_transform_i64(FIArrayI64::unit(N));
    let mut primes = vec![];
    for p in 2..=pi.isqrt {
        if pi.arr[p - 1] != pi.arr[p - 2] {
            primes.push(p);
        }
    }
    for i in 0..pi.arr.len() {
        pi.arr[i] *= 2;
    }

    let s2 = mult_correction(&pseudo_euler_transform_i64(pi), &primes, |pp, p, e| {
        e as i64 + 1
    });
    let end = start.elapsed();
    dbg!(end, s2[N]);
    let start = std::time::Instant::now();
    let mut pi = inverse_pseudo_euler_transform_fraction_i64(FIArrayI64::unit(N));
    let mut primes = vec![];
    for p in 2..=pi.isqrt {
        if pi.arr[p - 1] != pi.arr[p - 2] {
            primes.push(p);
        }
    }
    for i in 0..pi.arr.len() {
        pi.arr[i] *= 2;
    }

    let s2 = mult_correction(
        &pseudo_euler_transform_fraction_i64(pi),
        &primes,
        |pp, p, e| e as i64 + 1,
    );
    let end = start.elapsed();
    dbg!(end, s2[N]);

    //assert_eq!(s1, s2);
    println!("hello and goodbye");
    println!("counting sqf");
    let start = std::time::Instant::now();
    let s1 = sqf(N);
    let end = start.elapsed();
    dbg!(end, s1[N]);
    let start = std::time::Instant::now();
    let s2 = count_squarefree(N);
    let end = start.elapsed();
    dbg!(end, s2[N]);
    for i in 0..s2.arr.len() {
        assert_eq!(s1.arr[i] as usize, s2.arr[i]);
    }
    p300_399::e362::main();
    utils::primes::primecount::main();
    //utils::primes::prime_sieves::main();
    println!("Finished running at: {} ", Local::now().time());
}

#[must_use]
pub fn pseudo_euler_transform_i64(a: FIArrayI64) -> FIArrayI64 {
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
            r[(n / z) as _] += vi * (a[n_i / z] - a.arr[i - 2]);
        }
    }
    let mut r = DirichletFenwickI64::from(r);
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
pub fn inverse_pseudo_euler_transform_i64(a: FIArrayI64) -> FIArrayI64 {
    let n = a.x;
    let rt = a.isqrt;
    let len = a.arr.len();
    let mut r = FIArrayI64::new(n);
    let mut a_bit = DirichletFenwickI64::from(a);
    let x = iroot::<3>(n) + 1;
    for i in 2..x {
        let cur = a_bit.bit.sum(i - 1) - 1;
        if cur == 0 {
            continue;
        }
        r.arr[i - 1] = cur;
        a_bit.sparse_mul_at_most_one(i, cur);
    }
    let a = FIArrayI64::from(a_bit);
    for i in x..=len {
        r.arr[i - 1] = a.arr[i - 1] - a.arr[i - 2];
    }
    for i in x..=rt {
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

fn mult_correction(
    d: &FIArrayI64,
    primes: &[usize],
    f: impl Fn(usize, usize, u8) -> i64,
) -> FIArrayI64 {
    struct Correction(FIArrayI64);
    impl Correction {
        fn fill(
            &mut self,
            primes: &[usize],
            lim: usize,
            x: usize,
            y: i64,
            f: &impl Fn(usize, usize, u8) -> i64,
        ) {
            self.0[x as _] += y;
            for (i, &p) in primes.iter().enumerate() {
                if p > lim / p {
                    break;
                }
                let fp = f(p, p, 1);
                let mut prev = fp;
                let mut pp = p * p;
                let mut new_lim = lim / pp;
                for e in 2.. {
                    let cur = f(pp, p, e);
                    let hp = cur - fp * prev;
                    if hp != 0 {
                        self.fill(&primes[i + 1..], new_lim, x * pp, y * hp, f);
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
    let mut correction = Correction(FIArrayI64::new(d.x));
    correction.fill(primes, d.x, 1, 1, &f);
    for i in 1..correction.0.arr.len() {
        correction.0.arr[i] += correction.0.arr[i - 1];
    }
    mult_sparse_i64(d, &correction.0)
}
pub fn mult_sparse_with_buffer_i64(a: &FIArrayI64, b: &FIArrayI64, res: &mut FIArrayI64) {
    unsafe { core::hint::assert_unchecked(a.x == b.x && a.x == res.x) };
    res.arr.fill(0);
    let R2 = a.isqrt;
    let n = a.x;
    let s1 = |ds: &FIArrayI64| {
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
    for &(x, y) in va {
        res.arr[len - R2 / x] -= y * b.arr[R2 - 1];
    }
    for i in 1..len {
        res.arr[i] += res.arr[i - 1];
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
            res.arr[len - i] += y * a.arr[Nx / i - 1];
            i -= 1;
        }
        while i > 0 {
            res.arr[len - i] += y * a.arr[len - x * i];
            i -= 1;
        }
    }
    for &(x, y) in va {
        for j in 1..=R2 / x {
            res.arr[len - j] += y * b.arr[len - x * j];
        }
    }
}
#[must_use]
pub fn mult_sparse_i64(a: &FIArrayI64, b: &FIArrayI64) -> FIArrayI64 {
    let mut res = a.clone();
    mult_sparse_with_buffer_i64(a, b, &mut res);
    res
}
#[must_use]
pub fn mult_i64(a: &FIArrayI64, b: &FIArrayI64) -> FIArrayI64 {
    unsafe { core::hint::assert_unchecked(a.x == b.x) };
    let R2 = a.isqrt as usize;
    let n = a.x as usize;
    let mut res = FIArrayI64::new(n as _);

    let s1 = |ds: &FIArrayI64| {
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
        for &(x, y) in va {
            res[(x * x) as _] += y * y;

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
                res[(x * pa[r].0) as _] -= y * a.arr[pa[r - 1].0 - 1] * 2;
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
// faster by a log factor, but much more susceptible to overflow - multiplies input by a factor of 1296 before reducing
#[must_use]
pub fn pseudo_euler_transform_fraction_i64(a: FIArrayI64) -> FIArrayI64 {
    const INVS: [i64; 4] = [0, 6, 3, 2];

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

    let mut v = FIArrayI64::new(n);
    for i in x..=len {
        v.arr[i - 1] = v.arr[i - 2] + a.arr[i - 1];
    }

    let mut ret = v.clone();
    for e in &mut ret.arr {
        *e = (*e + INVS[1]) * 6 * INVS[1].pow(2);
    }

    let mut v_2 = mult_i64(&v, &v);
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
    let mut ret = DirichletFenwickI64::from(ret);
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
pub fn inverse_pseudo_euler_transform_fraction_i64(a: FIArrayI64) -> FIArrayI64 {
    const INVS: [i64; 4] = [0, 6, 3, 2];
    let mut a = DirichletFenwickI64::from(a);
    let rt = a.isqrt;
    let n = a.x;
    let len = a.bit.0.len();

    let mut ret = FIArrayI64::new(n);

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
    let mut a = FIArrayI64::from(a);

    // a now equals a_t - 1
    // compute log(a_t) using log(x + 1) = x^3 / 3 - x^2 / 2 + x
    // in order to not have to deal with rational numbers, we compute 6 * log(a_t)
    // and adjust later

    for i in x..=len {
        ret.arr[i - 1] = a.arr[i - 1] * INVS[1];
    }

    let mut a_2 = mult_i64(&a, &a);
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

fn mult_simple_i64(a: &FIArrayI64, b: &FIArrayI64) -> FIArrayI64 {
    unsafe { core::hint::assert_unchecked(a.x == b.x) };
    let R2 = a.isqrt;
    let n = a.x;
    let mut res = FIArrayI64::new(n);

    let s1 = |ds: &FIArrayI64| {
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

        for (a, b) in &pb[l..r] {
            res[x * a] += y * b;
        }
        if r != 0 && x * pb[r].0 <= n {
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

        for (a, b) in &pa[l..r] {
            res[x * a] += y * b;
        }
        if r != 0 && x * pa[r].0 <= n {
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

        for i in (1..=Nx / pb[r].0).rev() {
            res.arr[len - i] += y * b[Nx / i];
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

        for i in (1..=Nx / pa[r].0).rev() {
            res.arr[len - i] += y * a[Nx / i];
        }
    }

    res
}
