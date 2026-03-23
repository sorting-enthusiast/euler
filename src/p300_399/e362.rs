use itertools::Itertools;

use crate::{
    incremental_flattening::{
        DynamicPrefixSumUsize, count_squarefree_partial_flatten as count_squarefree, lucy_pet_slow,
    },
    utils::{
        FIArray::{DirichletFenwick, FIArray},
        math::iroot,
        multiplicative_function_summation::{
            pseudo_euler_transform, pseudo_euler_transform_fraction,
        },
    },
};
// 1e16: 2465164852430507540, 1336.4470363s
// 1e15: 190257704293010022, 286.3650447s
// 1e14: 14574188158034831, 56.9039742s
// 1e13: 1107277852610310, 11.0277353s
// 1e12: 83365737381734, 2.0323739s
// 1e11: 6213486362445, 391.9852ms
// 1e10: 457895958010, 79.9067ms
const N: usize = 1e14 as _;
const SQRT_N: usize = N.isqrt();
// fsf is just the pseudo-euler transform of sqf
// one of my favorite problems
pub fn main() {
    {
        let start = std::time::Instant::now();
        let sqf = count_squarefree(N);
        println!(
            "Finished counting squarefree integers: {:?}",
            start.elapsed()
        );
        let fsf = pseudo_euler_transform_fraction(sqf);
        let res = fsf[N] - 1;
        println!("res = {res}, took {:?}", start.elapsed());
    }
    dense_pseudo_euler_transform_based();
    {
        let start = std::time::Instant::now();
        let sqf = count_squarefree(N);
        println!(
            "Finished counting squarefree integers: {:?}",
            start.elapsed()
        );
        let fsf = pseudo_euler_transform_lucy_dense(sqf);
        let res = fsf[N] - 1;
        println!("lucy: res = {res}, took {:?}", start.elapsed());
    }
    {
        let start = std::time::Instant::now();
        let sqf = count_squarefree(N);
        println!(
            "Finished counting squarefree integers: {:?}",
            start.elapsed()
        );
        let fsf = lucy_pet_slow(sqf);
        let res = fsf[N] - 1;
        println!("lucy: res = {res}, took {:?}", start.elapsed());
    }
    /* {
        let start = std::time::Instant::now();
        let sqf = count_squarefree(N);
        println!(
            "Finished counting squarefree integers: {:?}",
            start.elapsed()
        );
        let fsf = pseudo_euler_transform(sqf);
        let res = fsf[N] - 1;
        println!("res = {res}, took {:?}", start.elapsed());
    } */
    initial_approach_fenwick();
    //initial_approach();
}

fn dense_pseudo_euler_transform_based() {
    const x: usize = 1 + SQRT_N.isqrt();
    const fn inv_odd(mut k: usize) -> usize {
        let mut exp = (1u64 << 63) - 1;

        let mut r: usize = 1;
        while exp > 1 {
            r = r.overflowing_mul(k).0;
            k = k.overflowing_mul(k).0;
            exp >>= 1;
        }
        r.overflowing_mul(k).0
    }

    const INVS: [usize; 4] = [0, 2, 1, inv_odd(3) << 1];
    const { assert!(x.pow(4) > N) };

    let start = std::time::Instant::now();
    let mut sqf = count_squarefree(N);
    println!(
        "Finished counting squarefree integers: {:?}",
        start.elapsed()
    );
    let len = sqf.arr.len();

    for i in (1..len).rev() {
        sqf.arr[i] -= sqf.arr[i - 1];
        sqf.arr[i] *= INVS[1];
    }
    sqf.arr[0] *= INVS[1]; // kinda pointless tbh
    for i in (x..=SQRT_N).rev() {
        let v = sqf.arr[i - 1];
        if v == 0 {
            continue;
        }
        let mut e = 1;
        let mut pi = i;
        while pi <= N / i {
            e += 1;
            pi *= i;
            sqf[pi] += INVS[e];
        }
    }
    println!(
        "Finished computing the contribution of powers: {:?}",
        start.elapsed()
    );

    let mut v = FIArray::new(N);
    for i in x..=len {
        v.arr[i - 1] += v.arr[i - 2] + sqf.arr[i - 1];
    }
    let v = v;

    let mut fsf = v.clone();
    for e in &mut fsf.arr {
        *e = (*e + INVS[1]) * 2 * INVS[1].pow(2);
    }

    let mut r = mult(&v, &v); //dirichlet_mul_zero_prefix(&v, &v, N, x, x);
    for i in x..=len {
        fsf.arr[i - 1] += r.arr[i - 1] * INVS[1];
    }
    r = mult_sparse(&v, &r); //dirichlet_mul_zero_prefix(&v, &r, N, x, SQRT_N);

    for i in 1..=len {
        fsf.arr[i - 1] += r.arr[i - 1] * const { inv_odd(3) };
        fsf.arr[i - 1] /= const { 2 * INVS[1].pow(3) };
    }
    println!("Finished computing exp(v): {:?}", start.elapsed());
    let mut fsf = DirichletFenwick::from(fsf);
    for q in (2..x).rev() {
        if sqf.arr[q - 1] == 0 {
            continue;
        }
        fsf.sparse_mul_unlimited(q, 1);
    }
    let res = fsf.bit.sum(len - 1) - 1;
    println!("res = {res}, took {:?}", start.elapsed());
}

// can optimize using fenwick trees to around O(n^3/4 logn)
fn initial_approach() {
    let start = std::time::Instant::now();
    let mut sqf = count_squarefree(N);
    let mut fsf = FIArray::new(N);
    sqf.adjacent_difference();
    fsf.arr[0] = 1;
    fsf.arr[SQRT_N..].copy_from_slice(&sqf.arr[SQRT_N..]);
    fsf.partial_sum();
    let sqf = sqf;
    let keys = FIArray::keys(N).collect_vec().into_boxed_slice();
    let len = keys.len();

    for q in 2..=SQRT_N {
        if sqf.arr[q - 1] == 0 {
            continue;
        }
        for (i, &v) in keys[q - 1..=len - q].iter().enumerate() {
            fsf.arr[i + q - 1] += fsf[v / q];
        }
        fsf.arr[len - 1] += fsf.arr[len - q];
    }

    let res = fsf.arr[len - 1] - 1;
    println!("res = {res}, took {:?}", start.elapsed());
}
// O(n^3/4 logn) time
fn initial_approach_fenwick() {
    let start = std::time::Instant::now();

    let mut sqf = count_squarefree(N);
    let mut fsf = FIArray::new(N);
    sqf.adjacent_difference();
    fsf.arr[0] = 1;
    fsf.arr[SQRT_N..].copy_from_slice(&sqf.arr[SQRT_N..]);
    fsf.partial_sum();
    let sqf = sqf;
    let mut fsf = DirichletFenwick::from(fsf);

    for q in (2..=SQRT_N).rev() {
        if sqf.arr[q - 1] != 0 {
            fsf.sparse_mul_unlimited(q, 1);
        }
    }
    let res = fsf.get_prefix(N) - 1;
    println!("res = {res}, took {:?}", start.elapsed());
}

// assumes a and b are distinct
pub fn mult_sparse_with_buffer(a: &FIArray, b: &FIArray, res: &mut FIArray) {
    unsafe { core::hint::assert_unchecked(a.x == b.x && a.x == res.x) };
    res.arr.fill(0);
    let R2 = a.isqrt;
    let n = a.x;
    let s1 = |ds: &FIArray| {
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
            res[x * pa[r].0] -= y * a.arr[pa[r - 1].0 - 1];
        }
    }
    for &(x, y) in va {
        res.arr[len - R2 / x] -= y * b.arr[R2 - 1];
    }
    res.partial_sum();
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
pub fn mult_sparse(a: &FIArray, b: &FIArray) -> FIArray {
    let mut res = a.clone();
    mult_sparse_with_buffer(a, b, &mut res);
    res
}
// credit to negiizhao - https://loj.ac/s/1214183
#[must_use]
pub fn mult(a: &FIArray, b: &FIArray) -> FIArray {
    unsafe { core::hint::assert_unchecked(a.x == b.x) };
    let R2 = a.isqrt;
    let n = a.x;
    let mut res = FIArray::new(n);

    let s1 = |ds: &FIArray| {
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
        res.partial_sum();
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
pub fn pseudo_euler_transform_lucy_dense(mut a: FIArray) -> FIArray {
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

    let mut a_fenwick = DynamicPrefixSumUsize(a.arr, len);
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
