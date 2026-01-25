use itertools::Itertools;

use crate::utils::{
    FIArray::{DirichletFenwick, FIArray, FIArrayU128},
    math::iroot,
    multiplicative_function_summation::{
        count_squarefree, dirichlet_mul_u128, dirichlet_mul_with_buffer_u128,
        pseudo_euler_transform,
    },
};
// 1e16: 2393996858318973775, 1424.4494686s
// 1e15: 190257704293010022, 286.3650447s
// 1e14: 14574188158034831, 56.9039742s
// 1e13: 1107277852610310, 11.0277353s
// 1e12: 83365737381734, 2.0323739s
// 1e11: 6213486362445, 391.9852ms
// 1e10: 457895958010, 79.9067ms
const N: usize = 1e10 as _;
const SQRT_N: usize = N.isqrt();
// fsf is just the pseudo-euler transform of sqf
// one of my favorite problems
pub fn main() {
    dense_pseudo_euler_transform_based();
    let start = std::time::Instant::now();
    let sqf = count_squarefree(N);
    println!(
        "Finished counting squarefree integers: {:?}",
        start.elapsed()
    );
    let fsf = pseudo_euler_transform(&sqf);
    let res = fsf[N] - 1;
    println!("res = {res}, took {:?}", start.elapsed());
    dense_pseudo_euler_transform_based_alt();
    initial_approach_fenwick();
    initial_approach();
}

// Also O(n^2/3) time, but makes fewer expensive calls to dirichlet_mul, and has no overflow issues.
// Instead of only excluding values up to n^1/6, and adding their contributions back the naive way in O(n^2/3) time,
// we can exclude all values up to n^1/4, and add their contributions back using a fenwick tree in O(n^5/8 logn) time
// Therefore, we only need to deal with fractions up to 1/3, instead of 1/5, and factorials up to 3!, instead of 5!,
// reducing the constant used for computing exp from 120*60^5 = 933120 to 6*6^3 = 1296, allowing us to compute much larger values
// without encountering 64-bit overflow
fn dense_pseudo_euler_transform_based() {
    const x: usize = 1 + SQRT_N.isqrt();
    const INVS: [usize; 4] = [0, 6, 3, 2];
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
        *e = (*e + INVS[1]) * 6 * INVS[1].pow(2);
    }

    let mut r = mult(&v, &v); //dirichlet_mul_zero_prefix(&v, &v, N, x, x);
    for i in x..=len {
        fsf.arr[i - 1] += r.arr[i - 1] * 3 * INVS[1];
    }
    r = mult_sparse(&r, &v); //dirichlet_mul_zero_prefix(&v, &r, N, x, SQRT_N);

    for i in 1..=len {
        fsf.arr[i - 1] += r.arr[i - 1];
        fsf.arr[i - 1] /= const { 6 * INVS[1].pow(3) };
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

// have to use u128 to stay in integer arithmetic, O(n^2/3) time solution
fn dense_pseudo_euler_transform_based_alt() {
    const x: usize = 1 + iroot::<3>(SQRT_N);
    const INVS: [u128; 6] = [0, 60, 30, 20, 15, 12];

    let start = std::time::Instant::now();
    let sqf_ = count_squarefree(N);
    let mut sqf = FIArrayU128::new(N as _);
    for (e, &q) in sqf.arr.iter_mut().zip(&sqf_.arr) {
        *e = q as u128;
    }

    let len = sqf.arr.len();
    let mut a_vals = sqf.clone();

    for i in (1..len).rev() {
        a_vals.arr[i] -= a_vals.arr[i - 1];
        a_vals.arr[i] *= INVS[1];
    }
    a_vals.arr[0] *= INVS[1]; // kinda pointless tbh
    for i in (x..=SQRT_N).rev() {
        let v = a_vals.arr[i - 1];
        if v == 0 {
            continue;
        }
        assert_eq!(a_vals.arr[i - 1], INVS[1]);
        let mut e = 1;
        let mut pi = i;
        while pi <= N / i {
            e += 1;
            pi *= i;
            a_vals[pi as _] += INVS[e];
        }
    }
    //println!("bello");
    let mut tmp = FIArrayU128::new(N as _);

    let mut v = FIArrayU128::new(N as _);
    for i in x..=len {
        v.arr[i - 1] += v.arr[i - 2] + a_vals.arr[i - 1];
    }
    let v = v;
    let mut fsf = v.clone();
    for e in &mut fsf.arr {
        *e = (*e + INVS[1]) * 120 * INVS[1].pow(4);
    }

    let mut r = dirichlet_mul_u128(&v, &v, N);
    for i in x..=len {
        fsf.arr[i - 1] += r.arr[i - 1] * 60 * INVS[1].pow(3);
    }
    //println!("bello");

    dirichlet_mul_with_buffer_u128(&r, &v, N, &mut tmp);
    core::mem::swap(&mut r.arr, &mut tmp.arr);

    for i in x..=len {
        fsf.arr[i - 1] += r.arr[i - 1] * 20 * INVS[1].pow(2);
    }
    //println!("bello");

    dirichlet_mul_with_buffer_u128(&r, &v, N, &mut tmp);
    core::mem::swap(&mut r.arr, &mut tmp.arr);
    for i in x..=len {
        fsf.arr[i - 1] += r.arr[i - 1] * 5 * INVS[1];
    }
    //println!("bello");

    dirichlet_mul_with_buffer_u128(&r, &v, N, &mut tmp);
    core::mem::swap(&mut r.arr, &mut tmp.arr);
    for i in 1..=len {
        fsf.arr[i - 1] += r.arr[i - 1];
        fsf.arr[i - 1] /= 120 * INVS[1].pow(5);
    }
    //println!("bello");
    let keys = FIArrayU128::keys(N as _).collect_vec().into_boxed_slice();

    for q in 2..x {
        if a_vals.arr[q - 1] == 0 {
            continue;
        }
        //dbg!(q);
        for (i, &v) in keys[q - 1..].iter().enumerate() {
            fsf.arr[i + q - 1] += fsf[v / q as u128];
        }
    }
    let res = fsf[N as _] - 1;
    println!("res = {res}, took {:?}", start.elapsed());
}

// can optimize using fenwick trees to around O(n^3/4 logn)
fn initial_approach() {
    let start = std::time::Instant::now();
    let sqf = count_squarefree(N);
    let mut fsf = FIArray::eps(N);
    let keys = FIArray::keys(N).collect_vec().into_boxed_slice();
    let len = keys.len();

    let mut res = 0;
    for q in 2..=SQRT_N {
        if sqf.arr[q - 1] == sqf.arr[q - 2] {
            continue;
        }
        for (i, &v) in keys[q - 1..=len - q].iter().enumerate() {
            fsf.arr[i + q - 1] += fsf[v / q];
        }
        res += fsf[N / q];
    }

    /*for q in SQRT_N + 1..=N {
        if squarefree[q as usize] == 0 {
            continue;
        }
        fsf[N] += fsf[N / q];
    }*/
    let mut q = SQRT_N + 1;
    while q <= N {
        let k = N / q;
        let q_max = N / k;

        let cnt = sqf[q_max] - sqf[q - 1];

        res += cnt * fsf[k];

        q = q_max + 1;
    }
    println!("res = {res}, took {:?}", start.elapsed());
}
// O(n^3/4 logn) time
fn initial_approach_fenwick() {
    let start = std::time::Instant::now();
    let sqf = count_squarefree(N);
    let mut fsf = DirichletFenwick::eps(N);

    for q in 2..=SQRT_N {
        if sqf.arr[q - 1] != sqf.arr[q - 2] {
            fsf.sparse_mul_unlimited(q, 1);
        }
    }
    let fsf = FIArray::from(fsf);
    let mut res = fsf[N] - 1;

    let mut q = SQRT_N + 1;
    while q <= N {
        let k = N / q;
        let q_max = N / k;

        let cnt = sqf[q_max] - sqf[q - 1];

        res += cnt * fsf[k];

        q = q_max + 1;
    }
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
        res.arr[len - R2 / x] -= y * b.arr[R2 - 1];
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
