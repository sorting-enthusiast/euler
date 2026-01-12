use itertools::Itertools;

use crate::utils::{
    FIArray::{FIArray, FIArrayU128},
    fenwick::FenwickTreeUsize,
    math::iroot,
    multiplicative_function_summation::{
        count_squarefree, dirichlet_mul_u128, dirichlet_mul_with_buffer_u128,
    },
    primes::log_zeta::dirichlet_mul_zero_prefix,
};
// 1e15: 190257704293010022, 572.8073284s
// 1e14: 14574188158034831, 115.9572641s
// 1e13: 1107277852610310, 21.7458745s
// 1e12: 83365737381734, 3.8677086s
// 1e11: 6213486362445, 735.9992ms
// 1e10: 457895958010, 162.7931ms
const N: usize = 1e11 as _;
const SQRT_N: usize = N.isqrt();
// fsf is just the pseudo-euler transform of sqf
// one of my favorite problems
pub fn main() {
    dense_pseudo_euler_transform_based();
    //dense_pseudo_euler_transform_based_alt();
    //initial_approach_fenwick();
    //initial_approach();
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

    let mut r = dirichlet_mul_zero_prefix(&v, &v, N, x, x);
    for i in x..=len {
        fsf.arr[i - 1] += r.arr[i - 1] * 3 * INVS[1];
    }

    r = dirichlet_mul_zero_prefix(&v, &r, N, x, SQRT_N);

    for i in 1..=len {
        fsf.arr[i - 1] += r.arr[i - 1];
        fsf.arr[i - 1] /= const { 6 * INVS[1].pow(3) };
    }
    println!("Finished computing exp(v): {:?}", start.elapsed());

    for i in (2..len).rev() {
        let j = i & (i + 1);
        if j != 0 {
            fsf.arr[i] -= fsf.arr[j - 1];
        }
    }
    let mut f = FenwickTreeUsize::new(0, 0);
    core::mem::swap(&mut fsf.arr, &mut f.0);

    for q in (2..x).rev() {
        if sqf.arr[q - 1] == 0 {
            continue;
        }
        let lim = N / q;
        // sparse_mul_unlimited(q,1)
        let mut prev = 0;
        /*for i in FIArray::keys(lim) {
            let cur = f.sum(sqf.get_index(i));
            if cur != prev {
                f.add(sqf.get_index(i * q), cur - prev);
                prev = cur;
            }
        }*/
        let mut i = 1;
        while i <= lim / i {
            let cur = f.sum(i - 1);
            if cur != prev {
                f.add(sqf.get_index(i * q), cur - prev);
                prev = cur;
            }
            i += 1;
        }
        for j in (1..=lim / i).rev() {
            let cur = f.sum(sqf.get_index(lim / j));
            if cur != prev {
                f.add(len - j, cur - prev);
                prev = cur;
            }
        }
    }
    fsf.arr = f.flatten();

    let res = fsf[N] - 1;
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

fn test() {
    const N: usize = 1e2 as _;
    const SQRT_N: usize = N.isqrt();

    let mut squarefree = vec![1u8; N + 1];
    squarefree[0] = 0;
    for q in 2..=SQRT_N {
        if squarefree[q * q] == 0 {
            continue;
        }
        for m in (q * q..=N as usize).step_by(q * q) {
            squarefree[m] = 0;
        }
    }
    let mut fsf = FIArray::eps(N);
    let keys = FIArray::keys(N).collect_vec().into_boxed_slice();
    for q in (2..=SQRT_N).rev() {
        if squarefree[q] == 0 {
            continue;
        }
        for (i, &v) in keys.iter().enumerate() {
            if q > v {
                continue;
            }
            fsf.arr[i] += fsf[v / q];
        }
    }
    for q in SQRT_N + 1..=N {
        if squarefree[q] == 0 {
            continue;
        }
        fsf[N] += fsf[N / q];
    }

    dbg!(fsf[N] - 1);
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
    let mut fsf = FIArray::new(N);
    //let keys = FIArray::keys(N).collect_vec().into_boxed_slice();
    let len = fsf.arr.len();

    fsf.arr[0] = 1;
    let mut f = FenwickTreeUsize::new(0, 0);
    core::mem::swap(&mut fsf.arr, &mut f.0);
    f.construct();

    for q in 2..=SQRT_N {
        if sqf.arr[q - 1] == sqf.arr[q - 2] {
            continue;
        }
        let lim = N / q;
        // sparse_mul_unlimited(q,1)
        let mut prev = 0;
        /*for i in FIArray::keys(lim) {
            let cur = f.sum(sqf.get_index(i));
            if cur != prev {
                f.add(sqf.get_index(i * q), cur - prev);
                prev = cur;
            }
        }*/
        let mut i = 1;
        while i <= lim / i {
            let cur = f.sum(i - 1);
            if cur != prev {
                f.add(sqf.get_index(i * q), cur - prev);
                prev = cur;
            }
            i += 1;
        }
        for j in (1..=lim / i).rev() {
            let cur = f.sum(sqf.get_index(lim / j));
            if cur != prev {
                f.add(len - j, cur - prev);
                prev = cur;
            }
        }
    }
    fsf.arr = f.flatten();
    let mut res = fsf[N] - 1;

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
