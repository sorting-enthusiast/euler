use itertools::Itertools;

use crate::utils::{
    FIArray::{FIArray, FIArrayU128},
    multiplicative_function_summation::{dirichlet_mul_u128, dirichlet_mul_with_buffer_u128},
};
// 1e13: 1107277852610310, 127.1574871s
const N: usize = 1e10 as _;
const SQRT_N: usize = N.isqrt();
const fn icbrt(x: usize) -> usize {
    let mut rt = 1 << (1 + x.ilog2().div_ceil(3));
    let mut x_div_rt2 = (x / rt) / rt;
    while rt > x_div_rt2 {
        rt = ((rt << 1) + x_div_rt2) / 3;
        x_div_rt2 = (x / rt) / rt;
    }
    rt
}

/* #[must_use]
fn sqf(x: u64) -> FIArrayU64 {
    let mut s = FIArrayU64::unit(x);
    let keys = FIArrayU64::keys(x).collect_vec().into_boxed_slice();
    let primes = wheel_sieve(x.isqrt() as _).into_boxed_slice();

    for p in primes {
        for (i, &v) in keys.iter().enumerate().rev() {
            if v < p * p {
                break;
            }
            s.arr[i] -= s[v / (p * p)];
        }
    }
    s
}
 */
#[must_use]
fn sqf(x: usize) -> FIArray {
    let mut Sqf = FIArray::eps(x);
    let xsqrt = Sqf.isqrt;
    let mut d2 = 1;
    for d in 2..=xsqrt.isqrt() {
        d2 += (d << 1) - 1;
        if Sqf.arr[d2 - 1] == 0 {
            continue;
        }
        for m in (d2..=xsqrt).step_by(d2) {
            Sqf.arr[m - 1] = 0;
        }
    }
    let mut sqrts = FIArray::unit(x);
    for v in &mut sqrts.arr {
        *v = v.isqrt();
    }
    for (i, v) in FIArray::keys(x).enumerate().skip(1) {
        if v <= xsqrt {
            Sqf.arr[i] += Sqf.arr[i - 1];
            continue;
        }
        let b = icbrt(v);
        let a = v / (b * b);
        let mut sqf = v + Sqf.arr[a - 1] * b - sqrts[v]; // v.isqrt();
        for i in 2..=a {
            if Sqf.arr[i - 1] != Sqf.arr[i - 2] {
                sqf -= sqrts[v / i]; //(v / i).isqrt();
            }
        }
        for i in 2..=b {
            sqf -= Sqf[v / (i * i)];
        }
        Sqf.arr[i] = sqf;
    }
    Sqf
}
// can optimize using fenwick trees, but the code is fast enough as is
pub fn main() {
    pseudo_euler_transform_based();

    let start = std::time::Instant::now();
    let sqf = sqf(N); //sqf(N as u64);
    dbg!(start.elapsed());
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
fn test1() {
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

// fsf is just the euler transform of sqf
fn pseudo_euler_transform_based() {
    const x: usize = 1 + icbrt(SQRT_N);
    const INVS: [u128; 6] = [0, 60, 30, 20, 15, 12];

    let start = std::time::Instant::now();
    let sqf_ = sqf(N);
    let mut sqf = FIArrayU128::new(N as _);
    for (e, &q) in sqf.arr.iter_mut().zip(&sqf_.arr) {
        *e = q as u128;
    }
    dbg!(start.elapsed());

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

    for q in (2..x).rev() {
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
