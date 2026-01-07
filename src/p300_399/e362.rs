use itertools::Itertools;

use crate::utils::{
    FIArray::FIArray,
    multiplicative_function_summation::{dirichlet_mul_usize, dirichlet_mul_with_buffer_usize},
    primes::log_zeta::dirichlet_mul_zero_prefix,
};
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
    const N: usize = 1e10 as _;
    const SQRT_N: usize = N.isqrt();
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
        for (i, &v) in keys[..=len - q].iter().enumerate().skip(q - 1) {
            fsf.arr[i] += fsf[v / q];
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
    const N: usize = 1e10 as _;
    const SQRT_N: usize = N.isqrt();
    let start = std::time::Instant::now();
    let sqf = sqf(N);
    dbg!(start.elapsed());

    const INVS: [usize; 6] = [0, 60, 30, 20, 15, 12];
    let len = sqf.arr.len();

    let mut a_vals = sqf.clone();
    for i in (1..len).rev() {
        a_vals[i] -= a_vals[i - 1];
    }
    const x: usize = 1 + icbrt(SQRT_N);
    for i in (x..=SQRT_N).rev() {
        let v = a_vals.arr[i - 1];
        if v == 0 {
            continue;
        }
        assert_eq!(a_vals.arr[i - 1], 1);
        a_vals.arr[i - 1] = INVS[1];
        let mut e = 1;
        let mut pi = i;
        while pi <= N / i {
            e += 1;
            pi *= i;
            a_vals[pi] += INVS[e];
        }
    }
    println!("bello");
    let mut tmp = FIArray::new(N);

    let mut v = FIArray::new(N);
    for i in x..=len {
        v.arr[i - 1] += v.arr[i - 2] + a_vals.arr[i - 1];
    }
    let v = v;
    let mut ret = v.clone();
    for e in &mut ret.arr {
        *e = (*e + INVS[1]) * 120;
    }
    let mut r = dirichlet_mul_zero_prefix(&v, &v, N, x, x);
    for i in x..=len {
        ret.arr[i - 1] += r.arr[i - 1] * 60;
    }
    println!("bello");

    dirichlet_mul_with_buffer_usize(&r, &v, N, &mut tmp);
    core::mem::swap(&mut r.arr, &mut tmp.arr);

    for i in x..=len {
        ret.arr[i - 1] += r.arr[i - 1] * 20;
    }
    println!("bello");

    dirichlet_mul_with_buffer_usize(&r, &v, N, &mut tmp);
    core::mem::swap(&mut r.arr, &mut tmp.arr);
    for i in x..=len {
        ret.arr[i - 1] += r.arr[i - 1] * 5;
    }
    println!("bello");

    dirichlet_mul_with_buffer_usize(&r, &v, N, &mut tmp);
    core::mem::swap(&mut r.arr, &mut tmp.arr);
    for i in x..=len {
        ret.arr[i - 1] += r.arr[i - 1];
    }
    println!("bello");

    for e in &mut r.arr {
        assert_eq!(*e % 120, 0);
        *e /= 120;
    }
    println!("bello");

    for i in (2..x).rev() {
        let ax = a_vals.arr[i - 1];
        if ax == 0 {
            continue;
        }
        for j in i..=len {
            r.arr[j - 1] += r[{
                if j <= SQRT_N {
                    j - 1
                } else {
                    N / (len - j + 1)
                }
            } / i];
        }
    }
    let res = r[N];
    println!("res = {res}, took {:?}", start.elapsed());
}
