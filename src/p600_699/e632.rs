use itertools::Itertools;

const N: usize = 1e16 as _;
const SQRT: usize = N.isqrt();
const B: usize = icbrt(N);
const A: usize = N / (B * B);
const MOD: usize = 1e9 as usize + 7;

const fn icbrt(x: usize) -> usize {
    let mut rt = 1 << (1 + x.ilog2().div_ceil(3));
    let mut x_div_rt2 = (x / rt) / rt;
    while rt > x_div_rt2 {
        rt = ((rt << 1) + x_div_rt2) / 3;
        x_div_rt2 = (x / rt) / rt;
    }
    rt
}
fn sqf_sieve(n: usize) -> Vec<usize> {
    unsafe { core::hint::assert_unchecked(n >= 1) };
    let mut sqf = vec![1; n];
    sqf[0] = 0;
    if n < 2 {
        return sqf;
    }
    let sqrtn = n.isqrt();
    let mut d2 = 1;
    for d in 2..=sqrtn {
        d2 += (d << 1) - 1;
        for m in (d2..n).step_by(d2) {
            sqf[m] = 0;
        }
    }
    sqf
}
fn omega_sieve(n: usize) -> Vec<u8> {
    unsafe { core::hint::assert_unchecked(n >= 1) };
    let mut omega = vec![0; n];
    if n < 2 {
        return omega;
    }
    for p in 2..n {
        if omega[p] == 0 {
            for m in (p..n).step_by(p) {
                omega[m] += 1;
            }
        }
    }
    omega
}
pub fn main() {
    let start = std::time::Instant::now();
    let mut small_sqf = sqf_sieve(A + 1);
    for i in 2..=A {
        small_sqf[i] += small_sqf[i - 1];
    }
    let sqrts = (2..=A)
        .map(|i| (N / i).isqrt())
        .collect_vec()
        .into_boxed_slice(); // precompute once, used very often
    let mut large_sqf = vec![0; B - 1].into_boxed_slice(); // indexed by denominator
    for d in (1..B).rev() {
        let v = N / (d * d);
        let b = icbrt(v);
        let a = v / (b * b);

        let mut sqf = v + small_sqf[a] * b - SQRT / d;
        for i in 2..=a {
            sqf -= (small_sqf[i] - small_sqf[i - 1]) * sqrts[i - 2] / d; //(v / i).isqrt();
        }
        for i in 2..=b {
            sqf -= if d * i < B {
                large_sqf[(i * d) - 1]
            } else {
                small_sqf[v / (i * i)]
            };
        }
        large_sqf[d - 1] = sqf;
    }
    println!("Counted squarefree numbers: {:?}", start.elapsed());
    let omega = omega_sieve(SQRT + 1);
    println!("Counted distinct prime factors: {:?}", start.elapsed());

    let mut C = [0; N.ilog2() as usize + 1];
    let mut ii = 0;
    for i in 1..B {
        ii += (i << 1) - 1;
        C[omega[i] as usize] += large_sqf[i - 1];
    }
    for i in B..=SQRT {
        ii += (i << 1) - 1;
        C[omega[i] as usize] += small_sqf[N / ii];
    }
    let mut res = 1;
    for (k, c) in C.into_iter().enumerate() {
        if c == 0 {
            break;
        }
        println!("C_{k} = {c}");
        res *= c % MOD;
        res %= MOD;
    }
    println!("res = {res}, took {:?}", start.elapsed());
}
