use crate::utils::{
    FIArray::{DirichletFenwick, FIArray},
    multiplicative_function_summation::{
        dirichlet_mul_single_usize, dirichlet_mul_usize, sum_n_usize,
    },
    primes::wheel_sieve,
};

const N: usize = 1e9 as _;

// https://en.wikipedia.org/wiki/Coprime_integers#Generating_all_coprime_pairs
fn calkin_wilf(m: usize, n: usize, sum_divisors: &FIArray) -> usize {
    let norm = m * m + n * n;
    if norm > N {
        return 0;
    }
    (m + n) * sum_divisors[N / norm]
        + calkin_wilf(m + n, m, sum_divisors)
        + calkin_wilf(m + n, n, sum_divisors)
}
fn calkin_wilf_iter(m: usize, n: usize, sum_divisors: &FIArray) -> usize {
    let mut res = 0;
    let mut stack = Vec::with_capacity(16);
    stack.push((m, n));
    while let Some((m, n)) = stack.pop() {
        let norm = m * m + n * n;
        if norm <= N {
            res += (m + n) * sum_divisors[N / norm];
            stack.push((m + n, m));
            stack.push((m + n, n));
        }
    }
    res
}
pub fn main() {
    min_25();
    let start = std::time::Instant::now();
    let sum_divisors = dirichlet_mul_usize(&FIArray::unit(N), &FIArray::id::<0>(N), N);
    let mut res = sum_divisors[N / 2] + calkin_wilf(2, 1, &sum_divisors);
    res <<= 1;
    res += sum_divisors[N];
    println!("res = {res}, took {:?}", start.elapsed());
}
fn min_25() {
    let start = std::time::Instant::now();

    let s = dirichlet_mul_usize(&FIArray::unit(N), &FIArray::id::<0>(N), N);
    let mut res = s[N];
    //dbg!(start.elapsed());

    let mut s = DirichletFenwick::from(s);
    for p in wheel_sieve(s.isqrt as u64) {
        let p = p as usize;
        s.sparse_mul_at_most_one(p * p, p);
    }
    let s = FIArray::from(s);
    //dbg!(start.elapsed());

    let mut C1 = FIArray::new(N);
    let c1 = |v: usize| {
        let s = (v / 2).isqrt();
        (s + 1..=v.isqrt())
            .map(move |x| {
                let u = (v - x * x).isqrt();
                x * u + sum_n_usize::<0>(u)
            })
            .sum::<usize>()
            + s * sum_n_usize::<0>(s)
    };
    for (i, v) in FIArray::keys(N).enumerate() {
        C1.arr[i] = c1(v);
    }
    //dbg!(start.elapsed());

    res += 2 * dirichlet_mul_single_usize(&s, &C1, N);
    println!("res = {res}, took {:?}", start.elapsed());
}
