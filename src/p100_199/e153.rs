use crate::utils::{FIArray::FIArray, multiplicative_function_summation::dirichlet_mul_usize};

const N: usize = 1e8 as _;

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

pub fn main() {
    let start = std::time::Instant::now();
    let sum_divisors =
        dirichlet_mul_usize(&FIArray::unit(N), &FIArray::id::<{ usize::MAX >> 1 }>(N), N);
    let mut res = sum_divisors[N / 2] + calkin_wilf(2, 1, &sum_divisors);
    res <<= 1;
    res += sum_divisors[N];
    println!("res = {res}, took {:?}", start.elapsed());
}
