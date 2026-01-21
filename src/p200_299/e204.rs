use crate::utils::{FIArray::DirichletFenwick, primes::wheel_sieve};

pub fn main() {
    const N: usize = 1e9 as _;
    let start = std::time::Instant::now();
    let mut smooth = DirichletFenwick::eps(N);
    for p in wheel_sieve(100) {
        smooth.sparse_mul_unlimited(p as _, 1);
    }
    let res = smooth.bit.sum(smooth.get_index(N));
    println!("res = {res}, took {:?}", start.elapsed());
}
