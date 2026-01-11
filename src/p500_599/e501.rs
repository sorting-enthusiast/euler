use crate::utils::primes::{
    count_signatures::count_signature, primecount::lucy_fenwick, wheel_sieve,
};
// ez: count composites with either 3 distinct prime factors,
// or 2 distinct prime factors, with one of them being of multplicity 3 and the other of 1,
// or p^7

pub fn main() {
    const N: usize = 1e12 as _;
    let start = std::time::Instant::now();
    let primes = wheel_sieve(N.isqrt() as u64);
    let pi = lucy_fenwick(N);
    let sum = count_signature(&[1, 1, 1], N, &pi, &primes)
        + count_signature(&[1, 3], N, &pi, &primes)
        + count_signature(&[3, 1], N, &pi, &primes)
        + count_signature(&[7], N, &pi, &primes);
    println!("res = {sum}, took {:?}", start.elapsed());
}
