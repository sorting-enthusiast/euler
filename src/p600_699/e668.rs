use crate::utils::primes::primecount::lucy_fastdivide_alt;

// n is not square root smooth iff x=py for some prime p and p>=y
// for given y <= N.isqrt(), subtract # of primes below N/y and above y: pi(x/y) - pi(y-1)
pub fn main() {
    const N: u64 = 1e10 as _;
    let start = std::time::Instant::now();
    let pis = lucy_fastdivide_alt(N);
    let mut count = N - pis[N];
    for i in 2..=N.isqrt() {
        count -= pis[N / i] - pis.arr[i as usize - 2];
    }
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");
}
