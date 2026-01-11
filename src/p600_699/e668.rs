use crate::utils::primes::primecount::lucy_fenwick;

// n is not square root smooth iff x=py for some prime p and p>=y
// for given y <= N.isqrt(), subtract # of primes below N/y and above y: pi(x/y) - pi(y-1)
pub fn main() {
    const N: usize = 1e10 as _;
    let start = std::time::Instant::now();
    let pis = lucy_fenwick(N);
    let len = pis.arr.len();
    let mut count = N - pis.arr[len - 1];
    for i in 2..=pis.isqrt {
        count -= pis.arr[len - i] - pis.arr[i - 2];
    }
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");
}
