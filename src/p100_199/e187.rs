use crate::utils::{prime_sieves::sift, primecount::lucy};
const N: usize = 1e8 as _;
pub fn main() {
    let start = std::time::Instant::now();
    let mut cnt = 0;
    let pi = lucy(N);
    for (i, p) in sift(N.isqrt() as u64).into_iter().enumerate() {
        let p = p as usize;
        cnt += pi[N / p] - i;
    }
    println!("res = {cnt}, took {:?}", start.elapsed());
}
