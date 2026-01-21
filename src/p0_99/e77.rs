use crate::utils::primes::prime_sieves::sift;
pub fn main() {
    const MAX: usize = 1e2 as _;

    let mut gf = vec![0; MAX + 1];
    let primes = sift(MAX as _);
    gf[0] = 1;
    for &p in &primes {
        let p = p as usize;
        for n in p..=MAX {
            gf[n] += gf[n - p];
        }
    }
    for n in 1..=MAX {
        if gf[n] >= 5000 {
            dbg!(n, gf[n]);
            break;
        }
    }
}
