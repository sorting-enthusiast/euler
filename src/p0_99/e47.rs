use crate::utils::primes::wheel_sieve;

pub fn main() {
    const LIM: usize = 2e5 as _;
    let primes = wheel_sieve(LIM as _);
    let mut fcs = vec![0u8; LIM + 1];
    for p in primes {
        let p = p as usize;
        fcs[p] = 1;

        for m in (2 * p..=LIM).step_by(p) {
            fcs[m] += 1;
        }
    }
    for (i, w) in fcs[1..].windows(4).enumerate() {
        if w == [4, 4, 4, 4] {
            dbg!(i + 1);
            break;
        }
    }
}
