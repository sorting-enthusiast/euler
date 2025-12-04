use crate::utils::primes::primality::is_prime;

pub fn main() {
    const N: u64 = 5e7 as _;
    let mut count = 0u64;
    let mut p = 1;
    for n in 2..=N {
        p += (n << 2) - 2;

        if is_prime(p) {
            count += 1;
            if count.trailing_zeros() >= 16 {
                println!("{count}: {n}, {p}");
            }
        }
    }
    dbg!(count);
}
