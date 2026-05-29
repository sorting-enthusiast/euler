use crate::utils::primes::{primality::powmod, wheel_sieve};

pub fn main() {
    assert_eq!(480, wheel_sieve(1e2 as _)[2..].iter().map(S).sum::<u64>());

    let start = std::time::Instant::now();
    let res = wheel_sieve(1e8 as _)[2..].iter().map(S).sum::<u64>();
    println!("res = {res}, took {:?}", start.elapsed());
}

const fn S(&p: &u64) -> u64 {
    // (p-1)! = -1
    // (p-2)! = -1/-1 = 1
    // (p-3)! = 1/-2
    // (p-4)! = 1/(-2*-3) = 1/6
    // (p-5)! = 1/(6*-4) = -1/24
    // -1/2 + 1/6 - 1/24 = (-12 + 4 - 1) / 24 = -9/24 = -3/8
    {
        // via Lucy_Hedgehog's forum post
        let k = (3 * p) & 7;
        (k * p - 3) >> 3
    }
    //((p - 3) * modinv(8, p)) % p
}

const fn modinv(x: u64, p: u64) -> u64 {
    powmod(x, p - 2, p)
}
