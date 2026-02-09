use crate::utils::primes::wheel_sieve;

const MOD: u64 = 1e9 as u64 + 9;
const N: u64 = 1e8 as _;
const fn powmod(mut x: u64, mut exp: u64) -> u64 {
    if exp == 0 {
        return 1;
    }
    let mut r = 1;
    x %= MOD;
    while exp > 1 {
        if exp & 1 == 1 {
            r = (r * x) % MOD;
        }
        x = (x * x) % MOD;
        exp >>= 1;
    }
    (r * x) % MOD
}
// https://oeis.org/A034676
pub fn main() {
    let mut res = 1;
    let primes = wheel_sieve(N);
    for p in primes {
        let mut cnt = N / p;
        let mut e = cnt;
        while p <= cnt {
            cnt /= p;
            e += cnt;
        }
        res *= powmod(p, e << 1) + 1;
        res %= MOD;
    }
    dbg!(res);
}
