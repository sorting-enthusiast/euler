use crate::utils::prime_sieves::sieve_it;

const MOD: i64 = 1e9 as i64 + 9;
const N: usize = 1e8 as _;

const fn powmod(x: i64, mut exp: i64) -> i64 {
    if exp == 0 {
        return 1;
    }
    let mut r = 1i64;
    let mut x = x as i64;
    x %= MOD as i64;
    while exp > 1 {
        if exp & 1 == 1 {
            r = (r * x) % MOD as i64;
        }
        x = (x * x) % MOD as i64;
        exp >>= 1;
    }
    ((r * x) % MOD as i64) as i64
}
const fn modinv(x: i64) -> i64 {
    powmod(x, MOD - 2)
}
const fn mulmod(x: i64, y: i64) -> i64 {
    (x % MOD * y % MOD) % MOD
}
pub fn main() {
    let start = std::time::Instant::now();
    let mut res = 0;
    let mut primes = sieve_it();
    let mut q2 = 2;
    let mut q3 = 3;
    const THREE_HALVES: i64 = (3 * modinv(2)) % MOD;
    let mut p = primes.next().unwrap();
    for n in 2..=N {
        if n.trailing_zeros() >= 23 {
            dbg!(n, res);
        }
        let invn = modinv(n as _);
        q2 *= invn;
        q2 %= MOD;
        q2 *= 4 * n as i64 - 2;
        q2 %= MOD;

        q3 *= THREE_HALVES;
        q3 %= MOD;
        q3 *= invn;
        q3 %= MOD;
        q3 *= ((3 * n as i64 - 1) * (3 * n as i64 - 2)) % MOD;
        q3 %= MOD;
        q3 *= modinv(2 * n as i64 - 1);
        q3 %= MOD;

        if n == p {
            res += (5 + invn * (q2 + q3 - 5)) % MOD;
            if res >= MOD {
                res -= MOD;
            }
            p = primes.next().unwrap();
            if p > N {
                break;
            }
        }
    }
    println!("{}, {:?}", res - 5, start.elapsed());
}
