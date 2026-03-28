use crate::utils::{polymul::NTT, primes::wheel_sieve};
const MOD: i64 = 1_004_535_809; // (497 << 21) + 1
type FFT = NTT<MOD, 702_606_812, 22>;
const N: usize = 20_000 as _;
const LEN: usize = (2 * N + 1).next_power_of_two();
// A-restricted composition - https://en.wikipedia.org/wiki/Composition_(combinatorics)
pub fn main() {
    let start = std::time::Instant::now();
    let primes = wheel_sieve(12 * N as u64);
    let mut x = vec![0; LEN].into_boxed_slice();
    x[0] = 1;
    for i in 1..=N {
        x[i] = (primes[i] - primes[i - 1]) as i32;
    }
    let mut r = vec![0; LEN].into_boxed_slice();
    r[0] = 1;
    let mut exp = N;
    while exp > 1 {
        FFT::ntt(&mut x, false);
        if exp & 1 == 1 {
            FFT::ntt(&mut r, false);
            for i in 0..LEN {
                r[i] = ((i64::from(r[i]) * i64::from(x[i])) % MOD) as i32;
            }
            FFT::ntt(&mut r, true);
            r[N + 1..].fill(0);
            //r = (r * x) % MOD;
        }
        for i in 0..LEN {
            x[i] = ((i64::from(x[i]) * i64::from(x[i])) % MOD) as i32;
        }
        FFT::ntt(&mut x, true);
        x[N + 1..].fill(0);
        //x = (x * x) % MOD;
        exp >>= 1;
    }
    FFT::ntt(&mut x, false);
    FFT::ntt(&mut r, false);
    for i in 0..LEN {
        r[i] = ((i64::from(r[i]) * i64::from(x[i])) % MOD) as i32;
    }
    FFT::ntt(&mut r, true);
    //r[N + 1..].fill(0);
    //(r * x) % MOD
    println!("res = {}, took {:?}", r[N], start.elapsed());
}
