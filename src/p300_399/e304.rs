use crate::utils::primes::wheel_sieve;

const MOD: u128 = 1_234_567_891_011;
const LOW: u128 = 1e14 as _;
const HIGH: u128 = LOW + 4e6 as u128;
pub fn main() {
    let start = std::time::Instant::now();
    let primes = wheel_sieve(HIGH.isqrt() as _);
    let mut is_prime = vec![true; (HIGH - LOW) as usize].into_boxed_slice();
    for p in primes {
        for m in (LOW.next_multiple_of(p as _) - LOW..HIGH - LOW).step_by(p as _) {
            is_prime[m as usize] = false;
        }
    }
    let mut pi = 0;
    let mut sum = 0;
    let mut f = fib(LOW - 1);
    let mut f1 = fib(LOW);
    for v in is_prime {
        (f, f1) = (f1, (f + f1) % MOD);
        if !v {
            continue;
        }
        sum += f;
        if sum >= MOD {
            sum -= MOD;
        }
        pi += 1;
        if pi >= 100_000 {
            println!("res = {sum}, took {:?}", start.elapsed());
            break;
        }
    }
}
const fn fib(mut n: u128) -> u128 {
    if n < 2 {
        return [0, 1][n as usize & 1];
    }
    let mut c = 3;
    let (mut a, mut b) = if n & 1 == 0 { (0, 1) } else { (1, MOD - 1) };
    n >>= 1;
    while n > 1 {
        if n & 1 == 0 {
            b = (a + b * c) % MOD;
        } else {
            a = (b + a * c) % MOD;
        }
        c = (c * c - 2) % MOD;
        n >>= 1;
    }
    (b + a * c) % MOD
}
