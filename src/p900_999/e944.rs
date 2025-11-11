use crate::utils::multiplicative_function_summation::sum_n_i64;
const N: i64 = 1e14 as _;
const SQRT_N: i64 = N.isqrt();
const MOD: i64 = 1_234_567_891;
const fn powmod(mut x: i64, mut exp: i64) -> i64 {
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

pub fn main() {
    const TWO_N_1: i64 = powmod(2, N - 1);
    const INV2: i64 = (MOD + 1) >> 1;

    let mut b = 0;
    for d in 1..=SQRT_N {
        b += (d * powmod(INV2, N / d)) % MOD;
        if b >= MOD {
            b -= MOD;
        }
    }
    for k in 1..N / SQRT_N {
        b += (powmod(INV2, k) * (sum_n_i64::<MOD>(N / k) + MOD - sum_n_i64::<MOD>(N / (k + 1)))
            % MOD)
            % MOD;
        if b >= MOD {
            b -= MOD;
        }
    }
    b *= (TWO_N_1 << 1) % MOD;
    b %= MOD;
    println!(
        "{}",
        (const { (sum_n_i64::<MOD>(N) * TWO_N_1) % MOD } + MOD - b) % MOD
    );
}
