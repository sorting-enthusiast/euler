#![warn(clippy::pedantic)]
#![allow(long_running_const_eval)]
#![allow(dead_code)]
#![allow(clippy::too_many_lines)]
#![allow(clippy::cast_possible_truncation)]
#![allow(clippy::inline_always)]
#![allow(clippy::cast_sign_loss)]
#![allow(clippy::cast_precision_loss)]
use crate::utils::{
    FIArray::FIArrayI64,
    multiplicative_function_summation::{
        dirichlet_div_i64, dirichlet_mulmod_i64, mertens, sum_n_i64, totient_sieve, totient_sum,
    },
    prime_sieves::sieve_it,
    primecount,
};

mod e107;
mod e153;
mod e156;
mod e169;
mod e175;
mod e193;
mod e214;
mod e233;
mod e245;
mod e249;
mod e250;
mod e255;
mod e258;
mod e273;
mod e302;
mod e351;
mod e355;
mod e401;
mod e432;
mod e439;
mod e448;
mod e484;
mod e501;
mod e508;
mod e512;
mod e521;
mod e530;
mod e606;
mod e625;
mod e639;
mod e668;
mod e708;
mod e745;
mod e810;
mod e955;
mod longest_collatz_chain;
mod utils;

const MOD: i64 = 1e9 as i64 + 7;
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
const fn modinv(x: i64) -> i64 {
    powmod(x, MOD - 2)
}

pub fn main() {
    /* const N: i64 = 1e10 as _;
    let start = std::time::Instant::now();
    let mob2 = dirichlet_div_i64(&FIArrayI64::eps(N as _), &FIArrayI64::unit(N as _), N as _);
    let end = start.elapsed();
    println!("{end:?}");

    let start = std::time::Instant::now();
    let mob = mertens(N as _);
    let end = start.elapsed();
    println!("{end:?}");

    assert_eq!(mob, mob2);

    let start = std::time::Instant::now();
    let mut res = dirichlet_div_i64(
        &FIArrayI64::id::<MOD>(N as _),
        &FIArrayI64::unit(N as _),
        N as _,
    );
    for v in &mut res.arr {
        *v %= MOD;
        if *v < 0 {
            *v += MOD;
        }
    }
    let end = start.elapsed();
    println!("{end:?}");

    let start = std::time::Instant::now();
    let tot = totient_sum::<MOD>(N as _);
    let end = start.elapsed();
    println!("{end:?}");

    assert_eq!(tot, res); */
    primecount::main();
    //const N: i64 = 2;
    //let kk = sum_f(N)[N];
    /* let mut res = (powmod(2, (N + 1) * (N + 1)) - (N + 1) * (N + 1) - 1) % MOD;
    if res < 0 {
        res += MOD;
    }
    res -= sum_n_i64::<MOD>((N + 1) * (N + 1) - 1);

    res %= MOD;
    if res < 0 {
        res += MOD;
    }
    res += 6 * totient_sum_single::<MOD>(N);
    res %= MOD;
    if res < 0 {
        res += MOD;
    }
    dbg!(res); */
}
const fn sum_squares(x: i64) -> i64 {
    let x = x % MOD;
    (((x * (x + 1)) % MOD * (2 * x + 1)) % MOD * const { modinv(6) }) % MOD
}

fn sum_f(x: i64) -> FIArrayI64 {
    let y = if x > 1023 {
        (1e8 as usize).min((x as f64).powf(2. / 3.) as usize >> 1)
    } else {
        x as usize
    };
    let mut small = totient_sieve(y + 1);
    for i in 2..=y {
        small[i] *= i as i64 % MOD;
        small[i] %= MOD;
        small[i] += small[i - 1];
        small[i] %= MOD;
    }

    let mut ret = FIArrayI64::new(x);

    for (i, v) in FIArrayI64::keys(x).enumerate() {
        if v as usize <= y {
            ret.arr[i] = small[v as usize];
            continue;
        }
        let mut f_v = sum_squares(v) - sum_n_i64::<MOD>(v);
        let vsqrt = v.isqrt();
        for i in 2..=vsqrt {
            f_v -= (i * ret[v / i]) % MOD;
            f_v -= ((ret.arr[i as usize - 1] - ret.arr[i as usize - 2]) * sum_n_i64::<MOD>(v / i))
                % MOD;
            f_v %= MOD;
        }
        f_v += (ret[vsqrt] * sum_n_i64::<MOD>(vsqrt)) % MOD;
        f_v %= MOD;
        if f_v < 0 {
            f_v += MOD;
        }
        ret.arr[i] = f_v;
    }

    ret
}
