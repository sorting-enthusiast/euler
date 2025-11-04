use std::collections::HashMap;

use itertools::Itertools;

pub fn main() {
    // S(n) = d(n^2) -> S(p^e) = 2e+1
    const N: usize = 1e7 as _;
    const MOD: usize = 1e9 as usize + 87;
    const fn powmod(mut x: usize, mut exp: usize) -> usize {
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
    const fn modinv(x: usize) -> usize {
        powmod(x, MOD - 2)
    }
    let mut lpf = (0..=N).collect_vec();
    for i in 2..=N {
        if lpf[i] == i {
            for m in (2 * i..=N).step_by(i) {
                lpf[m] = i;
            }
        }
    }
    let mut sum = 0;
    let mut prod = 1;
    let mut factors = HashMap::<usize, usize>::new();
    for i in 2..=N {
        let mut k = i;
        while lpf[k] != 1 {
            let p = lpf[k];
            k /= p;
            let mut count = 1;
            while lpf[k] == p {
                count += 1;
                k /= p;
            }
            let prev = factors.entry(p).or_insert(1);

            prod *= modinv(*prev);
            prod %= MOD;
            *prev += 2 * count;
            if *prev >= MOD {
                *prev -= MOD;
            }
            prod *= *prev;
            prod %= MOD;
        }
        sum += prod;
        if sum >= MOD {
            sum -= MOD;
        }
    }
    dbg!(sum);
}
