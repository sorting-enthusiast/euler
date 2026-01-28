use crate::utils::{
    FIArray::FIArrayI64, powerful_numbers::PowerfulExtAlt, primes::prime_sieves::sift,
};

const N: usize = 1e12 as _;
const MOD: i64 = 1e9 as i64 + 7;
const KMAX: usize = 50;

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

fn sum_k(x: i64, k: usize) -> i64 {
    // --- FACTORIALS ---
    const fn factorials() -> ([i64; KMAX + 3], [i64; KMAX + 3]) {
        let mut fact = [1i64; KMAX + 3];
        let mut invfact = [1i64; KMAX + 3];

        let mut i = 1;
        while i <= KMAX + 2 {
            fact[i] = (fact[i - 1] * i as i64) % MOD;
            i += 1;
        }

        invfact[KMAX + 2] = modinv(fact[KMAX + 2]);
        let mut i = KMAX + 1;
        while i > 0 {
            invfact[i] = (invfact[i + 1] * (i + 1) as i64) % MOD;
            i -= 1;
        }

        (fact, invfact)
    }

    // --- BINOMIAL ---
    const fn binom(n: usize, k: usize, fact: &[i64; KMAX + 3], invfact: &[i64; KMAX + 3]) -> i64 {
        ((fact[n] * invfact[k]) % MOD * invfact[n - k]) % MOD
    }

    // --- BERNOULLI NUMBERS MOD MOD ---
    const fn bernoulli() -> [i64; KMAX + 2] {
        let mut a = [0i64; KMAX + 2];
        let mut b = [0i64; KMAX + 2];

        let mut m = 0;
        while m <= KMAX + 1 {
            a[m] = modinv(m as i64 + 1);
            let mut j = m;
            while j > 0 {
                j -= 1;
                let tmp = (((j + 1) as i64) * (a[j] - a[j + 1] + MOD) % MOD) % MOD;
                a[j] = tmp;
            }
            b[m] = a[0];
            m += 1;
        }
        //b[1] = (MOD - b[1]) % MOD;
        b
    }

    // --- FAULHABER COEFFICIENTS ---
    const fn faulhaber_coeffs() -> [[i64; KMAX + 2]; KMAX + 1] {
        let bern = bernoulli();
        let (fact, invfact) = factorials();
        let mut coeffs = [[0i64; KMAX + 2]; KMAX + 1];

        let mut k = 1;
        while k <= KMAX {
            let inv_kp1 = modinv(k as i64 + 1);
            let mut j = 0;
            while j <= k {
                let b = bern[j];
                if b != 0 {
                    let bin = binom(k + 1, j, &fact, &invfact);
                    let deg = k + 1 - j;
                    let term = bin * b % MOD * inv_kp1 % MOD;
                    coeffs[k][deg] = (coeffs[k][deg] + term + MOD) % MOD;
                }
                j += 1;
            }
            k += 1;
        }
        coeffs
    }

    static FAULHABER_COEFFS: [[i64; KMAX + 2]; KMAX + 1] = faulhaber_coeffs();

    let x = x % MOD;
    let mut xp = x;
    let mut sum = 0;
    for c in &FAULHABER_COEFFS[k][1..=k + 1] {
        sum += (c * xp) % MOD;
        xp = (x * xp) % MOD;
    }
    sum % MOD
}

// Kinda dumb, powerful number trick with g = N^k and h(p^e) = p^k (1-p^k)
// Precalculate faulhaber coefficients at compile time
// At runtime, iterate over every k and just sum using powerful number trick with g and h
// Accelerated by caching partial sums of g in a FIArray,
// as there are few possible values for floor(N/i) - reduces time from 12 seconds to 2
// O(sqrt(n) + k^2) space, O(k^2*sqrt(n)) time at runtime
pub fn main() {
    let start = std::time::Instant::now();
    let ps = sift(N.isqrt() as u64 + 1);

    let mut sum = 0;
    let mut cache = FIArrayI64::new(N);
    for k in 1..=KMAX {
        cache.arr.fill(0);
        let h = |p, e| {
            unsafe { core::hint::assert_unchecked(e > 1) };
            let pk = powmod(p, k as _);
            (pk * (MOD + 1 - pk)) % MOD
        };
        for (n, hn) in PowerfulExtAlt::<_, MOD>::new(N as _, h, &ps).filter(|&(_, hn)| hn != 0) {
            let i = cache.get_index(N / n as usize);
            if cache.arr[i] == 0 {
                cache.arr[i] = sum_k(N as i64 / n, k);
            }
            sum += (hn * cache.arr[i]) % MOD;
            sum %= MOD;
        }
        println!("k = {k}, took {:?}", start.elapsed());
    }

    let end = start.elapsed();
    println!("res = {sum}, took {end:?}");
}
