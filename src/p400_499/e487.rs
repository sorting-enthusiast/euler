use crate::utils::primes::wheel_sieve;

pub fn main() {
    const N: usize = 1e12 as _;
    const K: usize = 1e4 as _;
    const LOW: usize = 2e9 as _;
    const HIGH: usize = LOW + 2000;

    let start = std::time::Instant::now();
    let primes = wheel_sieve(HIGH.isqrt() as _);
    let mut is_prime = vec![true; HIGH - LOW].into_boxed_slice();
    for p in primes {
        for m in (LOW.next_multiple_of(p as _) - LOW..HIGH - LOW).step_by(p as _) {
            is_prime[m] = false;
        }
    }
    let mut pts = vec![0; K + 3].into_boxed_slice();
    let mut sum = 0;
    for (i, v) in is_prime.iter().enumerate() {
        if !v {
            continue;
        }
        let p = i + LOW;
        for i in 1..=K + 2 {
            pts[i] = pts[i - 1] + powmod(i, K, p); // can optimize using linear sieve
            if pts[i] >= p {
                pts[i] -= p;
            }
        }
        for i in 1..=K + 2 {
            pts[i] += pts[i - 1];
            if pts[i] >= p {
                pts[i] -= p;
            }
        }
        sum += interpolate_evaluate(&pts, N, p);
    }
    println!("res = {sum}, took {:?}", start.elapsed());
}

fn interpolate_evaluate(y: &[usize], x_eval: usize, p: usize) -> usize {
    let n = y.len() - 1;
    let x_mod = x_eval % p;

    if x_mod <= n {
        return y[x_mod];
    }

    let mut fact = vec![1; n + 1];
    let mut inv_fact = vec![1; n + 1];

    for i in 1..=n {
        fact[i] = (fact[i - 1] * i) % p;
    }

    inv_fact[n] = modinv(fact[n], p);
    for i in (1..=n).rev() {
        inv_fact[i - 1] = (inv_fact[i] * i) % p;
    }

    let mut pref = vec![1; n + 1];
    let mut suff = vec![1; n + 1];

    pref[0] = 1;
    for i in 1..=n {
        let term = (x_mod + p - (i - 1)) % p;
        pref[i] = (pref[i - 1] * term) % p;
    }

    suff[n] = 1;
    for i in (0..n).rev() {
        let term = (x_mod + p - (i + 1)) % p;
        suff[i] = (suff[i + 1] * term) % p;
    }

    let mut result = 0;

    for i in 0..=n {
        // product of (x - j) for all j != i
        let num = (pref[i] * suff[i]) % p;

        // inverse of (i! * (n-i)!)
        let den = (inv_fact[i] * inv_fact[n - i]) % p;

        let sign = if (n - i) & 1 == 1 { p - 1 } else { 1 };

        // Term = y[i] * numerator * denominator * sign
        let mut term = y[i];
        term = (term * num) % p;
        term = (term * den) % p;
        term = (term * sign) % p;

        result = (result + term) % p;
    }

    result
}

fn powmod(mut x: usize, mut exp: usize, p: usize) -> usize {
    if exp == 0 {
        return 1;
    }
    let mut r = 1;
    x %= p;
    while exp > 1 {
        if exp & 1 == 1 {
            r = (r * x) % p;
        }
        x = (x * x) % p;
        exp >>= 1;
    }
    (r * x) % p
}

fn modinv(n: usize, p: usize) -> usize {
    powmod(n, p - 2, p)
}
