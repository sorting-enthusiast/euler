use crate::utils::multiplicative_function_summation::mobius_sieve;
const N: usize = 1e14 as _;
const SQRT_N: usize = N.isqrt();

const MOD: u64 = 1e9 as u64 + 9;
const SQRT_5: u64 = 383_008_016;
const PHI_PLUS: u64 = (((MOD + 1) >> 1) * (1 + SQRT_5)) % MOD;
const PHI_MINUS: u64 = (((MOD + 1) >> 1) * (1 + MOD - SQRT_5)) % MOD;
// the DGF of G(n) is \frac{\zeta(s)L(s;\chi_5)}{\zeta(2s)}
// golden ratio has exact repr mod 10^9 + 9, use that and dirichlet hyperbola to compute the sum
pub fn main() {
    let start = std::time::Instant::now();
    let mu = mobius_sieve(SQRT_N + 1);
    let mut res = 0;
    for a in 1..=SQRT_N {
        if mu[a] == 0 {
            continue;
        }
        let aa = a * a;
        let naa = N / aa;
        let aa = aa as u64;
        let mut coeff = sum_chi_zeta(powmod(PHI_PLUS, aa), naa) + MOD
            - sum_chi_zeta(powmod(PHI_MINUS, aa), naa);
        if coeff >= MOD {
            coeff -= MOD;
        }
        res += if mu[a] == 1 { coeff } else { MOD - coeff };
        if res >= MOD {
            res -= MOD;
        }
    }
    res *= const { modinv(SQRT_5) };
    res %= MOD;
    println!("res = {res}, took {:?}", start.elapsed());
}
// assumes k != 1
const fn sum_geometric(k: u64, n: u64) -> u64 {
    ((modinv(k - 1) * (powmod(k, n) + MOD - 1)) % MOD * k) % MOD
}
fn sum_chi_zeta(k: u64, n: usize) -> u64 {
    const fn sum_chi_x(k: u64, n: usize) -> u64 {
        let (q, r) = (n / 5, n % 5);
        if k == 1 {
            return [0, 1, 0, MOD - 1, 0][r];
        }
        let mut p = [0; 5];
        p[1] = k;
        let mut ret = k;
        let mut kk = (k * k) % MOD;
        ret += MOD - kk;
        if ret >= MOD {
            ret -= MOD;
        }
        p[2] = ret;
        kk *= k;
        kk %= MOD;
        ret += MOD - kk;
        if ret >= MOD {
            ret -= MOD;
        }
        p[3] = ret;
        kk *= k;
        kk %= MOD;
        ret += kk;
        if ret >= MOD {
            ret -= MOD;
        }
        p[4] = ret;

        kk *= k;
        kk %= MOD;
        let k5q = powmod(kk, q as _);
        ret *= MOD + 1 - k5q;
        ret %= MOD;
        ret *= modinv(MOD + 1 - kk);
        ret %= MOD;

        ret += (k5q * p[r]) % MOD;
        if ret >= MOD {
            ret -= MOD;
        }
        ret
    }
    let mut ret = 0;
    let mut pow_k = 1;
    let sqrt = n.isqrt();
    for i in 1..=sqrt {
        let chi5 = [0, 1, MOD - 1, MOD - 1, 1][i % 5];
        pow_k *= k;
        pow_k %= MOD;
        let mut coeff = sum_geometric(pow_k, (n / i) as _) + MOD - sum_geometric(pow_k, sqrt as _);
        if coeff >= MOD {
            coeff -= MOD;
        }
        ret += chi5 * coeff;
        ret %= MOD;
        ret += sum_chi_x(pow_k, n / i);
        if ret >= MOD {
            ret -= MOD;
        }
    }
    ret
}
fn dfs2(omega: u64, acc: u64, lim: u64, primes: &[u64]) -> u64 {
    let mut sum = (omega * fib(acc)) % MOD;
    let mut new_omega = omega;
    new_omega <<= 1;
    if new_omega >= MOD {
        new_omega -= MOD;
    }
    for (i, &p) in primes.iter().enumerate() {
        if p > lim {
            break;
        }
        let mut new_lim = lim;
        let mut pp = 1;
        loop {
            pp *= p;
            new_lim /= p;
            sum += dfs2(new_omega, acc * pp, new_lim, &primes[i + 1..]);
            if sum >= MOD {
                sum -= MOD;
            }
            if p > new_lim {
                break;
            }
        }
    }
    sum
}

fn dfs(omega: u64, acc: u64, lim: u64, primes: &[u64]) -> u64 {
    let mut sum = (omega * (powmod(PHI_PLUS, acc) + MOD - powmod(PHI_MINUS, acc))) % MOD;
    let mut new_omega = omega;
    new_omega <<= 1;
    if new_omega >= MOD {
        new_omega -= MOD;
    }
    for (i, &p) in primes.iter().enumerate() {
        if p > lim {
            break;
        }
        let mut new_lim = lim;
        let mut pp = 1;
        loop {
            pp *= p;
            new_lim /= p;
            sum += dfs(new_omega, acc * pp, new_lim, &primes[i + 1..]);
            if sum >= MOD {
                sum -= MOD;
            }
            if p > new_lim {
                break;
            }
        }
    }
    sum
}
const fn fib(mut n: u64) -> u64 {
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
const fn modinv(x: u64) -> u64 {
    powmod(x, MOD - 2)
}
