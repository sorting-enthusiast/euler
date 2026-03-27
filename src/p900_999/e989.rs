use crate::utils::{math::sum_geometric_mod, multiplicative_function_summation::mobius_sieve};
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
fn sum_chi_zeta(k: u64, n: usize) -> u64 {
    // coefficient extraction from \frac{kx(1-kx)^2(1+kx)}{(1-x)(1-(kx)^5)}
    const fn sum_chi_x_alt(k: u64, mut n: usize) -> u64 {
        let mut p = [0; 6];
        let mut i = 1;
        let mut kk = k;
        while i < 5 {
            p[i] = kk;
            kk *= k;
            kk %= MOD;
            i += 1;
        }
        p[2] = MOD - p[2];
        p[3] = MOD - p[3];
        while n > 0 {
            if n & 1 == 0 {
                let mut p0p1 = p[0] + p[1];
                if p0p1 >= MOD {
                    p0p1 -= MOD;
                }
                let mut p2p3 = p[2] + p[3];
                if p2p3 >= MOD {
                    p2p3 -= MOD;
                }
                let mut p4p5 = p[4] + p[5];
                if p4p5 >= MOD {
                    p4p5 -= MOD;
                }
                p[1] += p[2];
                if p[1] >= MOD {
                    p[1] -= MOD;
                }
                p[2] = p[3] + p[4];
                if p[2] >= MOD {
                    p[2] -= MOD;
                }
                p[3] = (p[5] + kk * p0p1) % MOD;
                p[4] = (p2p3 * kk) % MOD;
                p[5] = (p4p5 * kk) % MOD;
            } else {
                let mut p1p2 = p[1] + p[2];
                if p1p2 >= MOD {
                    p1p2 -= MOD;
                }
                let mut p3p4 = p[3] + p[4];
                if p3p4 >= MOD {
                    p3p4 -= MOD;
                }
                let p0kk = (p[0] * kk) % MOD;
                p[0] += p[1];
                if p[0] >= MOD {
                    p[0] -= MOD;
                }
                p[1] = p[2] + p[3];
                if p[1] >= MOD {
                    p[1] -= MOD;
                }
                p[2] = (p[4] + p[5] + p0kk) % MOD;
                p[3] = (p1p2 * kk) % MOD;
                p[4] = (p3p4 * kk) % MOD;
                p[5] *= kk;
                p[5] %= MOD;
            }
            kk *= kk;
            kk %= MOD;
            n >>= 1;
        }
        p[0]
    }
    // decomposition of the sum into full periods and one partial period
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
        ret *= if q > 0 {
            sum_geometric_mod::<MOD>(kk, q - 1)
        } else {
            0
        };
        ret %= MOD;
        let k5q = powmod(kk, q as _);
        /* ret *= k5q + MOD - 1;
        ret %= MOD;
        ret *= modinv(kk + MOD - 1);
        ret %= MOD; */

        ret += (k5q * p[r]) % MOD;
        if ret >= MOD {
            ret -= MOD;
        }
        ret
    }
    let mut ret = 0;
    let mut pow_k = 1;
    let sqrt = n.isqrt();
    let k_sqrt = powmod(k, sqrt as u64 + 1);
    let mut pow_k_sqrt = 1;
    for i in 1..=sqrt {
        let chi5 = [0, 1, MOD - 1, MOD - 1, 1][i % 5];
        pow_k *= k;
        pow_k %= MOD;
        pow_k_sqrt *= k_sqrt;
        pow_k_sqrt %= MOD;

        if n / i != sqrt {
            /* let mut coeff = sum_geometric_mod::<MOD>(pow_k, n / i) + MOD
                - sum_geometric_mod::<MOD>(pow_k, sqrt);
            if coeff >= MOD {
                coeff -= MOD;
            } */
            let coeff = (pow_k_sqrt * sum_geometric_mod::<MOD>(pow_k, (n / i) - sqrt - 1)) % MOD;
            ret += chi5 * coeff;
            ret %= MOD;
        }
        ret += sum_chi_x_alt(pow_k, n / i);
        if ret >= MOD {
            ret -= MOD;
        }
    }
    ret
}
const fn powmod(mut x: u64, mut exp: u64) -> u64 {
    exp %= MOD - 1;
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
