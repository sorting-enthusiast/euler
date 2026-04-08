use itertools::Itertools;

fn powmod(mut base: u128, mut exp: usize, m: u128) -> u128 {
    let mut res = 1;
    base %= m;
    while exp > 0 {
        if exp & 1 == 1 {
            res = (res * base) % m;
        }
        base = (base * base) % m;
        exp >>= 1;
    }
    res
}

fn modinv(a: i128, m: i128) -> i128 {
    let (mut t, mut newt) = (0, 1);
    let (mut r, mut newr) = (m, a);
    while newr != 0 {
        let quotient = r / newr;
        (t, newt) = (newt, t - quotient * newt);
        (r, newr) = (newr, r - quotient * newr);
    }
    if t < 0 {
        t += m;
    }
    t
}

fn solve_mod_pk(n: usize, p: u128, k: u32) -> u128 {
    let pk = p.pow(k);

    // The sum truncates when v_p((n)_j) >= k.
    // For k=3, j_max is roughly 3*p.
    let j_max = (p as usize) * (k as usize);

    // 1. Compute A_j = j! * S(n, j) mod p^k using a difference table
    // A_j = Delta^j [i^n]_{i=0}
    let mut diffs = (0..=j_max).map(|i| powmod(i as u128, n, pk)).collect_vec();

    let mut a_vals = Vec::with_capacity(j_max + 1);
    a_vals.push(diffs[0]);

    for j in 1..=j_max {
        for i in 0..(j_max - j + 1) {
            diffs[i] = (diffs[i + 1] + pk - diffs[i]) % pk;
        }
        a_vals.push(diffs[0]);
    }

    // 2. Compute the sum: Sum_{j=0 to j_max} binom(n, j) * A_j * 2^(n-j)
    let mut total_sum = 0;

    for j in 0..=j_max {
        if a_vals[j] == 0 {
            continue;
        }

        // Compute binom(n, j) mod p^k
        // binom(n, j) = (n)_j / j!
        let mut num = 1;
        let mut den = 1;
        let mut v_p_total = 0;

        for i in 0..j {
            let mut val = (n - i) as u128;
            while val > 0 && val.is_multiple_of(p) {
                val /= p;
                v_p_total += 1;
            }
            num = (num * (val % pk)) % pk;

            let mut d_val = i as u128 + 1;
            while d_val > 0 && d_val.is_multiple_of(p) {
                d_val /= p;
                v_p_total -= 1;
            }
            den = (den * (d_val % pk)) % pk;
        }
        if v_p_total >= k {
            continue; // This term is 0 mod p^k
        }

        let binom = (num * modinv(den as _, pk as _) as u128) % pk;
        let binom_final = (binom * p.pow(v_p_total)) % pk;

        let p2 = powmod(2, n - j, pk);
        let term = (binom_final * a_vals[j]) % pk;
        total_sum = (total_sum + term * p2) % pk;
    }

    total_sum
}
const N: u64 = 1e18 as _;
// https://web.archive.org/web/20210120202131/https://min-25.hatenablog.com/entry/2019/10/08/214743
// (xD)^n (1+x)^n = \sum_{k=0}^n {\binom{n}{k}k^n x^k}
// evaluating at x=1 mod p^3 for each of the 3 primes + CRT gives the answer
// evaluation based on the form from min25's blog and noting that (n)_j is a multiple of p^3 for large enough j,
// another cool way to approach the problem is to try to group k's which are equivalent mod p^3,
// which can be done by computing (1 + x)^n mod x^{p^3} - 1
pub fn main() {
    let primes = [83, 89, 97];
    let k = 3;

    let mut results = vec![];
    let mut moduli = vec![];

    for &p in &primes {
        let res = solve_mod_pk(N as _, p, k);
        results.push(res);
        moduli.push(p.pow(k));
    }

    let mut m_total = 1;
    for m in &moduli {
        m_total *= m;
    }

    let mut res = 0;
    for i in 0..results.len() {
        let mi = m_total / moduli[i];
        let inv = modinv(mi as _, moduli[i] as _) as u128;
        let term = (results[i] * mi) % m_total;
        res = (res + (term * inv) % m_total) % m_total;
    }

    println!("res = {res}");
}
