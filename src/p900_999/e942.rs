use std::time::Instant;

const MOD: usize = 1_000_000_007;

// Legendre symbol (a|p) for odd prime p: a^((p-1)/2) mod p (returns 1 or p-1 since p!|a)
const fn legendre_symbol(a: usize, p: usize) -> i32 {
    let mut base = a % p;
    let mut exp = (p - 1) / 2;
    let mut ls = 1;
    while exp > 0 {
        if exp & 1 == 1 {
            ls = (ls * base) % p;
        }
        base = (base * base) % p;
        exp >>= 1;
    }

    if ls == 1 { 1 } else { -1 }
}

// https://arxiv.org/pdf/1309.2831
// https://en.wikipedia.org/wiki/Quadratic_Gauss_sum
pub fn main() {
    const Q: usize = 74_207_281;
    const CHI_Q_MINUS_2: i32 = legendre_symbol(Q - 2, Q); // = (q-2 | q)

    // Two-level table for powers of 2 modulo MOD:
    // 2^e = big[e / B] * small[e % B]  (mod MOD)
    // big[i] = 2^(i*B) mod MOD for 0 <= i <= floor(q / B)
    // small[j] = 2^j mod MOD for 0 <= j < B

    const B: usize = Q.isqrt();
    let start = Instant::now();

    let mut small = vec![0; B];
    let mut big = vec![0; B + 1];
    small[0] = 1;
    big[0] = 1;

    for j in 1..B {
        small[j] = (small[j - 1] << 1) % MOD;
    }
    let step = (small[B - 1] << 1) % MOD;

    for i in 1..=B {
        big[i] = (big[i - 1] * step) % MOD;
    }

    let pow2_at = |e| (big[e / B] * small[e % B]) % MOD;

    let p = pow2_at(Q) - 1;

    // 2 is an order q root of unity modulo p, and p = 1 mod 4, therefore
    // sqrt(q) = sum_{a=1}^{q-1} (a|q)2^a = 2 * sum_{t=1}^{(q-1)/2} 2^{t^2 mod q} - (2^q - 2)
    // proven by splitting sum over quadratic residues and non-residues modulo q
    // notice that 2^q - 2 is exactly -1 mod p, so we can compute 1 + sum_{t=1}^{(q-1)/2} 2^{t^2 + 1 mod q} instead
    let root = ((1..=(Q - 1) / 2)
        .map(|t| (1 + t * t) % Q)
        .map(pow2_at)
        .sum::<usize>()
        + 1)
        % MOD;

    // Decide minimal root: check (q-2 / q)
    let minimal_root = if CHI_Q_MINUS_2 == -1 {
        root
    } else {
        // S >= p/2
        if p >= root { p - root } else { p + MOD - root }
    };

    let end = start.elapsed();
    println!("res = {minimal_root}, took {end:?}");
}
