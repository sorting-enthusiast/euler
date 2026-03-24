use crate::utils::primes::wheel_sieve;

// the ogf of S is \prod_p\frac{1}{1-px^p}
pub fn main() {
    const N: usize = 46368;
    const MOD: usize = 1e9 as _;
    let start = std::time::Instant::now();
    let mut gf = [0; N + 1];
    gf[0] = 1;
    for p in wheel_sieve(N as _) {
        let p = p as usize;
        for i in p..=N {
            gf[i] += p * gf[i - p];
            gf[i] %= MOD;
        }
    }
    let mut res = 0;
    let [mut fk, mut fk1] = [1, 1];
    for _ in 1..24 {
        res += gf[fk1];
        if res >= MOD {
            res -= MOD;
        }
        [fk, fk1] = [fk1, fk + fk1];
    }
    println!("res = {res}, took {:?}", start.elapsed());
}
