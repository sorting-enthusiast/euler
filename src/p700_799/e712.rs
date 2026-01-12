use crate::utils::primes::primecount::lucy_fenwick;

const N: usize = 1e12 as _;
const SQRT_N: usize = N.isqrt();
const MOD: usize = 1e9 as usize + 7;

pub fn main() {
    const N_MOD: usize = N % MOD;
    let start = std::time::Instant::now();
    let pi = lucy_fenwick(N);
    let primes_done = start.elapsed();
    let mut res = 0;
    for p in 2..=SQRT_N {
        if pi.arr[p - 1] != pi.arr[p - 2] {
            let mut lim = N;
            while lim >= p {
                lim /= p;
                res += ((lim % MOD) * (N_MOD + MOD - lim % MOD)) % MOD;
                if res >= MOD {
                    res -= MOD;
                }
            }
        }
    }
    for d in (1..SQRT_N).chain((SQRT_N != N / SQRT_N).then_some(SQRT_N)) {
        res += (((d % MOD) * (N_MOD + MOD - d % MOD)) % MOD) * (pi[N / d] - pi[N / (d + 1)]) % MOD;
        if res >= MOD {
            res -= MOD;
        }
    }
    res <<= 1;
    if res >= MOD {
        res -= MOD;
    }
    println!("res = {res}, took {:?}", start.elapsed());
    dbg!(primes_done);
}
