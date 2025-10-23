use crate::utils::{FIArray::FIArrayU64, prime_sieves::sift};
const MOD: u64 = 1e9 as u64 + 7;

pub fn main() {
    // the first player to move loses in factorization nim iff the xor sum of all
    // prime factors of n with multiplicity is 0. This is trivially the case for all squares,
    // and for numbers such as 70. For a given n, only primes that appear an odd number of
    // times in the factorization of n affect the xor sum.
    // all such n are of the form n = s^2 * m, with m squarefree, such that the xorsum of the primes of m is 0
    // so, for each possible s, sum all squarefree m such that m <= N/(s^2) and m is special.
    // moreover, all factors are less than sqrt(N), since otherwise the xorsum would not be able to be 0
    const N: u64 = 1e11 as u64;
    const SQRT_N: u64 = N.isqrt();
    let start = std::time::Instant::now();
    let primes = sift(SQRT_N);
    println!("primes generated: {:?}", start.elapsed());
    dbg!(primes.len());
    let start = std::time::Instant::now();

    let mut sums = FIArrayU64::new(N);
    sums.arr.fill(!0);
    let mut xx_mod = 1;
    let mut res = 0;
    for x in 2..=SQRT_N {
        xx_mod += (2 * x - 1) % MOD;
        if xx_mod >= MOD {
            xx_mod -= MOD;
        }
        let i = sums.get_index(N / (x * x));
        if sums.arr[i] == !0 {
            sums.arr[i] = 1 + sum(N / (x * x), 0, 1, &primes);
            println!("{x}, {:?}", start.elapsed());

            if sums.arr[i] >= MOD {
                sums.arr[i] -= MOD;
            }
        }
        res += (xx_mod * sums.arr[i]) % MOD;
        if res >= MOD {
            res -= MOD;
        }
    }
    res += 1 + sum(N, 0, 1, &primes);
    res %= MOD;
    println!("res = {res}, took {:?}", start.elapsed());

    dbg!(
        sums.arr.len(),
        sums.arr.iter().filter(|&&e| e != !0).count()
    );
    /* dbg!(
        (1..=SQRT_N)
            .map(|x| x * x * (1 + sum(N / (x * x), 0, 1, &primes)))
            .sum::<u64>()
            % MOD
    ); */
}
fn sum(hi: u64, xor: u64, acc: u64, primes: &[u64]) -> u64 {
    let mut res = 0;
    for (i, &p) in primes.iter().enumerate() {
        let prod = acc * p;
        if prod > hi {
            break;
        }
        let x = xor ^ p;
        res += sum(hi, x, prod, &primes[i + 1..]);
        if res >= MOD {
            res -= MOD;
        }
        if x == 0 {
            res += prod;
            if res >= MOD {
                res -= MOD;
            }
        }
    }
    res
}
