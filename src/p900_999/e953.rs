use crate::utils::prime_sieves::sift;

const MOD: u64 = 1e9 as u64 + 7;
const N: u64 = 1e14 as _;
const SQRT_N: u64 = N.isqrt();

// I was a massive idiot, original code would have run just fine had I stopped the recursion at the right time
pub fn main() {
    // the first player to move loses in factorization nim iff the xor sum of all
    // prime factors of n with multiplicity is 0. This is trivially the case for all squares,
    // and for numbers such as 70. For a given n, only primes that appear an odd number of
    // times in the factorization of n affect the xor sum.
    // all such n are of the form n = s^2 * m, with m squarefree, such that the xorsum of the primes of m is 0
    // so, each possible m contributes m * sum_squares(isqrt(N / m)) to the total.
    // moreover, all factors are less than sqrt(N), since otherwise the xorsum would not be able to be 0
    let primes = sift(SQRT_N);

    let start = std::time::Instant::now();
    let (res, count, total) = sum(0, 1, &primes[1..] /*, &mut solutions*/);
    println!(
        "res = {}, took {:?}",
        (res + sum_squares(SQRT_N)) % MOD,
        start.elapsed()
    );
    dbg!(count, total);

    let start = std::time::Instant::now();
    let mut primes_by_msb = [const { vec![] }; 23];
    for &p in &primes[..] {
        primes_by_msb[p.ilog2() as usize - 1].push(p);
    }
    let (res, count, total) = sum3(0, 1, 23, &primes_by_msb, SQRT_N, false);
    println!(
        "res = {}, took {:?}",
        (res + sum_squares(SQRT_N)) % MOD,
        start.elapsed()
    );
    dbg!(count, total);
}
const fn sum_squares(x: u64) -> u64 {
    let x = (x % (6 * MOD)) as u128;
    (((x * (x + 1) * (2 * x + 1)) / 6) % MOD as u128) as u64
}
fn sum(xor: u64, acc: u64, primes: &[u64]) -> (u64, u64, u64) {
    let mut res = 0;
    let mut count = 0;
    let mut total = 2;
    let lim = N / acc;
    if xor <= lim && primes.binary_search(&xor).is_ok() {
        count += 1;
        res += (((acc * xor) % MOD) * sum_squares((lim / xor).isqrt())) % MOD;
        if res >= MOD {
            res -= MOD;
        }
    }
    if (xor ^ 2) * 2 <= lim && primes.binary_search(&(xor ^ 2)).is_ok() {
        count += 1;
        res +=
            (((acc * ((xor ^ 2) * 2)) % MOD) * sum_squares((lim / (2 * (xor ^ 2))).isqrt())) % MOD;
        if res >= MOD {
            res -= MOD;
        }
    }
    for (i, &p) in primes.iter().enumerate() {
        if p * p > lim || (xor & 1 == 1 && p * p * p > lim) {
            break;
        }
        let (a, b, c) = sum(xor ^ p, acc * p, &primes[i + 1..]);
        res += a;
        count += b;
        total += c;
        if res >= MOD {
            res -= MOD;
        }
    }
    (res, count, total)
}

fn sum3(
    xor: u64,
    acc: u64,
    b: usize,
    primes_by_msb: &[Vec<u64>], /*, sol: &mut HashMap<u64,u8>*/
    lpf: u64,
    check: bool,
) -> (u64, u64, u64) {
    let mut res = 0;
    let mut count = 0;
    let mut total = 0;
    let lim = N / acc;
    if check {
        total += 1;
        if xor == 0 {
            count += 1;
            res += ((acc % MOD) * sum_squares(lim.isqrt())) % MOD;
            if res >= MOD {
                res -= MOD;
            }
        }
    }
    if b > 0 {
        if (xor >> b) & 1 == 0 {
            if !primes_by_msb[b - 1].is_empty()
                && primes_by_msb[b - 1][0] * primes_by_msb[b - 1][0] <= lim
            {
                for (i, &p1) in primes_by_msb[b - 1].iter().enumerate() {
                    if p1 >= lpf || p1 > lim {
                        break;
                    }
                    for &p2 in &primes_by_msb[b - 1][i + 1..] {
                        if p2 >= lpf || (p1 * p2) > lim {
                            break;
                        }

                        let (a, b, c) =
                            sum3(xor ^ p1 ^ p2, acc * p1 * p2, b, primes_by_msb, p1, true);
                        res += a;
                        count += b;
                        total += c;
                        if res >= MOD {
                            res -= MOD;
                        }
                    }
                }
            }
            let (a, b, c) = sum3(xor, acc, b - 1, primes_by_msb, lpf, false);
            res += a;
            count += b;
            total += c;
            if res >= MOD {
                res -= MOD;
            }
        } else {
            for &p in &primes_by_msb[b - 1] {
                if p >= lpf || p > lim {
                    break;
                }
                let (a, b, c) = sum3(xor ^ p, acc * p, b, primes_by_msb, p, true);
                res += a;
                count += b;
                total += c;
                if res >= MOD {
                    res -= MOD;
                }
            }
        }
    }
    (res, count, total)
}
