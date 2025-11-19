use itertools::{Either, Itertools};
use std::collections::HashSet;

use crate::utils::{
    FIArray::FIArrayU64,
    bit_array::BitArray,
    primes::{primality::powmod, prime_sieves::sift},
};

const fn modexp(mut x: usize, mut exp: usize, modulo: usize) -> usize {
    if exp == 0 {
        return 1;
    }
    let mut r = 1;
    while exp > 1 {
        if exp & 1 == 1 {
            r = (r * x) % modulo;
        }
        x = (x * x) % modulo;
        exp >>= 1;
    }
    (r * x) % modulo
}
fn find_divisor(p: usize) -> (usize, usize) {
    assert_eq!(
        p & 3,
        1,
        "cannot find divisor, prime is not equivalent to 1 mod 4"
    );
    let r = (2..p)
        .map(|i| modexp(i, (p - 1) >> 2, p))
        .find(|r| p - 1 == (r * r) % p)
        .expect("divisor was not found, something is very wrong");
    let mut a = p;
    let mut b = r;

    // Euclidean algorithm
    while b * b > p {
        let tmp = a % b;
        a = b;
        b = tmp;
    }

    (b, (p - b * b).isqrt())
}
fn sieve(n: usize) -> Vec<usize> {
    let mut res = vec![0; n];
    let mut pow = vec![0; n];
    //let mut composite = BitArray::zeroed(n);
    let mut primes = vec![];
    let mut prime_power_divisors = vec![];
    res[1] = 1;
    for i in 2..n {
        if res[i] == 0
        /* !composite.get(i) */
        {
            //dbg!(!composite.get(i) && (res[i] == 0));
            primes.push(i);
            pow[i] = i;
            res[i] = i
                + 1
                + match i & 3 {
                    1 => {
                        let (x, y) = find_divisor(i);
                        assert_eq!(x * x + y * y, i);
                        if i * i < n {
                            prime_power_divisors.push(Some((
                                (x as i64, y as i64),
                                vec![
                                    (x as i64, y as i64),
                                    (x as i64, -(y as i64)),
                                    (y as i64, x as i64),
                                    (y as i64, -(x as i64)),
                                ],
                            )));
                        }
                        (x + y) << usize::from(x != y)
                    }
                    2 => {
                        if i * i < n {
                            prime_power_divisors.push(None); /* 
                            prime_power_divisors
                            .push(Some(((1, 1), vec![(1, 0), (1, 1), (1, -1), (2, 0)]))); */
                        }
                        2
                    }
                    3 => {
                        if i * i < n {
                            prime_power_divisors.push(None);
                        }
                        0
                    }
                    _ => unsafe { core::hint::unreachable_unchecked() },
                }
        }
        for (j, &p) in primes.iter().enumerate() {
            if i * p >= n {
                break;
            }
            //composite.set(i * p);
            if i.is_multiple_of(p) {
                pow[i * p] = pow[i] * p;
                let v = i / pow[i];
                if v != 1 {
                    assert_eq!(v * pow[i], i);
                    assert_ne!(v % p, 0);
                    if p & 3 == 3 {
                        res[i * p] = res[v] * res[p * pow[i]];
                    } else {
                        todo!()
                    }
                } else
                //prime power
                if p & 3 != 3 {
                    //has complex divisors :(
                    if p == 2 {
                        res[i << 1] = res[i] + (i << 2);
                    } else {
                        let Some(((x, y), ref mut v)) = prime_power_divisors[j] else {
                            unsafe { core::hint::unreachable_unchecked() }
                        };
                        println!("{p}: {x} {y}, {}", i * p);
                        for t in v.iter() {
                            print!("{t:?} ");
                        }
                        println!();
                        let len = v.len();
                        for (a, b) in [
                            /* (1, 0),  */ (x, y),
                            (x, -y),
                            (y, x),
                            (y, -x),
                            (p as i64, 0),
                        ] {
                            for k in 0..len {
                                let (c, d) = v[k];
                                let acbd = a * c + b * d;
                                let bcad = b * c - a * d;
                                if acbd <= 0 || v.contains(&(acbd, bcad)) {
                                    continue;
                                }
                                v.push((acbd, bcad));
                            }
                        }
                        v[len..].sort_by_key(|(a, b)| a.pow(2) + b.pow(2));

                        for t in &v[len..] {
                            print!("{t:?} ");
                        }
                        println!("\n");
                        assert_eq!(v.iter().map(|t| t.1).sum::<i64>(), 0);
                        assert_eq!(res[i], v[..len].iter().map(|t| t.0 as usize).sum::<usize>());
                        res[i * p] = res[i] + v[len..].iter().map(|t| t.0 as usize).sum::<usize>();
                    }
                } else {
                    res[i * p] = res[i] + i * p;
                }
                break;
            }
            assert_eq!(res[i * p], 0);
            res[i * p] = res[i] * res[p];
            pow[i * p] = p;
        }
    }
    res
}

const fn sum_n(x: u64) -> u64 {
    if x & 1 == 0 {
        (x / 2) * (x + 1)
    } else {
        x.div_ceil(2) * x
    }
}
const fn sum_divisor_summatory(n: u64) -> u64 {
    let mut sum = 0;
    let sqrtn = n.isqrt();
    let mut i = 1;
    while i <= sqrtn {
        let ni = n / i;
        sum += i * ni + sum_n(ni);
        i += 1;
    }
    sum - sqrtn * sum_n(sqrtn)
}

const N: u64 = 7 as u64;
// sum all numbers below lim s.t. all their prime factors are from the prime set given
fn dfs(acc: u64, lim: u64, primes: &[u64]) -> u64 {
    let mut sum = acc;
    for (i, p2) in primes.iter().map(|p| p.pow(2)).enumerate() {
        if p2 * acc > lim {
            break;
        }
        let mut power = p2;
        loop {
            sum += dfs(acc * power, lim, &primes[i + 1..]);
            power *= p2;
            if power * acc > lim {
                break;
            }
        }
    }
    sum
}
fn wrapper(lim: u64, primes: &[u64]) -> u64 {
    let mut acc = 1;
    let mut sum = 0;
    while acc <= lim {
        sum += dfs(acc, lim, primes);
        acc <<= 1;
    }
    sum
}

fn sum_of_squares(p: u64) -> (u64, u64) {
    assert_eq!(
        p & 3,
        1,
        "cannot find divisor, prime is not equivalent to 1 mod 4"
    );
    let r = (2..p)
        .map(|i| powmod(i, (p - 1) >> 2, p))
        .find(|r| p - 1 == (r * r) % p)
        .expect("divisor was not found, something is very wrong");
    let mut a = p;
    let mut b = r;

    // Euclidean algorithm
    while b * b > p {
        let tmp = a % b;
        a = b;
        b = tmp;
    }

    (b, (p - b * b).isqrt())
}
fn count(
    acc: u64,
    lim: u64,
    primes: &[u64],
    points: &[(i64, i64)],
    factors: &mut Vec<((i64, i64), u8)>,
    dfs_cache: &FIArrayU64,
) -> u64 {
    //dbg!(acc);
    let mut res = /* dbg! */(sum_divisors(factors)) * /* dbg! */(dfs_cache[N / acc]);
    for (i, (&p, &point)) in primes.iter().zip(points.iter()).enumerate() {
        if p * acc > lim {
            break;
        }
        let mut power = p;
        factors.push((point, 0));
        loop {
            factors.last_mut().unwrap().1 += 1;
            res += count(
                acc * power,
                lim,
                &primes[i + 1..],
                &points[i + 1..],
                factors,
                dfs_cache,
            );
            power *= p;
            if power * acc > lim {
                break;
            }
        }
        factors.pop();
    }
    res
}
const fn add_gaussian((a, b): (i64, i64), (c, d): (i64, i64)) -> (i64, i64) {
    (a + c, b + d)
}
const fn mul_gaussian((a, b): (i64, i64), (c, d): (i64, i64)) -> (i64, i64) {
    (a * c - b * d, a * d + b * c)
}
fn sum_divisors(factors: &[((i64, i64), u8)]) -> u64 {
    fn dfs(remaining: &[((i64, i64), u8)], acc: (i64, i64), seen: &mut HashSet<(i64, i64)>) {
        if remaining.is_empty() {
            if acc.0 > 0 && acc.1 > 0 {
                seen.insert(acc);
            }
            return;
        }

        let &[((a, b), e), ref rest @ ..] = remaining else {
            unsafe { core::hint::unreachable_unchecked() }
        };
        let pi = (a, b);
        let pi_conj = (a, -b);

        let mut pow_pi = (1, 0); // pi^a
        for a in 0..=e {
            let mut pow_conj = (1, 0); // (conj pi)^b
            for b in 0..=e {
                let next = mul_gaussian(mul_gaussian(acc, pow_pi), pow_conj);
                //dbg!(next);
                dfs(rest, next, seen);
                pow_conj = mul_gaussian(pow_conj, pi_conj);
            }
            pow_pi = mul_gaussian(pow_pi, pi);
        }
    }
    let mut seen = HashSet::new();
    dfs(factors, (1, 0), &mut seen);
    //dbg!(&seen);
    seen.into_iter()
        .fold(0, |acc, (a, b)| acc + ((a + b) << i64::from(a != b))) as u64
}
// can probably do smth like 233, and precompute partial sums to store in an FIArray, since I only need the values at N/i for some i
fn solve() {
    let start = std::time::Instant::now();

    let mut sum = const { sum_divisor_summatory(N) };
    dbg!(sum);
    let primes = sift(N);
    let (primes_1mod4, rest): (Vec<_>, Vec<_>) = primes.into_iter().partition_map(|p| {
        if p & 3 == 1 {
            Either::Left(p)
        } else {
            Either::Right(p)
        }
    });
    dbg!(primes_1mod4.len());
    let points = primes_1mod4
        .iter()
        .copied()
        .map(sum_of_squares)
        .map(|(a, b)| (a as i64, b as i64))
        .collect_vec();
    dbg!(&rest[1..]);
    let mut dfs_cache = FIArrayU64::new(N);
    for (i, v) in FIArrayU64::keys(N).enumerate() {
        dfs_cache.arr[i] = wrapper(v, &rest[1..]);
    }
    //dbg!(&dfs_cache);
    sum += dbg!(count(1, N, &primes_1mod4, &points, &mut vec![], &dfs_cache));
    let k = N.ilog2() as u64;

    sum += 2 * ((1 << (k + 1)) - 2 - k);
    let end = start.elapsed();
    println!("res = {sum}, took {end:?}");
    println!("{sum}");
    println!("17924657155");
    dbg!(wrapper(30, &[3, 5]));
    dbg!(sum_divisors(&[((2, 1), 1), ((3, 2), 1)]));
    // 2*(2^(k+1)-2-k)
}
pub fn main() {
    solve();
    /* let expected = [
        1, 6, 10, 23, 35, 55, 63, 92, 105, 161, 173, 225, 249, 289, 337, 398, 426, 491, 511, 655,
        687, 747, 771, 887, 968, 1080, 1120, 1224, 1268, 1492, 1524, 1649, 1697, 1833, 1929, 2098,
        2150, 2250, 2346, 2666, 2726, 2886, 2930, 3086, 3242, 3362, 3410, 3654, 3711,
    ];
    for t in measured {
        print!("{t} ");
    }
    println!();
    for t in expected {
        print!("{t} ");
    }
    println!(); */
}
