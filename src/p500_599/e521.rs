use itertools::Itertools;

use crate::utils::FIArray::FIArray;

//euler 521
//could do incl excl principle
//will use modified segmented sieve of eratosthenes
const WHEEL_2_3_5: [u8; 8] = [6, 4, 2, 4, 2, 4, 6, 2];
const MOD30_TO_BIT8: [u64; 30] = [
    !0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7,
];
const MOD30_TO_MASK: [u8; 30] = [
    0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 4, 4, 8, 8, 8, 8, 16, 16, 32, 32, 32, 32, 64, 64, 64, 64, 64,
    64, 128,
];

#[must_use]
pub fn segmented_sieve_of_eratosthenes(n: usize) -> u128 {
    let limit = (n as f64).sqrt().floor() as usize + 1;
    let mut primes = Vec::with_capacity((limit as f64 / (limit as f64).log(3.0)) as usize);
    let bitmap_size = (((limit / 30) >> 3) + 1) << 3;
    let mut sieve = vec![0u8; bitmap_size].into_boxed_slice();
    let bitmap = sieve.as_mut_ptr();

    let mut num = 1;
    let mut wheel_incr = WHEEL_2_3_5.into_iter().cycle();
    let mut cnt =
        (n / 2 * 2 + n / 3 * 3 + n / 5 * 5 - n / 6 * 3 - n / 10 * 5 - n / 15 * 5 + n / 30 * 5) as _;

    loop {
        num += wheel_incr.next().unwrap() as usize;
        if num * num > limit {
            break;
        }
        unsafe {
            if *bitmap.add(num / 30) & MOD30_TO_MASK[num % 30] == 0 {
                cnt += num as u128;
                primes.push(num);
                let mut multiple = num * num;
                let precomp = [0, 0, num << 1, 0, num << 2, 0, num * 6]; //to reduce bit complexity, multiply only twice instead of every iteration of inner loop

                for incr in wheel_incr.clone() {
                    if multiple > limit {
                        break;
                    }
                    if *bitmap.add(multiple / 30) & MOD30_TO_MASK[multiple % 30] == 0 {
                        cnt += num as u128;
                    }

                    *bitmap.add(multiple / 30) |= MOD30_TO_MASK[multiple % 30];
                    multiple += precomp[incr as usize];
                }
            }
        }
    }
    for incr in wheel_incr {
        if num > limit {
            break;
        }
        unsafe {
            if *bitmap.add(num / 30) & MOD30_TO_MASK[num % 30] == 0 {
                primes.push(num);
                cnt += num as u128;
            }
        }
        num += incr as usize;
    }

    let mut low = limit;
    let mut high = limit << 1;

    let mut first_multiple = Vec::with_capacity(primes.len());
    for &prime in &primes {
        let mut multiple = prime * prime;
        let precomp = [0, 0, prime << 1, 0, prime << 2, 0, prime * 6]; //to reduce bit complexity, multiply only twice instead of every iteration of inner loop

        if multiple <= low {
            for incr in WHEEL_2_3_5
                .into_iter()
                .cycle()
                .skip(MOD30_TO_BIT8[prime % 30] as usize)
            {
                let tmp = multiple + precomp[incr as usize];
                if tmp > high {
                    break;
                }
                multiple += precomp[incr as usize];
                if tmp > low {
                    break;
                }
            }
        }
        first_multiple.push(multiple);
    }
    //let mut i = 1;
    while low < n {
        /* if i % 10000 == 0 {
            println!("{i}");
        }
        i += 1; */
        unsafe {
            core::ptr::write_bytes(bitmap, 0, bitmap_size);
        }
        if high > n {
            high = n;
        }
        for (&prime, multiple) in primes.iter().zip(first_multiple.iter_mut()) {
            let mut to_erase = *multiple;
            if to_erase > high {
                continue;
            }
            let precomp = [0, 0, prime << 1, 0, prime << 2, 0, prime * 6];
            for incr in WHEEL_2_3_5
                .into_iter()
                .cycle()
                .skip(MOD30_TO_BIT8[(to_erase / prime) % 30] as usize)
            {
                unsafe {
                    if *bitmap.add((to_erase - low) / 30) & MOD30_TO_MASK[to_erase % 30] == 0 {
                        cnt += prime as u128;
                    }
                    *bitmap.add((to_erase - low) / 30) |= MOD30_TO_MASK[to_erase % 30];
                }
                to_erase += precomp[incr as usize];
                if to_erase > high {
                    break;
                }
            }
            *multiple = to_erase;
        }
        let mut num = low - low % 30 + 1;

        for &incr in WHEEL_2_3_5.iter().cycle() {
            if num <= low {
                num += incr as usize;
                continue;
            }
            if num > high {
                break;
            }
            unsafe {
                if *bitmap.add((num - low) / 30) & MOD30_TO_MASK[num % 30] == 0 {
                    cnt += num as u128;
                }
            }

            num += incr as usize;
        }
        low += limit;
        high += limit;
    }
    cnt
}

const MOD: usize = 1e9 as usize;
pub fn main() {
    let n = 1e12 as _;
    //let first_primes = sift(n as _);
    //dbg!(incl_excl(n as _, 1, 1, &first_primes));
    let start = std::time::Instant::now();
    //dbg!(segmented_sieve_of_eratosthenes(n) % MOD as u128);
    let sum = lucy_based(n);
    let end = start.elapsed();
    println!("{sum}, {end:?}");
}
const fn sum_n<const MOD: usize>(x: usize) -> usize {
    let x = x % (MOD << 1);
    (if x & 1 == 0 {
        (x / 2) * (x + 1)
    } else {
        x.div_ceil(2) * x
    }) % MOD
}

fn lucy_based(x: usize) -> usize {
    let mut sum = 0;
    let mut s = FIArray::new(x);
    let mut sums = FIArray::new(x);
    let keys = FIArray::keys(x).collect_vec();
    for (i, &v) in keys.iter().enumerate() {
        s.arr[i] = v - 1;
        sums.arr[i] = (sum_n::<MOD>(v) - 1 + MOD) % MOD;
    }
    for p in 2..=x.isqrt() {
        if s.arr[p - 1] == s.arr[p - 2] {
            continue;
        }
        let sp = s.arr[p - 2];
        let sp2 = sums.arr[p - 2];
        for (i, &v) in keys.iter().enumerate().rev() {
            if v < p * p {
                break;
            }
            s.arr[i] -= s[v / p] - sp;
            sums.arr[i] += MOD - (p * (sums[v / p] + MOD - sp2)) % MOD;
            sums.arr[i] %= MOD;
            if v == x {
                sum += (p * (s[v / p] - sp)) % MOD; // p = lpf
                if sum >= MOD {
                    sum -= MOD;
                }
            }
        }
    }
    (sum + sums[x]) % MOD
}
/* use crate::sieve_of_pritchard::sift;

fn incl_excl(limit: u128, acc: u128, largest_pmf: u128, primes: &[u64]) -> u128 {
    let mut res = 0;
    for i in 0..primes.len() {
        let p = primes[i] as u128;
        let new_largest_pmf = if p > largest_pmf { p } else { largest_pmf };
        let prod = acc * p;
        if prod > limit {
            break;
        }
        res += limit / prod * new_largest_pmf;
        res -= incl_excl(limit, prod, new_largest_pmf, &primes[i + 1..]);
    }
    res
}
 */
