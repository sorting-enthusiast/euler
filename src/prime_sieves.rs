/*
let C(n) be the number of squarefree integers of the form x^2+1 such that 1 <= x <= n.
For example, C(10) = 9 and C(1000) = 895.
Find C(123567101113).
*/
//const N: usize = 123_567_101_113;
use crate::bit_array::BitArray;
const WHEEL_INCR_2: [u8; 1] = [1];
const WHEEL_INCR_2_3: [u8; 2] = [1, 4];
const WHEEL_INCR_2_3_5: [u8; 8] = [1, 6, 4, 2, 4, 2, 4, 6];
const WHEEL_INCR_2_3_5_7: [u8; 48] = [
    1, 10, 2, 4, 2, 4, 6, 2, 6, 4, 2, 4, 6, 6, 2, 6, 4, 2, 6, 4, 6, 8, 4, 2, 4, 2, 4, 8, 6, 4, 6,
    2, 4, 6, 2, 6, 6, 4, 2, 4, 6, 2, 6, 4, 2, 4, 2, 10,
];
pub fn wheel_factorized_sieve_of_eratosthenes(n: usize) -> Vec<usize> {
    if n < 7 {
        let mut ret = vec![2];
        if n > 2 {
            ret.push(3);
        }
        if n > 4 {
            ret.push(5)
        }
        return ret;
    }
    let mut primes = Vec::with_capacity((n as f64 / (n as f64).log(3.0)) as usize);
    primes.extend([2, 3, 5, 7]);
    let mut sieve = BitArray::zeroed((n + 1) / 2);
    'outer: for w in 0.. {
        let mut num = 210 * w;
        for incr in WHEEL_INCR_2_3_5_7 {
            num += incr as usize;
            if num < 3 {
                continue;
            }
            if num * num > n {
                break 'outer;
            }
            if !sieve.get(num >> 1) {
                primes.push(num);
                'inner_outer: for w in num / 210.. {
                    let mut multiple = 210 * w * num;
                    for incr in WHEEL_INCR_2_3_5_7 {
                        multiple += incr as usize * num;
                        if multiple > n {
                            break 'inner_outer;
                        }
                        sieve.set(multiple >> 1);
                    }
                }
            }
        }
    }
    let sqrt_n = (n as f64).sqrt() as usize;
    'outer: for w in sqrt_n / 210.. {
        let mut num = 210 * w;
        for incr in WHEEL_INCR_2_3_5_7 {
            num += incr as usize;
            if num < sqrt_n {
                continue;
            }
            if num > n {
                break 'outer;
            }
            if !sieve.get(num >> 1) {
                primes.push(num);
            }
        }
    }
    primes
}
pub fn sieve_of_eratosthenes(n: usize) -> Vec<usize> {
    let mut primes = Vec::with_capacity((n as f64 / (n as f64).log(3.0)) as usize);
    primes.push(2);
    let mut sieve = BitArray::zeroed((n + 1) / 2);
    let sqrt_n = ((n as f64).sqrt().ceil() as usize) & !1;
    for num in (3..=sqrt_n).step_by(2) {
        if !sieve.get(num >> 1) {
            primes.push(num);
            for multiple in (num * num..=n).step_by(num << 1) {
                sieve.set(multiple >> 1);
            }
        }
    }
    primes.extend(
        ((sqrt_n + 1)..=n)
            .step_by(2)
            .filter(|&i| !sieve.get(i >> 1)),
    );
    primes
}
pub fn sundaram_sieve(n: usize) -> Vec<usize> {
    let k = (n - 1) / 2;
    let mut sieve = BitArray::zeroed(k);
    for i in 0..=((((n as f64).sqrt().ceil() as usize) - 3) / 2) {
        let p = 2 * i + 3;
        let s = (p * p - 3) / 2;
        for j in (s..k).step_by(p) {
            sieve.set(j);
        }
    }
    let mut primes = Vec::with_capacity((n as f64 / (n as f64).log(3.0)) as usize);
    primes.push(2);
    primes.extend((0..k).filter_map(|i| if sieve.get(i) { None } else { Some(2 * i + 3) }));
    primes
}
pub fn segmented_sieve_of_eratosthenes(n: usize) -> Vec<usize> {
    if n < 7 {
        let mut ret = vec![2];
        if n > 2 {
            ret.push(3);
        }
        if n > 4 {
            ret.push(5)
        }
        return ret;
    }
    let limit = (n as f64).sqrt().floor() as usize + 1;
    let mut primes = Vec::with_capacity((n as f64 / (n as f64).log(3.0)) as usize);
    primes.extend([2, 3, 5, 7]);

    let mut sieve = BitArray::zeroed((limit + 1) / 2);
    'outer: for w in 0.. {
        let mut num = 210 * w;
        for incr in WHEEL_INCR_2_3_5_7 {
            num += incr as usize;
            if num < 3 {
                continue;
            }
            if num * num > limit {
                break 'outer;
            }
            if !sieve.get(num >> 1) {
                primes.push(num);
                'inner_outer: for w in num / 210.. {
                    let mut multiple = 210 * w * num;
                    for incr in WHEEL_INCR_2_3_5_7 {
                        multiple += incr as usize * num;
                        if multiple > limit {
                            break 'inner_outer;
                        }
                        sieve.set(multiple >> 1);
                    }
                }
            }
        }
    }
    let sqrt_limit = (limit as f64).sqrt() as usize;
    'outer: for w in sqrt_limit / 210.. {
        let mut num = 210 * w;
        for incr in WHEEL_INCR_2_3_5_7 {
            num += incr as usize;
            if num < sqrt_limit {
                continue;
            }
            if num > limit {
                break 'outer;
            }
            if !sieve.get(num >> 1) {
                primes.push(num);
            }
        }
    }

    let first_primes_count = primes.len();

    let mut low = limit;
    let mut high = limit << 1;

    while low < n {
        sieve.zero();
        if high > n {
            high = n;
        }
        for &prime in &primes[1..first_primes_count] {
            let tmp = low / (6 * prime);
            'outer: for w in tmp.. {
                let mut multiple = 6 * w * prime;
                for incr in WHEEL_INCR_2_3 {
                    multiple += incr as usize * prime;
                    if multiple < low {
                        continue;
                    }
                    if multiple > high {
                        break 'outer;
                    }
                    sieve.set((multiple - low) >> 1);
                }
            }
        }
        'outer: for w in low / 210.. {
            let mut num = 210 * w;
            for incr in WHEEL_INCR_2_3_5_7 {
                num += incr as usize;
                if num < low {
                    continue;
                }
                if num > high {
                    break 'outer;
                }
                if !sieve.get((num - low) >> 1) {
                    primes.push(num);
                }
            }
        }
        low += limit;
        high += limit;
    }
    primes
}
pub fn sieve_of_atkin(limit: usize) -> Vec<usize> {
    let end = limit + 1;
    let mut sieve = BitArray::zeroed(end / 2);
    const SET0: u64 = 0
        | (1 << 1)
        | (1 << 13)
        | (1 << 17)
        | (1 << 29)
        | (1 << 37)
        | (1 << 41)
        | (1 << 49)
        | (1 << 53);
    const SET1: u64 = 0 | (1 << 7) | (1 << 19) | (1 << 31) | (1 << 43);
    const SET2: u64 = 0 | (1 << 11) | (1 << 23) | (1 << 47) | (1 << 59);
    const WHEEL_INCR: [u8; 16] = [1, 6, 4, 2, 4, 2, 4, 6, 2, 6, 4, 2, 4, 2, 4, 6];
    let sqrt_end = ((end as f64).sqrt() + 1.01) as usize;
    for x in 1..((sqrt_end as f64 / 2.0 + 2.01) as usize) {
        let xx4 = 4 * x * x;
        for y in (1..sqrt_end).step_by(2) {
            let n = xx4 + y * y;
            if n >= end {
                break;
            }
            if (SET0 >> (n % 60)) & 1 == 1 {
                sieve.flip(n >> 1);
            }
        }
    }
    for x in (1..((sqrt_end as f64 / 3.0f64.sqrt() + 2.01) as usize)).step_by(2) {
        let xx3 = 3 * x * x;
        for y in (2..sqrt_end).step_by(2) {
            let n = xx3 + y * y;
            if n >= end {
                break;
            }
            if (SET1 >> (n % 60)) & 1 == 1 {
                sieve.flip(n >> 1);
            }
        }
    }
    for x in 1..(sqrt_end as f64 / 2.0f64.sqrt() + 2.01) as usize {
        let xx3 = 3 * x * x;
        for y in (0..x).rev().step_by(2) {
            let n = xx3 - y * y;
            if n >= end {
                break;
            }
            if (SET2 >> (n % 60)) & 1 == 1 {
                sieve.flip(n >> 1);
            }
        }
    }
    'outer: for w in 0.. {
        let mut n = 60 * w;
        for incr in WHEEL_INCR {
            n += incr as usize;
            if n < 7 {
                continue;
            }
            if n > sqrt_end {
                break 'outer;
            }
            if sieve.get(n >> 1) {
                let n_squared = n * n;
                'inner_outer: for w in 0.. {
                    let mut c = 60 * w * n_squared;
                    for incr in WHEEL_INCR {
                        c += incr as usize * n_squared;
                        if c > end {
                            break 'inner_outer;
                        }
                        sieve.clear(c >> 1);
                    }
                }
            }
        }
    }
    let mut primes = Vec::with_capacity((limit as f64 / (limit as f64).log(3.0)) as usize);
    primes.extend([2, 3, 5]);
    'outer: for w in 0.. {
        let mut n = 60 * w;
        for incr in WHEEL_INCR {
            n += incr as usize;
            if n < 7 {
                continue;
            }
            if n > end {
                break 'outer;
            }
            if sieve.get(n >> 1) {
                primes.push(n);
            }
        }
    }
    primes
}

pub fn linear_segmented_wheel_sieve(n: usize) -> Vec<usize> {
    _ = n;
    todo!();
}

pub fn main() {
    use std::time::Instant;
    assert_eq!(
        sieve_of_atkin(169),
        wheel_factorized_sieve_of_eratosthenes(169)
    );
    dbg!(segmented_sieve_of_eratosthenes(500));
    const N: usize = 49; //1e9 as usize + 7;

    let start = Instant::now();
    dbg!(segmented_sieve_of_eratosthenes(N).len());
    let end = start.elapsed();
    println!("{:?}", end);

    let start = Instant::now();
    dbg!(wheel_factorized_sieve_of_eratosthenes(N).len());
    let end = start.elapsed();
    println!("{:?}", end);

    let start = Instant::now();
    dbg!(sieve_of_atkin(N).len());
    let end = start.elapsed();
    println!("{:?}", end);
}
