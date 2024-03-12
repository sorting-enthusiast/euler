/*
let C(n) be the number of squarefree integers of the form x^2+1 such that 1 <= x <= n.
For example, C(10) = 9 and C(1000) = 895.
Find C(123567101113).
*/
//const N: usize = 123_567_101_113;
use crate::bit_array::BitArray;
use crate::eratosthenes_variants::{
    _wheel_factorized_sieve_of_eratosthenes, sieve_of_eratosthenes,
    simple_wheel_factorized_sieve_of_eratosthenes,
};
use crate::sieve_of_pritchard::sift;
const WHEEL_2_3_5_7: [u8; 48] = [
    10, 2, 4, 2, 4, 6, 2, 6, 4, 2, 4, 6, 6, 2, 6, 4, 2, 6, 4, 6, 8, 4, 2, 4, 2, 4, 8, 6, 4, 6, 2,
    4, 6, 2, 6, 6, 4, 2, 4, 6, 2, 6, 4, 2, 4, 2, 10, 2,
];
const WHEEL_2_3_5: [u8; 8] = [6, 4, 2, 4, 2, 4, 6, 2];

const WHEEL_2_3: [u8; 2] = [4, 2];

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
    let bitmap_size = (((n / 30) >> 3) + 1) << 3;
    let mut sieve = vec![0u8; bitmap_size].into_boxed_slice();
    let bitmap = sieve.as_mut_ptr();

    let mut num = 1;
    let mut wheel_incr = WHEEL_2_3_5_7.iter().cycle();
    primes.extend([2, 3, 5, 7]);

    loop {
        num += *wheel_incr.next().unwrap() as usize;
        if num * num > n {
            break;
        }
        unsafe {
            if *bitmap.add(num / 30) & (1 << ((num % 30) * 8 / 30)) == 0 {
                primes.push(num);
                let mut multiple = num * num;
                for &incr in wheel_incr.clone() {
                    if multiple > n {
                        break;
                    }
                    *bitmap.add(multiple / 30) |= 1 << ((multiple % 30) * 8 / 30);
                    multiple += incr as usize * num;
                }
            }
        }
    }
    for &incr in wheel_incr {
        if num > n {
            break;
        }
        unsafe {
            if *bitmap.add(num / 30) & (1 << ((num % 30) * 8 / 30)) == 0 {
                primes.push(num);
            }
        }
        num += incr as usize;
    }

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
    let mut sieve = BitArray::zeroed((limit + 1) / 2);

    let mut wheel_incr = WHEEL_2_3_5_7.iter().cycle();
    primes.extend([2, 3, 5, 7]);
    let mut num = 1;
    loop {
        let incr = *wheel_incr.next().unwrap();
        num += incr as usize;
        if num * num > limit {
            break;
        }
        if !sieve.get(num >> 1) {
            primes.push(num);
            let mut multiple = num * num;
            for &incr in wheel_incr.clone() {
                if multiple > limit {
                    break;
                }
                sieve.set(multiple >> 1);
                multiple += incr as usize * num;
            }
        }
    }
    for &incr in wheel_incr {
        if num > limit {
            break;
        }
        if !sieve.get(num >> 1) {
            primes.push(num);
        }
        num += incr as usize;
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
            let mut multiple = low - low % (prime * 6) + prime;
            for &incr in WHEEL_2_3.iter().cycle() {
                if multiple < low {
                    //println!("{prime}, {multiple}, {low}");
                    multiple += incr as usize * prime;
                    continue;
                }
                if multiple > high {
                    break;
                }
                sieve.set((multiple - low) >> 1);
                multiple += incr as usize * prime;
            }
        }

        let mut num = low - low % 210 + 1;

        for &incr in WHEEL_2_3_5_7.iter().cycle() {
            if num < low {
                num += incr as usize;
                continue;
            }
            if num > high {
                break;
            }
            if !sieve.get((num - low) >> 1) {
                primes.push(num);
            }
            num += incr as usize;
        }

        low += limit;
        high += limit;
    }
    primes
}
pub fn sieve_of_atkin(limit: usize) -> Vec<usize> {
    let end = limit + 1;
    let mut primes = Vec::with_capacity((limit as f64 / (limit as f64).log(3.0)) as usize);
    primes.extend([2, 3, 5]);
    let mut sieve = BitArray::zeroed(1 + (limit / 30) << 3);
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
    let sqrt_end = ((end as f64).sqrt() + 1.01) as usize;
    for x in 1..((sqrt_end as f64 / 2.0 + 2.01) as usize) {
        let xx4 = 4 * x * x;
        for y in (1..sqrt_end).step_by(2) {
            let n = xx4 + y * y;
            if n >= end {
                break;
            }
            if (SET0 >> (n % 60)) & 1 == 1 {
                let ind = ((n / 30) << 3) + ((n % 30) << 3) / 30;
                sieve.flip(ind);
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
                let ind = ((n / 30) << 3) + ((n % 30) << 3) / 30;
                sieve.flip(ind);
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
                let ind = ((n / 30) << 3) + ((n % 30) << 3) / 30;
                sieve.flip(ind);
            }
        }
    }
    let mut n = 1;
    for (&incr, ind) in WHEEL_2_3_5.iter().cycle().zip(0..) {
        if n > sqrt_end {
            break;
        }
        if sieve.get(ind) {
            let n_squared = n * n;
            let mut c = n_squared;
            for &incr in WHEEL_2_3_5.iter().cycle() {
                if c > end {
                    break;
                }
                let ind = ((c / 30) << 3) + ((c % 30) << 3) / 30;
                sieve.clear(ind);
                c += incr as usize * n_squared;
            }
        }
        n += incr as usize;
    }
    n = 1;
    for (&incr, ind) in WHEEL_2_3_5.iter().cycle().zip(1..) {
        n += incr as usize;
        if n > end {
            break;
        }
        if sieve.get(ind) {
            primes.push(n);
        }
    }
    primes
}

pub fn main() {
    use std::time::Instant;
    const N: usize = 1e9 as usize + 7;

    assert_eq!(
        sieve_of_atkin(500),
        sift(500)
            .into_iter()
            .map(|x| x as usize)
            .collect::<Vec<usize>>()
    );
    dbg!(sift(500));

    let start = Instant::now();
    dbg!(sift(N as u64).len());
    let end = start.elapsed();
    println!("{:?}", end);

    let start = Instant::now();
    dbg!(segmented_sieve_of_eratosthenes(N).len());
    let end = start.elapsed();
    println!("{:?}", end);

    let start = Instant::now();
    dbg!(simple_wheel_factorized_sieve_of_eratosthenes(N).len());
    let end = start.elapsed();
    println!("{:?}", end);

    let start = Instant::now();
    dbg!(wheel_factorized_sieve_of_eratosthenes(N).len());
    let end = start.elapsed();
    println!("{:?}", end);

    let start = Instant::now();
    dbg!(_wheel_factorized_sieve_of_eratosthenes(N).len());
    let end = start.elapsed();
    println!("{:?}", end);

    let start = Instant::now();
    dbg!(sieve_of_atkin(N).len());
    let end = start.elapsed();
    println!("{:?}", end);

    let start = Instant::now();
    dbg!(sieve_of_eratosthenes(N).len());
    let end = start.elapsed();
    println!("{:?}", end);

    let start = Instant::now();
    dbg!(sundaram_sieve(N).len());
    let end = start.elapsed();
    println!("{:?}", end);

    /* let mut num = 1;
    for (ind, &i) in WHEEL_2_3_5_7.iter().enumerate() {
        println!("{ind}: {} % 30 = {}, ", num, num % 30);
        num += i as usize;
    } */
}
