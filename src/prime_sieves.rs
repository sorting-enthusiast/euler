use std::ptr::write_bytes;
use std::vec;

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
const WHEEL_2_3_5_VALS: [u8; 8] = [1, 7, 11, 13, 17, 19, 23, 29];
#[inline(always)]
pub const fn num2bits(num: usize) -> u8 {
    (num % 30 * 8 / 30) as u8
}
//TODO: fix
pub fn hprime_sieve(max: usize) -> Vec<usize> {
    let mut adj_max = (max - 1) / 30 * 30 + 1;
    let sqrt_max = (max as f64).sqrt() as usize;
    let mut sieve = vec![0u8; max / 30 + 2 * 32 * 1024 + sqrt_max + 1].into_boxed_slice();
    let mut AS = vec![0u32; sqrt_max / 4 + 1000].into_boxed_slice();
    let mut BMPS = vec![0u32; sqrt_max / 4 + 1000].into_boxed_slice();
    let mut OFFS = vec![0u32; 8 * (sqrt_max / 4 + 1000)].into_boxed_slice();

    let offs = OFFS.as_mut_ptr();
    let bmps = BMPS.as_mut_ptr();
    let as_ptr = AS.as_mut_ptr();

    let mut cnt = 0;
    let mut a_i = 0;
    while adj_max + WHEEL_2_3_5[a_i & 7] as usize <= max {
        adj_max += WHEEL_2_3_5[a_i & 7] as usize;
        a_i += 1;
    }
    a_i = 0;
    let mut pattern = [[0u8; 16 * 8]; 3];
    let mut primes = vec![2, 3, 5];
    let bitmap = sieve.as_mut_ptr();

    for i in 1..4 {
        let mut r = WHEEL_2_3_5_VALS[i];
        while r < WHEEL_2_3_5_VALS[i] * 30 * 8 {
            pattern[i - 1][r as usize / 30] |= 1 << num2bits(r as usize);
            r += WHEEL_2_3_5_VALS[WHEEL_2_3_5[a_i & 7] as usize];
            a_i += 1;
        }
    }
    unsafe {
        let mut bitmap_end = bitmap.add(sqrt_max + 1);
        for i in 1..4 {
            let mut bmp = bitmap;
            while bmp < bitmap_end {
                for j in 0..WHEEL_2_3_5_VALS[i] as usize {
                    *(bmp as *mut u64).add(j) |= *((pattern[i - 1].as_ptr() as *const u64).add(j));
                }
                bmp = bmp.add(WHEEL_2_3_5_VALS[i] as usize * 8);
            }
        }
        a_i = 4;
        let mut a = 17u32;
        while a as usize <= sqrt_max {
            if *bitmap.add(a as usize / 30) & 1 << (a_i & 7) == 0 {
                *as_ptr.add(cnt) = a;
                for i in 0..8 {
                    *offs.add(
                        8 * cnt + num2bits((a * WHEEL_2_3_5_VALS[i] as u32) as usize) as usize,
                    ) = a * WHEEL_2_3_5_VALS[i] as u32 / 30;
                }

                let mut bmp = bitmap.add(a as usize * (a / 30) as usize);
                while bmp < bitmap_end {
                    for i in 0..8 {
                        *bmp.add(*offs.add(cnt * 8 + i) as usize) |= 1 << i;
                    }
                    bmp = bmp.add(a as usize);
                }
                *bmps.add(cnt) = bmp.offset_from(bitmap) as u32;
                cnt += 1;
            }

            a += WHEEL_2_3_5_VALS[a_i & 7] as u32;
            a_i += 1;
        }

        bitmap_end = bitmap.add(32 * 1024);
        for cur in (0..=max).step_by(32 * 1024 * 30) {
            for i in 1..4 {
                let mut bmp = bitmap.add(
                    cur / 30 / (8 * WHEEL_2_3_5_VALS[i] as usize)
                        * (8 * WHEEL_2_3_5_VALS[i] as usize),
                );
                while bmp < bitmap_end {
                    for j in 0..WHEEL_2_3_5_VALS[i] as usize {
                        *(bmp as *mut u64).add(j) |=
                            *((pattern[i - 1].as_ptr() as *const u64).add(j));
                    }
                    bmp = bmp.add(WHEEL_2_3_5_VALS[i] as usize * 8);
                }
            }

            for i in 0..cnt {
                while (*bmps.add(i) as usize) < cur / 30 + 32 * 1024 {
                    for j in 0..8 {
                        *bitmap.add((*bmps.add(i) + *offs.add(8 * i + j)) as usize) |= 1 << j;
                    }
                }
            }

            bitmap_end = bitmap_end.add(32 * 1024);
        }

        *bitmap = 0x1;

        let mut num = 1;
        for &incr in WHEEL_2_3_5.iter().cycle() {
            num += incr as usize;
            if num > max {
                break;
            }
            if *bitmap.add(num / 30) & 1 << num2bits(num) == 0 {
                primes.push(num);
            }
        }
    }

    primes
}

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
    let bitmap_size = (((limit / 30) >> 3) + 1) << 3;
    let mut sieve = vec![0u8; bitmap_size].into_boxed_slice();
    let bitmap = sieve.as_mut_ptr();

    let mut num = 1;
    let mut wheel_incr = WHEEL_2_3_5_7.iter().cycle();
    primes.extend([2, 3, 5, 7]);

    loop {
        num += *wheel_incr.next().unwrap() as usize;
        if num * num > limit {
            break;
        }
        unsafe {
            if *bitmap.add(num / 30) & (1 << ((num % 30) * 8 / 30)) == 0 {
                primes.push(num);
                let mut multiple = num * num;
                for &incr in wheel_incr.clone() {
                    if multiple > limit {
                        break;
                    }
                    *bitmap.add(multiple / 30) |= 1 << ((multiple % 30) * 8 / 30);
                    multiple += incr as usize * num;
                }
            }
        }
    }
    for &incr in wheel_incr {
        if num > limit {
            break;
        }
        unsafe {
            if *bitmap.add(num / 30) & (1 << ((num % 30) * 8 / 30)) == 0 {
                primes.push(num);
            }
        }
        num += incr as usize;
    }
    let first_primes_count = primes.len();
    let mut low = limit;
    let mut high = limit << 1;

    let mut first_multiple = Vec::with_capacity(first_primes_count);
    for &prime in &primes[3..] {
        let mut multiple = prime;
        if multiple * prime <= low {
            for &incr in WHEEL_2_3_5.iter().cycle().skip(((prime % 30) * 8) / 30) {
                let tmp = (multiple + incr as usize) * prime;
                if tmp > high {
                    break;
                }
                multiple += incr as usize;
                if tmp > low {
                    break;
                }
            }
        }
        first_multiple.push(multiple);
    }

    while low < n {
        unsafe {
            write_bytes(bitmap, 0, bitmap_size);
        }
        if high > n {
            high = n;
        }
        for (&prime, multiple) in primes[3..first_primes_count]
            .iter()
            .zip(first_multiple.iter_mut())
        {
            let mut to_erase = (*multiple) * prime;
            if to_erase > high {
                continue;
            }
            for &incr in WHEEL_2_3_5
                .iter()
                .cycle()
                .skip((((*multiple) % 30) * 8) / 30)
            {
                unsafe {
                    *bitmap.add((to_erase - low) / 30) |= 1 << (((to_erase % 30) * 8) / 30);
                }
                *multiple += incr as usize;
                to_erase += incr as usize * prime;
                if to_erase > high {
                    break;
                }
            }
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
                if *bitmap.add((num - low) / 30) & (1 << (((num % 30) * 8) / 30)) == 0 {
                    primes.push(num);
                }
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

    let bitmap_size = (((limit / 30) >> 3) + 1) << 3;
    let mut sieve = vec![0u8; bitmap_size].into_boxed_slice();
    let bitmap = sieve.as_mut_ptr();

    const SET0: u64 = (1 << 1)
        | (1 << 13)
        | (1 << 17)
        | (1 << 29)
        | (1 << 37)
        | (1 << 41)
        | (1 << 49)
        | (1 << 53);
    const SET1: u64 = (1 << 7) | (1 << 19) | (1 << 31) | (1 << 43);
    const SET2: u64 = (1 << 11) | (1 << 23) | (1 << 47) | (1 << 59);
    let sqrt_end = ((end as f64).sqrt() + 1.01) as usize;
    for x in 1..((sqrt_end as f64 / 2.0 + 2.01) as usize) {
        let xx4 = 4 * x * x;
        for y in (1..sqrt_end).step_by(2) {
            let n = xx4 + y * y;
            if n >= end {
                break;
            }
            if (SET0 >> (n % 60)) & 1 == 1 {
                unsafe {
                    *bitmap.add(n / 30) ^= 1 << ((n % 30) * 8 / 30);
                }
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
                unsafe {
                    *bitmap.add(n / 30) ^= 1 << ((n % 30) * 8 / 30);
                }
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
                unsafe {
                    *bitmap.add(n / 30) ^= 1 << ((n % 30) * 8 / 30);
                }
            }
        }
    }
    let mut n = 1;
    for &incr in WHEEL_2_3_5.iter().cycle() {
        if n > sqrt_end {
            break;
        }
        unsafe {
            if *bitmap.add(n / 30) & (1 << (((n % 30) * 8) / 30)) != 0 {
                let n_squared = n * n;
                let mut c = n_squared;
                for &incr in WHEEL_2_3_5.iter().cycle() {
                    if c > end {
                        break;
                    }

                    *bitmap.add(c / 30) &= !(1 << ((c % 30) * 8 / 30));

                    c += incr as usize * n_squared;
                }
            }
        }
        n += incr as usize;
    }
    n = 1;
    for &incr in WHEEL_2_3_5.iter().cycle() {
        n += incr as usize;
        if n > end {
            break;
        }
        unsafe {
            if *bitmap.add(n / 30) & (1 << (((n % 30) * 8) / 30)) != 0 {
                primes.push(n);
            }
        }
    }
    primes
}

pub fn rolling_sieve(n: usize) -> Vec<usize> {
    sieve_it().take_while(|&p| p <= n).collect()
}

/// Sorenson's rolling sieve - return iterator with no end.
/// ```
/// let l: Vec<usize> = sieve_it().take_while(|&p| p < 10).collect();
/// assert_eq!(l, vec![2, 3, 5, 7]);
/// ```
pub fn sieve_it() -> impl Iterator<Item = usize> {
    let start = 100;
    let mut r = f64::sqrt(start as f64) as usize + 1;
    let mut s = r * r;
    let mut delta = r + 2;
    let mut t = vec![vec![]; delta + 1];
    let primes = [
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89,
        97,
    ];

    for p in primes.into_iter().take_while(|&p| p < r) {
        let j = (p - (start % p)) % p;
        t[j].push(p);
    }
    let mut pos = 0;

    primes.into_iter().chain((start..).filter(move |&n| {
        let mut is_prime = true;
        while let Some(p) = t[pos].pop() {
            t[(pos + p) % delta].push(p);
            is_prime = false;
        }
        if n == s {
            if is_prime {
                t[(pos + r) % delta].push(r);
                is_prime = false;
            }
            r += 1;
            s = r * r;
        }
        pos = (pos + 1) % delta;
        if pos == 0 {
            delta += 2;
            t.extend([vec![], vec![]]);
        }
        is_prime
    }))
}
pub fn test(limit: usize) -> Vec<usize> {
    let end = limit + 1;
    let mut primes = Vec::with_capacity((limit as f64 / (limit as f64).log(3.0)) as usize);
    primes.extend([2, 3, 5]);

    let bitmap_size = (((limit / 30) >> 3) + 1) << 3;
    let mut sieve = vec![0u8; bitmap_size].into_boxed_slice();
    let bitmap = sieve.as_mut_ptr();

    const SET0: u64 = (1 << 1)
        | (1 << 13)
        | (1 << 17)
        | (1 << 29)
        | (1 << 37)
        | (1 << 41)
        | (1 << 49)
        | (1 << 53);
    const SET1: u64 = (1 << 7) | (1 << 19) | (1 << 31) | (1 << 43);
    const SET2: u64 = (1 << 11) | (1 << 23) | (1 << 47) | (1 << 59);
    let sqrt_end = ((end as f64).sqrt() + 1.01) as usize;
    for x in 1..((sqrt_end as f64 / 2.0 + 2.01) as usize) {
        let xx4 = 4 * x * x;
        for y in (1..sqrt_end).step_by(2) {
            let n = xx4 + y * y;
            if n >= end {
                break;
            }
            if (SET0 >> (n % 60)) & 1 == 1 {
                unsafe {
                    *bitmap.add(n / 30) ^= 1 << ((n % 30) * 8 / 30);
                }
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
                unsafe {
                    *bitmap.add(n / 30) ^= 1 << ((n % 30) * 8 / 30);
                }
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
                unsafe {
                    *bitmap.add(n / 30) ^= 1 << ((n % 30) * 8 / 30);
                }
            }
        }
    }
    let mut n = 1;
    for &incr in WHEEL_2_3_5.iter().cycle() {
        if n > sqrt_end {
            break;
        }
        unsafe {
            if *bitmap.add(n / 30) & (1 << (((n % 30) * 8) / 30)) != 0 {
                let n_squared = n * n;
                let mut c = n_squared;
                for &incr in WHEEL_2_3_5.iter().cycle() {
                    if c > end {
                        break;
                    }

                    *bitmap.add(c / 30) &= !(1 << ((c % 30) * 8 / 30));

                    c += incr as usize * n_squared;
                }
            }
        }
        n += incr as usize;
    }
    n = 1;
    for &incr in WHEEL_2_3_5.iter().cycle() {
        n += incr as usize;
        if n > end {
            break;
        }
        unsafe {
            if *bitmap.add(n / 30) & (1 << (((n % 30) * 8) / 30)) != 0 {
                primes.push(n);
            }
        }
    }
    primes
}
pub fn main() {
    use std::time::Instant;
    const N: usize = 1e6 as usize + 7;
    //dbg!(hprime_sieve(500));
    assert_eq!(sieve_of_atkin(210), segmented_sieve_of_eratosthenes(210));
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
    dbg!(rolling_sieve(N).len());
    let end = start.elapsed();
    println!("{:?}", end);

    /* let mut num = 1;
    for (ind, &i) in WHEEL_2_3_5_7.iter().enumerate() {
        println!("{ind}: {} % 30 = {}, ", num, num % 30);
        num += i as usize;
    } */
}
