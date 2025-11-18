use std::ptr::write_bytes;
use std::vec;

use itertools::Itertools;

//const N: usize = 123_567_101_113;
use super::bit_array::BitArray;
use super::eratosthenes_variants::sieve_of_eratosthenes;
pub use super::sieve_of_pritchard::{BIT64TOVAL240, sift};
pub const WHEEL_2_3_5_7: [u8; 48] = [
    10, 2, 4, 2, 4, 6, 2, 6, 4, 2, 4, 6, 6, 2, 6, 4, 2, 6, 4, 6, 8, 4, 2, 4, 2, 4, 8, 6, 4, 6, 2,
    4, 6, 2, 6, 6, 4, 2, 4, 6, 2, 6, 4, 2, 4, 2, 10, 2,
];
pub const WHEEL_2_3_5: [u8; 8] = [6, 4, 2, 4, 2, 4, 6, 2];
pub const WHEEL_2_3_5_VALS: [u8; 8] = [1, 7, 11, 13, 17, 19, 23, 29];
pub const MOD30_TO_BIT8: [u64; 30] = [
    !0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7,
];
pub const MOD30_TO_MASK: [u8; 30] = [
    0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 4, 4, 8, 8, 8, 8, 16, 16, 32, 32, 32, 32, 64, 64, 64, 64, 64,
    64, 128,
];
#[inline(always)]
#[must_use]
pub const fn num2bits(num: usize) -> u8 {
    MOD30_TO_BIT8[num % 30] as u8
}
//TODO: fix
/* #[must_use]
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
            pattern[i - 1][r as usize / 30] |= MOD30_TO_MASK[r as usize % 30];
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
            if *bitmap.add(a as usize / 30) & (1 << (a_i & 7)) == 0 {
                *as_ptr.add(cnt) = a;
                for i in 0..8 {
                    *offs.add(
                        8 * cnt + num2bits((a * WHEEL_2_3_5_VALS[i] as u32) as usize) as usize,
                    ) = a * u32::from(WHEEL_2_3_5_VALS[i]) / 30;
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

            a += u32::from(WHEEL_2_3_5_VALS[a_i & 7]);
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
                        *(bmp.cast::<u64>()).add(j) |=
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
            if *bitmap.add(num / 30) & MOD30_TO_MASK[num % 30] == 0 {
                primes.push(num);
            }
        }
    }

    primes
}
 */
#[must_use]
pub fn linear_sieve(n: u64) -> Vec<u64> {
    let mut p = 7;
    let mut primes = Vec::with_capacity((n as f64 / (n as f64).log(3.0)) as usize);
    primes.extend([2, 3, 5]);

    let bitmap_size = (((n as usize / 30) >> 3) + 1) << 3;
    let mut sieve = vec![0; bitmap_size].into_boxed_slice();
    let bitmap = sieve.as_mut_ptr();
    let mut s = vec![0; bitmap_size << 3].into_boxed_slice();
    for (x, incr) in s.iter_mut().zip(WHEEL_2_3_5.into_iter().cycle()) {
        *x = u16::from(incr);
    }
    while p * p <= n {
        primes.push(p);
        let mut q = p;
        while p * q <= n {
            let mut x = p * q;
            while x <= n {
                unsafe {
                    *bitmap.add(x as usize / 30) |= MOD30_TO_MASK[x as usize % 30];
                }
                x *= p;
            }

            let prev = q;
            q += u64::from(s[((q / 30) * 8 + MOD30_TO_BIT8[q as usize % 30]) as usize]);
            while unsafe { *bitmap.add(q as usize / 30) & MOD30_TO_MASK[q as usize % 30] != 0 } {
                q += u64::from(s[((q / 30) * 8 + MOD30_TO_BIT8[q as usize % 30]) as usize]);
            }
            s[((prev / 30) * 8 + MOD30_TO_BIT8[prev as usize % 30]) as usize] = (q - prev) as u16;
        }

        p += u64::from(s[((p / 30) * 8 + MOD30_TO_BIT8[p as usize % 30]) as usize]);
        while unsafe { *bitmap.add(p as usize / 30) & MOD30_TO_MASK[p as usize % 30] != 0 } {
            p += u64::from(s[((p / 30) * 8 + MOD30_TO_BIT8[p as usize % 30]) as usize]);
        }
    }
    /* while p <= n {
        if unsafe { *bitmap.add(p as usize / 30) & MOD30_TO_MASK[p as usize % 30] == 0 } {
            primes.push(p);
        } else {
            println!("oops");
        }
        p += s[((p / 30) * 8 + MOD30_TO_BIT8[p as usize % 30]) as usize];
    } */
    drop(s);
    for incr in WHEEL_2_3_5
        .into_iter()
        .cycle()
        .skip(MOD30_TO_BIT8[p as usize % 30] as usize)
    {
        if p > n {
            break;
        }
        unsafe {
            if *bitmap.add(p as usize / 30) & MOD30_TO_MASK[p as usize % 30] == 0 {
                primes.push(p);
            }
        }
        p += u64::from(incr);
    }
    primes
}

#[must_use]
pub fn erat_wheel210(n: usize) -> Vec<usize> {
    if n < 7 {
        let mut ret = vec![2];
        if n > 2 {
            ret.push(3);
        }
        if n > 4 {
            ret.push(5);
        }
        return ret;
    }
    let bitmap_size = (((n / 30) >> 3) + 1) << 3;
    let mut sieve = vec![0u8; bitmap_size].into_boxed_slice();
    let bitmap = sieve.as_mut_ptr();
    let mut primes = Vec::with_capacity((n as f64 / (n as f64).log(3.0)) as usize);
    primes.extend([2, 3, 5, 7]);
    let mut num = 1;
    let mut wheel_incr = WHEEL_2_3_5_7.into_iter().cycle();

    loop {
        num += unsafe { wheel_incr.next().unwrap_unchecked() } as usize;
        if num * num > n {
            break;
        }
        if unsafe { *bitmap.add(num / 30) } & MOD30_TO_MASK[num % 30] == 0 {
            primes.push(num);
            let mut multiple = num * num;
            /*for incr in wheel_incr.clone() {
                if multiple > n {
                    break;
                }
                *bitmap.add(multiple / 30) |= MOD30_TO_MASK[num % 30];
                multiple += incr as usize * num;
            } */
            let precomp = [
                0,
                0,
                num << 1,
                0,
                num << 2,
                0,
                (num << 2) + (num << 1),
                0,
                num << 3,
                0,
                (num << 3) + (num << 1),
            ]; //to reduce bit complexity, only use adds and shifts instead of a multiplication every iteration of inner loop
            for incr in wheel_incr.clone() {
                if multiple > n {
                    break;
                }
                unsafe { *bitmap.add(multiple / 30) |= MOD30_TO_MASK[multiple % 30] };
                multiple += precomp[incr as usize];
            }
        }
    }
    for incr in wheel_incr {
        if num > n {
            break;
        }
        if unsafe { *bitmap.add(num / 30) & MOD30_TO_MASK[num % 30] } == 0 {
            primes.push(num);
        }
        num += incr as usize;
    }

    primes
}

#[must_use]
pub fn erat_wheel30(n: usize) -> Vec<usize> {
    if n < 7 {
        let mut ret = vec![2];
        if n > 2 {
            ret.push(3);
        }
        if n > 4 {
            ret.push(5);
        }
        return ret;
    }
    let mut primes = Vec::with_capacity((n as f64 / (n as f64).log(3.0)) as usize);
    let bitmap_size = (((n / 30) >> 3) + 1) << 3;
    let mut sieve = vec![0u8; bitmap_size].into_boxed_slice();
    let bitmap = sieve.as_mut_ptr();
    sieve[0] = 1;
    let mut num = 1;
    let mut wheel_incr = WHEEL_2_3_5.into_iter().cycle();
    primes.extend([2, 3, 5]);

    loop {
        num += unsafe { wheel_incr.next().unwrap_unchecked() } as usize;
        if num * num > n {
            break;
        }
        if unsafe { *bitmap.add(num / 30) } & MOD30_TO_MASK[num % 30] == 0 {
            //primes.push(num);
            let mut multiple = num * num;
            let precomp = [
                0,
                0,
                num << 1,
                0,
                num << 2,
                0,
                (num << 2) + (num << 1),
                /* 0,
                num << 3,
                0,
                num * 10, */
            ]; //to reduce bit complexity, multiply only twice instead of every iteration of inner loop
            for incr in wheel_incr.clone() {
                if multiple > n {
                    break;
                }
                unsafe { *bitmap.add(multiple / 30) |= MOD30_TO_MASK[multiple % 30] };
                multiple += precomp[incr as usize];
            }
        }
    }
    /* for incr in wheel_incr {
        if num > n {
            break;
        }
        unsafe {
            if *bitmap.add(num / 30) & MOD30_TO_MASK[num % 30] == 0 {
                primes.push(num);
            }
        }
        num += incr as usize;
    } */

    let bitmap64: *const u64 = bitmap.cast();
    let mut base = 0;
    for k in 0..bitmap_size / 8 - 1 {
        let mut bitset = unsafe { *bitmap64.add(k) };
        bitset = !bitset;
        while bitset != 0 {
            let r = bitset.trailing_zeros() as usize;
            primes.push(base + BIT64TOVAL240[r] as usize);
            bitset &= bitset - 1;
        }

        base += 240;
    }
    let mut bitset = unsafe { *bitmap64.add(bitmap_size / 8 - 1) };
    bitset = !bitset;
    while bitset != 0 {
        let r = bitset.trailing_zeros() as usize;
        let prime_cand = base + BIT64TOVAL240[r] as usize;
        if prime_cand > n {
            break;
        }
        primes.push(prime_cand);
        bitset &= bitset - 1;
    }
    primes
}
#[must_use]
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
// TODO: fix
#[must_use]
pub fn segmented_sieve_of_eratosthenes(n: usize) -> Vec<usize> {
    if n < 7 {
        let mut ret = vec![2];
        if n > 2 {
            ret.push(3);
        }
        if n > 4 {
            ret.push(5);
        }
        return ret;
    }
    let limit = (n as f64).sqrt().floor() as usize + 1;
    let mut primes = Vec::with_capacity((n as f64 / (n as f64).log(3.0)) as usize);
    let bitmap_size = (((limit / 30) >> 3) + 1) << 3;
    let mut sieve = vec![0u8; bitmap_size].into_boxed_slice();
    let bitmap = sieve.as_mut_ptr();

    let mut num = 1;
    let mut wheel_incr = WHEEL_2_3_5_7.into_iter().cycle();
    primes.extend([2, 3, 5, 7]);

    loop {
        num += unsafe { wheel_incr.next().unwrap_unchecked() } as usize;
        if num * num > limit {
            break;
        }
        unsafe {
            if *bitmap.add(num / 30) & MOD30_TO_MASK[num % 30] == 0 {
                primes.push(num);
                let mut multiple = num * num;
                let precomp = [
                    0,
                    0,
                    num << 1,
                    0,
                    num << 2,
                    0,
                    num * 6,
                    0,
                    num << 3,
                    0,
                    num * 10,
                ];
                for incr in wheel_incr.clone() {
                    if multiple > limit {
                        break;
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
            }
        }
        num += incr as usize;
    }
    let first_primes_count = primes.len();
    let mut low = limit;
    let mut high = limit << 1;

    let mut first_multiple = Vec::with_capacity(first_primes_count);
    for &prime in &primes[3..] {
        let mut multiple = prime * prime;
        let precomp = [0, 0, num << 1, 0, num << 2, 0, num * 6];
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
                multiple = tmp;
                if tmp > low {
                    break;
                }
            }
        }
        first_multiple.push(multiple);
    }

    while low < n {
        unsafe { write_bytes(bitmap, 0, bitmap_size) };
        if high > n {
            high = n;
        }
        for (&prime, multiple) in primes[3..first_primes_count]
            .iter()
            .zip(first_multiple.iter_mut())
        {
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

        for incr in WHEEL_2_3_5.into_iter().cycle() {
            if num <= low {
                num += incr as usize;
                continue;
            }
            if num > high {
                break;
            }
            unsafe {
                if *bitmap.add((num - low) / 30) & MOD30_TO_MASK[num % 30] == 0 {
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
#[must_use]
pub fn sieve_of_atkin(limit: usize) -> Vec<usize> {
    const SET0: u64 = (1 << 1)
        | (1 << 13)
        | (1 << 17)
        | (1 << 29)
        | (1 << 37)
        | (1 << 41)
        | (1 << 49)
        | (1 << 53);
    const SET1: u64 = (1 << 7) | (1 << 19) | (1 << 31) | (1 << 43)
    //|(1<<1)|(1<<13)|(1<<37) |(1<<49)
    ;
    const SET2: u64 = (1 << 11) | (1 << 23) | (1 << 47) | (1 << 59);
    let end = limit + 1;
    let mut primes = Vec::with_capacity((limit as f64 / (limit as f64).ln()) as usize);
    primes.extend([2, 3, 5]);

    let bitmap_size = (((limit / 30) >> 3) + 1) << 3;
    let mut sieve = vec![0u8; bitmap_size].into_boxed_slice();
    let bitmap = sieve.as_mut_ptr();

    let sqrt_end = limit.isqrt() + 1; //((end as f64).sqrt() + 1.01) as usize;
    //n = 4x^2+y^2
    let mut xx4 = 4;
    //((sqrt_end as f64 / 2.0 + 2.01) as usize)
    for x in 1..=(limit >> 1).isqrt() {
        let mut n = xx4 + 1;
        for y in (1..sqrt_end).step_by(2) {
            if n >= end {
                break;
            }
            if (SET0 >> (n % 60)) & 1 == 1 {
                unsafe { *bitmap.add(n / 30) ^= MOD30_TO_MASK[n % 30] };
            }
            n += (y + 1) << 2;
        }
        xx4 += (x << 1 | 1) << 2;
    }
    //n = 3x^2+y^2
    let mut xx3 = 3;
    //(12 * (sqrt_end as f64 / 3.0f64.sqrt() + 2.01) as usize)

    for dxx in (0..=12 * (limit / 3).isqrt()).step_by(24) {
        xx3 += dxx;
        let mut n = xx3 + 4;
        for y in (2..sqrt_end).step_by(2) {
            if n >= end {
                break;
            }
            if (SET1 >> (n % 60)) & 1 == 1 {
                unsafe { *bitmap.add(n / 30) ^= MOD30_TO_MASK[n % 30] };
            }
            n += (y + 1) << 2;
        }
    }
    let mut n_start = 3;
    //(sqrt_end as f64 / 2.0f64.sqrt() + 2.01) as usize
    for x in 1..=(limit >> 1).isqrt() {
        let mut n = n_start;
        for y in (0..x).rev().step_by(2) {
            if n >= end {
                break;
            }
            if (SET2 >> (n % 60)) & 1 == 1 {
                unsafe { *bitmap.add(n / 30) ^= MOD30_TO_MASK[n % 30] };
            }
            n += (y - 1) << 2;
        }
        n_start += (x + 1) << 2;
    }
    let mut n = 1;
    for incr in WHEEL_2_3_5.into_iter().cycle() {
        if n > sqrt_end {
            break;
        }
        if unsafe { *bitmap.add(n / 30) } & MOD30_TO_MASK[n % 30] != 0 {
            let n_squared = n * n;
            let mut c = n_squared;
            let precomp = [0, 0, n_squared << 1, 0, n_squared << 2, 0, n_squared * 6];
            for incr in WHEEL_2_3_5.into_iter().cycle() {
                if c > end {
                    break;
                }
                unsafe { *bitmap.add(c / 30) &= !MOD30_TO_MASK[c % 30] };
                c += precomp[incr as usize];
            }
        }
        n += incr as usize;
    }
    /* n = 1;
    for incr in WHEEL_2_3_5.into_iter().cycle() {
        n += incr as usize;
        if n > end {
            break;
        }
        unsafe {
            if *bitmap.add(n / 30) & MOD30_TO_MASK[n % 30] != 0 {
                primes.push(n);
            }
        }
    } */
    let bitmap64: *const u64 = bitmap.cast();
    let mut base = 0;
    for k in 0..bitmap_size / 8 {
        unsafe {
            let mut bitset = *bitmap64.add(k);
            while bitset != 0 {
                let r = bitset.trailing_zeros() as usize;
                primes.push(base + BIT64TOVAL240[r] as usize);
                bitset &= bitset - 1;
            }
        }
        base += 240;
    }
    primes
}

#[must_use]
pub fn sieve_of_atkin_alt(limit: usize) -> Vec<usize> {
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
    let end = limit + 1;
    let mut primes = Vec::with_capacity((limit as f64 / (limit as f64).ln()) as usize);
    primes.extend([2, 3, 5]);

    let bitmap_size = (((limit / 30) >> 3) + 1) << 3;
    let mut sieve = vec![0u8; bitmap_size].into_boxed_slice();
    let bitmap = sieve.as_mut_ptr();
    let sqrt_end = limit.isqrt();
    let mut xx = 1;
    for x in 1..=sqrt_end {
        let mut yy = 1;
        for y in 1..=sqrt_end {
            let xx3 = (xx << 1) + xx;
            let n = xx3 + xx + yy;
            if n <= limit && (SET0 >> (n % 60)) & 1 == 1 {
                unsafe { *bitmap.add(n / 30) ^= MOD30_TO_MASK[n % 30] };
            }
            let n = xx3 + yy;
            if n <= limit && (SET1 >> (n % 60)) & 1 == 1 {
                unsafe { *bitmap.add(n / 30) ^= MOD30_TO_MASK[n % 30] };
            }
            if y < x {
                let n = xx3 - yy;
                if n <= limit && (SET2 >> (n % 60)) & 1 == 1 {
                    unsafe { *bitmap.add(n / 30) ^= MOD30_TO_MASK[n % 30] };
                }
            }
            yy += y << 1 | 1;
        }
        xx += x << 1 | 1;
    }
    let mut n = 1;
    for incr in WHEEL_2_3_5.into_iter().cycle() {
        if n > sqrt_end {
            break;
        }
        if unsafe { *bitmap.add(n / 30) } & MOD30_TO_MASK[n % 30] != 0 {
            let n_squared = n * n;
            let mut c = n_squared;
            let precomp = [0, 0, n_squared << 1, 0, n_squared << 2, 0, n_squared * 6];
            for incr in WHEEL_2_3_5.into_iter().cycle() {
                if c > end {
                    break;
                }
                unsafe { *bitmap.add(c / 30) &= !MOD30_TO_MASK[c % 30] };
                c += precomp[incr as usize];
            }
        }
        n += incr as usize;
    }
    /* n = 1;
    for incr in WHEEL_2_3_5.into_iter().cycle() {
        n += incr as usize;
        if n > end {
            break;
        }
        unsafe {
            if *bitmap.add(n / 30) & MOD30_TO_MASK[n % 30] != 0 {
                primes.push(n);
            }
        }
    } */
    let bitmap64: *const u64 = bitmap.cast();
    let mut base = 0;
    for k in 0..bitmap_size / 8 {
        unsafe {
            let mut bitset = *bitmap64.add(k);
            while bitset != 0 {
                let r = bitset.trailing_zeros() as usize;
                primes.push(base + BIT64TOVAL240[r] as usize);
                bitset &= bitset - 1;
            }
        }
        base += 240;
    }
    primes
}

#[must_use]
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
    let mut r = (start as f64).sqrt() as usize + 1;
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

// TODO:
// reduce the space required for storing the primes eightfold by saving p_i - p_i-1 instead of p_i, allowing for a Vec<u8> instead of Vec<u64/usize>
// use priority queue in a segmented sieve to store the next multiple of to erase
// size: pi(sqrt(n)) == sqrt(n)/log(sqrt(n))
// decrease_key on root: log(size)
// times called: n - sqrt(n) - pi(n) + pi(sqrt(n)) == n - n/log(n) - (sqrt(n) - sqrt(n)/log(sqrt(n)))
// total complexity:
//

//TODO: fix
pub fn dijkstra_sieve(limit: usize) -> Vec<usize> {
    let mut primes = vec![2, 3, 5];
    let mut multiples = vec![];
    let mut incrs: Vec<std::iter::Cycle<std::array::IntoIter<u8, 8>>> = vec![];
    let mut num = 1;
    let mut wheel_incr = WHEEL_2_3_5.into_iter().cycle();

    loop {
        num += unsafe { wheel_incr.next().unwrap_unchecked() } as usize;
        if num > limit {
            break;
        }
        if let Some(pos) = multiples.iter().position(|&e| e == num) {
            let incr = unsafe { incrs[pos].next().unwrap_unchecked() } as usize;
            multiples[pos] += incr * primes[pos + 3];
        } else {
            primes.push(num);
            if num * num <= limit {
                multiples.push(num * num);
                incrs.push(wheel_incr.clone());
            }
        }
    }
    primes
}

pub fn main() {
    use std::time::Instant;
    const N: usize = 1e9 as usize;
    //dbg!(hprime_sieve(500));
    /* {
         use rand::distr::uniform::{UniformSampler, UniformUsize};
           let mut rng = rand::rng();
           let range = UniformUsize::new_inclusive(100, 1e9 as usize).unwrap();
           for i in 0..32 {
               if i & 7 == 0 {
                   print!("{i} ");
               }
               let n = range.sample(&mut rng);
               assert_eq!(
                   sieve_of_atkin(n),
                   wheel_factorized_sieve_of_eratosthenes(n),
                   "{n}"
               );
           }
           println!();
       }
    */
    /* {
        let mut i = 1;
        dbg!(
            WHEEL_2_3_5_7
                .into_iter()
                .cycle()
                .take(96)
                .filter(move |&incr| {
                    let ret = i & 3 == 1;
                    i += u32::from(incr);
                    ret
                })
                .count()
        );
        println!();
        let mut i = 1;
        dbg!(
            WHEEL_2_3_5_7
                .into_iter()
                .cycle()
                .take(96)
                .filter(move |&incr| {
                    let ret = i % 6 == 1;
                    i += u32::from(incr);
                    ret
                })
                .count()
        );
        println!();
        let mut i = 1;

        dbg!(
            WHEEL_2_3_5_7
                .into_iter()
                .cycle()
                .take(96)
                .filter(move |&incr| {
                    let ret = i % 12 == 11;
                    i += u32::from(incr);
                    ret
                })
                .count()
        );
        println!();

        let mut i = 1;
        dbg!(
            WHEEL_2_3_5
                .into_iter()
                .cycle()
                .take(16)
                .filter(move |&incr| {
                    let ret = i & 3 == 1;
                    i += u32::from(incr);
                    ret
                })
                .count()
        );
        println!();
        let mut i = 1;
        dbg!(
            WHEEL_2_3_5
                .into_iter()
                .cycle()
                .take(16)
                .filter(move |&incr| {
                    let ret = i % 6 == 1;
                    i += u32::from(incr);
                    ret
                })
                .count()
        );
        println!();
        let mut i = 1;

        dbg!(
            WHEEL_2_3_5
                .into_iter()
                .cycle()
                .take(16)
                .filter(move |&incr| {
                    let ret = i % 12 == 11;
                    i += u32::from(incr);
                    ret
                })
                .count()
        );
    }
     */
    assert_eq!(sieve_of_atkin(240), sieve_of_atkin_alt(240));
    assert_eq!(
        sieve_of_atkin(29)
            .into_iter()
            .map(|p| p as u64)
            .collect_vec(),
        sift(29)
    );
    let start = Instant::now();
    dbg!(sift(N as u64).len());
    let end = start.elapsed();
    println!("{end:?}");

    /* let start = Instant::now();
    dbg!(segmented_sieve_of_eratosthenes(N).len());
    let end = start.elapsed();
    println!("{:?}", end); */
    let start = Instant::now();
    dbg!(erat_wheel210(N).len());
    let end = start.elapsed();
    println!("{end:?}");

    let start = Instant::now();
    dbg!(erat_wheel30(N).len());
    let end = start.elapsed();
    println!("{end:?}");

    let start = Instant::now();
    dbg!(dijkstra_sieve(N).len());
    let end = start.elapsed();
    println!("{end:?}");

    let start = Instant::now();
    dbg!(linear_sieve(N as u64).len());
    let end = start.elapsed();
    println!("{end:?}");

    let start = Instant::now();
    dbg!(sieve_of_atkin(N).len());
    let end = start.elapsed();
    println!("{end:?}");
    /* let start = Instant::now();
    dbg!(sieve_of_atkin_alt(N).len());
    let end = start.elapsed();
    println!("{end:?}"); */
    let start = Instant::now();
    dbg!(sieve_of_eratosthenes(N).len());
    let end = start.elapsed();
    println!("{end:?}");

    let start = Instant::now();
    dbg!(sundaram_sieve(N).len());
    let end = start.elapsed();
    println!("{end:?}");

    let start = Instant::now();
    dbg!(rolling_sieve(N).len());
    let end = start.elapsed();
    println!("{end:?}");
}
