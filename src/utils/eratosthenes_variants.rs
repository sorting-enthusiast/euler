use super::bit_array::BitArray;

const WHEEL_2_3_5_7: [u8; 48] = [
    10, 2, 4, 2, 4, 6, 2, 6, 4, 2, 4, 6, 6, 2, 6, 4, 2, 6, 4, 6, 8, 4, 2, 4, 2, 4, 8, 6, 4, 6, 2,
    4, 6, 2, 6, 6, 4, 2, 4, 6, 2, 6, 4, 2, 4, 2, 10, 2,
];
const WHEEL_2_3_5: [u8; 8] = [6, 4, 2, 4, 2, 4, 6, 2];

const MOD30_TO_MASK: [u8; 30] = [
    0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 4, 4, 8, 8, 8, 8, 16, 16, 32, 32, 32, 32, 64, 64, 64, 64, 64,
    64, 128,
];
#[must_use]
pub fn _wheel_factorized_sieve_of_eratosthenes(n: usize) -> Vec<usize> {
    let mut primes = Vec::with_capacity((n as f64 / (n as f64).log(3.0)) as usize);
    let bitmap_size = (((n / 30) >> 3) + 1) << 3;
    let mut sieve = vec![0u8; bitmap_size].into_boxed_slice();
    let bitmap = sieve.as_mut_ptr();
    let mut num = 1;
    primes.extend([2, 3, 5]);
    let mut wheel_incr = WHEEL_2_3_5.iter().cycle();
    loop {
        num += unsafe { *wheel_incr.next().unwrap_unchecked() } as usize;
        if num * num > n {
            break;
        }
        unsafe {
            if *bitmap.add(num / 30) & MOD30_TO_MASK[num % 30] == 0 {
                primes.push(num);
                let mut multiple = num * num;
                for &incr in wheel_incr.clone() {
                    if multiple > n {
                        break;
                    }
                    *bitmap.add(multiple / 30) |= MOD30_TO_MASK[multiple % 30];
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
            if *bitmap.add(num / 30) & MOD30_TO_MASK[num % 30] == 0 {
                primes.push(num);
            }
        }
        num += incr as usize;
    }

    primes
}
#[must_use]
pub fn sieve_of_eratosthenes(n: usize) -> Vec<usize> {
    let mut primes = Vec::with_capacity((n as f64 / (n as f64).log(3.0)) as usize);
    primes.push(2);
    let mut sieve = BitArray::zeroed(n.div_ceil(2));
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
#[must_use]
pub fn simple_wheel_factorized_sieve_of_eratosthenes(n: usize) -> Vec<usize> {
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
    let mut sieve = BitArray::zeroed(n.div_ceil(2));

    let mut num = 1;
    let mut wheel_incr = WHEEL_2_3_5_7.iter().cycle();
    primes.extend([2, 3, 5, 7]);

    loop {
        num += unsafe { *wheel_incr.next().unwrap_unchecked() } as usize;
        if num * num > n {
            break;
        }
        if !sieve.get(num >> 1) {
            primes.push(num);
            let mut multiple = num * num;
            for &incr in wheel_incr.clone() {
                if multiple > n {
                    break;
                }
                sieve.set(multiple >> 1);
                multiple += incr as usize * num;
            }
        }
    }
    for &incr in wheel_incr {
        if num > n {
            break;
        }
        if !sieve.get(num >> 1) {
            primes.push(num);
        }
        num += incr as usize;
    }

    primes
}
