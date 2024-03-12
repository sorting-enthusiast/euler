use crate::bit_array::BitArray;

const WHEEL_2_3_5_7: [u8; 48] = [
    10, 2, 4, 2, 4, 6, 2, 6, 4, 2, 4, 6, 6, 2, 6, 4, 2, 6, 4, 6, 8, 4, 2, 4, 2, 4, 8, 6, 4, 6, 2,
    4, 6, 2, 6, 6, 4, 2, 4, 6, 2, 6, 4, 2, 4, 2, 10, 2,
];
const WHEEL_2_3_5: [u8; 8] = [6, 4, 2, 4, 2, 4, 6, 2];
pub fn _wheel_factorized_sieve_of_eratosthenes(n: usize) -> Vec<usize> {
    let mut primes = Vec::with_capacity((n as f64 / (n as f64).log(3.0)) as usize);
    let bitmap_size = (((n / 30) >> 3) + 1) << 3;
    let mut sieve = vec![0u8; bitmap_size].into_boxed_slice();
    let bitmap = sieve.as_mut_ptr();
    let mut num = 1;
    primes.extend([2, 3, 5]);
    let mut wheel_incr = WHEEL_2_3_5.iter().cycle();
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
pub fn sieve(n: usize) -> Vec<usize> {
    let mut p = 2;
    let mut sieve = BitArray::zeroed(n + 1);
    let mut S = vec![0u64; n + 2].into_boxed_slice();
    let mut i = 1;
    while i < n + 1 {
        i += 1;
        S[i] = 1;
    }

    while p * p <= n {
        let mut q = p;
        while p * q <= n {
            let mut x = p * q;
            while x <= n {
                sieve.set(x);
                x *= p;
            }
            q = {
                let k = q;
                let mut t = k + S[k] as usize;
                while sieve.get(t) {
                    S[k] += S[t];
                    t += S[t] as usize;
                }
                t
            };
        }
        p = {
            let k = p;
            let mut t = k + S[k] as usize;
            while sieve.get(t) {
                S[k] += S[t];
                t += S[t] as usize;
            }
            t
        };
    }
    drop(S);
    let mut primes = Vec::with_capacity((n as f64 / (n as f64).log(3.0)) as usize);
    primes.extend((2..=n).filter(|&x| !sieve.get(x)));
    primes
}
pub fn simple_wheel_factorized_sieve_of_eratosthenes(n: usize) -> Vec<usize> {
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
    let mut sieve = BitArray::zeroed((n + 1) / 2);

    let mut num = 1;
    let mut wheel_incr = WHEEL_2_3_5_7.iter().cycle();
    primes.extend([2, 3, 5, 7]);

    loop {
        num += *wheel_incr.next().unwrap() as usize;
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
