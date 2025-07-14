use crate::utils::{bit_array::BitArray, multiplicative_function_summation::totient_sieve};

pub fn main() {
    const LIM: usize = 2e10 as usize;
    const BLOCK: usize = LIM.isqrt() + 1;

    let mut sum = 0;
    let mut phis = vec![0; BLOCK];
    let mut composite = BitArray::zeroed(BLOCK);
    let mut primes = vec![];
    phis[1] = 1;
    for i in 2..BLOCK {
        if !composite.get(i) {
            if i * i <= LIM {
                primes.push(i);
            }
            phis[i] = i - 1;
        }
        for &p in &primes {
            let ip = i * p;
            if ip >= BLOCK {
                break;
            }
            composite.set(ip);
            if i % p == 0 {
                phis[ip] = phis[i] * p;
                break;
            }
            phis[ip] = phis[i] * (p - 1);
        }
    }
    let mut count = 0;
    for i in 2..BLOCK {
        let cophi = i - phis[i];
        if cophi != 1 && (i - 1) % cophi == 0 {
            println!("{i}: {}, {}", i - 1, cophi);
            sum += i;
            count += 1;
        }
    }
    dbg!(sum);

    let primes = primes;
    dbg!(primes.len());
    let mut l = BLOCK;
    let mut large_prime_factor = vec![0; BLOCK];
    while l < LIM {
        let r = LIM.min(l + BLOCK);

        phis.iter_mut().enumerate().for_each(|(i, e)| *e = i + l);
        large_prime_factor
            .iter_mut()
            .enumerate()
            .for_each(|(i, e)| *e = i + l);
        for &p in &primes {
            if p * p >= r {
                break;
            }
            let mut multiple = l.next_multiple_of(p);
            while multiple < r {
                //dbg!((p, multiple));
                phis[multiple - l] /= p;
                phis[multiple - l] *= p - 1;

                let mut k = large_prime_factor[multiple - l];
                //dbg!(k);

                while k % p == 0 {
                    k /= p;
                }
                large_prime_factor[multiple - l] = k;
                multiple += p;
            }
        }
        for i in l..r {
            let k = large_prime_factor[i - l];
            if k != 1 {
                phis[i - l] /= k;
                phis[i - l] *= k - 1;
            }
            let cophi = i - phis[i - l];
            if cophi != 1 && (i - 1) % cophi == 0 {
                println!("{i}: {}, {}", i - 1, cophi);
                sum += i;
                count += 1;
            }
        }
        dbg!((l, r, sum, count));
        l = r;
    }
    dbg!(sum);
    /* sum = 0;
    let phis = totient_sieve(LIM + 1);
    for i in 2..=LIM {
        let cophi = i - phis[i] as usize;
        if cophi != 1 && (i - 1) % cophi == 0 {
            //println!("{i}: {}, {}", i - 1, cophi);
            sum += i;
        }
    }
    dbg!(sum); */
}
