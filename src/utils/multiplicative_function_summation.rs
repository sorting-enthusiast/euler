use crate::utils::{FIArray::FIArrayI64, bit_array::BitArray};

pub const fn sum_n<const MOD: i64>(x: i64) -> i64 {
    let x = x % (MOD << 1);
    (if x & 1 == 0 {
        (x / 2) * (x + 1)
    } else {
        ((x + 1) / 2) * x
    }) % MOD
}
pub fn totient_sieve(n: usize) -> Vec<i64> {
    let mut res = vec![0; n];
    if n < 2 {
        return res;
    }
    let mut composite = BitArray::zeroed(n);
    let mut primes = vec![];
    res[1] = 1;
    for i in 2..n {
        if !composite.get(i) {
            if (i << 1) < n {
                primes.push(i);
            }
            res[i] = i as i64 - 1;
        }
        for &p in &primes {
            if i * p >= n {
                break;
            }
            composite.set(i * p);
            if i % p == 0 {
                res[i * p] = res[i] * p as i64;
                break;
            }
            res[i * p] = res[i] * (p as i64 - 1);
        }
    }
    res
}
pub fn mobius_sieve(n: usize) -> Vec<i64> {
    let mut res = vec![0; n];
    if n < 2 {
        return res;
    }
    let mut composite = BitArray::zeroed(n);
    let mut primes = vec![];
    res[1] = 1;
    for i in 2..n {
        if !composite.get(i) {
            primes.push(i);
            res[i] = -1;
        }
        for &p in &primes {
            if i * p >= n {
                break;
            }
            composite.set(i * p);
            if i % p == 0 {
                res[i * p] = 0;
                break;
            }
            res[i * p] = -res[i];
        }
    }
    res
}
pub fn totient_sum<const MOD: i64>(x: i64) -> FIArrayI64 {
    let y = if x > 64 {
        (1e8 as usize).min((x as f64).powf(2. / 3.) as usize >> 1)
    } else {
        x as usize
    };
    let mut small_phi = totient_sieve(y + 1);
    for i in 2..y {
        small_phi[i] += small_phi[i - 1];
        small_phi[i] %= MOD;
    }
    let mut Phi = FIArrayI64::new(x);

    for v in FIArrayI64::keys(x) {
        if v as usize <= y {
            Phi[v] = small_phi[v as usize];
            continue;
        }
        let mut phi_v = sum_n::<MOD>(v);
        let vsqrt = v.isqrt();
        for i in 1..=vsqrt {
            phi_v -= ((small_phi[i as usize] - small_phi[i as usize - 1]) * (v / i)) % MOD;
            phi_v -= Phi[v / i];
            phi_v %= MOD;
        }
        phi_v += Phi[vsqrt] * vsqrt;
        phi_v %= MOD;
        if phi_v < 0 {
            phi_v += MOD;
        }
        Phi[v] = phi_v;
    }
    Phi
}
pub const fn divisor_summatory(x: i64) -> i64 {
    let mut sum = 0;
    let sqrt = x.isqrt();
    let mut n = 1;
    while n <= sqrt {
        sum += x / n;
        n += 1;
    }
    sum <<= 1;
    sum - sqrt * sqrt
}
pub fn totient_sum_alt<const MOD: i64>(x: i64) -> FIArrayI64 {
    let y = if x > 64 {
        (1e8 as usize).min((x as f64).powf(2. / 3.) as usize >> 1)
    } else {
        x as usize
    };
    let mut small_phi = totient_sieve(y + 1);
    for i in 2..y {
        small_phi[i] += small_phi[i - 1];
        small_phi[i] %= MOD;
    }
    let mut Phi = FIArrayI64::new(x);

    for (ind, v) in FIArrayI64::keys(x).enumerate() {
        if v as usize <= y {
            Phi.set(ind, small_phi[v as usize]);
            continue;
        }
        let mut phi_v = sum_n::<MOD>(v);
        let vsqrt = v.isqrt();
        for i in 1..=vsqrt {
            phi_v -= ((small_phi[i as usize] - small_phi[i as usize - 1]) * (v / i)) % MOD;
            phi_v -= Phi[v / i];
            phi_v %= MOD;
        }
        phi_v += Phi[vsqrt] * vsqrt;
        phi_v %= MOD;
        if phi_v < 0 {
            phi_v += MOD;
        }
        Phi.set(ind, phi_v);
    }
    Phi
}
