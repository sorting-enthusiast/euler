use crate::squarefree_square_plus_one::BitArray;
use std::{
    arch::x86_64::{_mm_clmulepi64_si128, _mm_cvtsi64_si128, _mm_extract_epi64},
    collections::HashMap,
    time::Instant,
};
const fn xor_product(x: u64, y: u64) -> u64 {
    let mut res = 0;
    let mut a = if x < y { x } else { y };
    let mut b = if x > y { x } else { y };
    a >>= a.trailing_zeros();
    while a > 0 {
        res ^= b;
        a >>= 1;
        let tmp = a.trailing_zeros();
        a >>= tmp;
        b <<= 1 + tmp;
    }
    res
}
#[inline(always)]
fn xor_mul(a: u64, b: u64) -> u64 {
    unsafe {
        _mm_extract_epi64(
            _mm_clmulepi64_si128(_mm_cvtsi64_si128(a as i64), _mm_cvtsi64_si128(b as i64), 0),
            0,
        ) as u64
    }
}
#[inline(never)]
fn xor_power(base: u64, mut exp: u64) -> u64 {
    unsafe {
        let mut x = _mm_cvtsi64_si128(base as i64);
        let mut y = _mm_cvtsi64_si128(1);
        while exp > 1 {
            if exp & 1 == 1 {
                y = _mm_clmulepi64_si128(x, y, 0);
            }
            x = _mm_clmulepi64_si128(x, x, 0);
            exp >>= 1;
        }
        _mm_extract_epi64(_mm_clmulepi64_si128(x, y, 0), 0) as u64
    }
}

fn sieve_of_eratosthenes_xor_primes(limit: u64) -> Option<u64> {
    let mut sieve = BitArray::zeroed(((limit + 1) >> 1) as usize);
    let new_lim = 1 << deg(limit) as u64;
    let mut counter = 1;
    for num in (3..=new_lim).step_by(2) {
        if !sieve.get((num >> 1) as usize) {
            counter += 1;
            if counter == 5_000_000 {
                return Some(num);
            }

            for i in (num..S(limit, num)).step_by(2) {
                let multiple = xor_mul(num, i);
                if multiple > new_lim {
                    continue;
                }
                sieve.set((multiple >> 1) as usize);
            }
        }
    }
    None
}
#[inline(always)]
const fn deg(a: u64) -> i64 {
    (u64::BITS - 1) as i64 - a.leading_zeros() as i64
}
#[inline]
const fn rem(mut n: u64, a: u64) -> u64 {
    let p = deg(n);
    let q = deg(a);
    let mut i = p - q;
    while i >= 0 {
        if n >= 1 << (i + q) {
            n ^= a << i;
        }
        i -= 1;
    }
    n
}
//len([b for b in 1..=n if xor_mul(a,b) <= n])
#[inline(always)]
const fn S(n: u64, a: u64) -> u64 {
    (n >> deg(a)) - ((n ^ rem(n, a) > n) as u64)
}
fn mu(n: i64) -> i64 {
    let mut r = (n == 1) as i64;
    for k in 1..n {
        if n % k == 0 {
            r -= mu(k);
        }
    }
    r
}
//adapted from stimmer's solution, with minor changes
fn find_nth_xorprime(T: u64) -> u64 {
    let mut pp = 0;
    let mut m = 1;
    loop {
        let mut a = 0;
        for d in 1..=m {
            if m % d == 0 {
                a += mu(d) << (m / d);
            }
        }
        if pp + a / m >= T as i64 {
            break;
        }
        pp += a / m;
        m += 1;
    }
    let h = (m + 2) / 2;
    let limit: u64 = 1 << h;
    let mut sieve = BitArray::zeroed(limit as usize);
    let mut primes = Vec::new();
    for x in 2..limit {
        if !sieve.get(x as usize) {
            primes.push(x);
            for y in x..(limit >> deg(x)) {
                sieve.set(xor_mul(x, y) as usize);
            }
        }
    }
    let k = (m + 1) / 4;
    let mut d = HashMap::new();
    d.insert(1, 1);
    for &p in primes.iter().rev() {
        let e = d.clone();
        for (q, t) in e {
            if deg(p) + deg(q) > m {
                break;
            }
            let mut r = xor_mul(p, q);
            if deg(r) >= k {
                r &= ((1 << k) - 1) << (deg(r) - k + 1);
            }
            d.entry(r).and_modify(|v| *v -= t);
        }
    }
    let w = 1 << (m - k + 1);
    let mut s = 1 << m;
    while s < (1 << (m + 1)) {
        let n = s + w - 1;
        let mut z = 0;
        for (&x, &t) in &d {
            z += (S(n, x) as i64) * t;
        }
        if z + primes.len() as i64 > T as i64 {
            break;
        }
        pp = z + primes.len() as i64 - 1;
        s += w;
    }
    let ss = 1 << h;
    let ds = deg(ss);
    let mut l = BitArray::zeroed(ss as usize);
    loop {
        for &p in primes.iter().rev() {
            let dp = deg(p);
            let mut i = rem(s, p);
            for c in 1..=(1 << (ds - dp)) {
                l.set(i as usize);
                i ^= p * (((c ^ (c - 1)) + 1) >> 1);
            }
        }
        for c in 0..ss {
            if !l.get(c as usize) {
                pp += 1;
                if pp == T as i64 {
                    return s + c;
                }
            }
        }
        l.zero();
        s += ss;
    }
}
pub fn main() {
    dbg!(xor_power(11, 20));
    dbg!(xor_product(3, 3));
    dbg!(xor_mul(3, 3));
    let now = Instant::now();
    dbg!(find_nth_xorprime(5_000_000));
    let elapsed = now.elapsed();
    println!("{:?}", elapsed);
    let now = Instant::now();
    dbg!(sieve_of_eratosthenes_xor_primes(150_000_000));
    let elapsed = now.elapsed();
    println!("{:?}", elapsed);
    dbg!((!0u64) as i64);
}
