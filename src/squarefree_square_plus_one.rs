/*
let C(n) be the number of squarefree integers of the form x^2+1 such that 1 <= x <= n.
For example, C(10) = 9 and C(1000) = 895.
Find C(123567101113).
*/
//const N: usize = 123_567_101_113;
const ELEM_BIT_WIDTH: u8 = (7 - (usize::BITS as u8).leading_zeros()) as u8;
const MASK: usize = usize::BITS as usize - 1;
pub struct BitArray {
    bits: Box<[usize]>,
}
impl BitArray {
    pub fn zeroed(n: usize) -> BitArray {
        BitArray {
            bits: vec![0; (n >> ELEM_BIT_WIDTH) + 1].into_boxed_slice(),
        }
    }
    pub fn zero(&mut self) {
        self.bits.iter_mut().for_each(|x| *x = 0);
    }
    pub fn set(&mut self, i: usize) {
        self.bits[i >> ELEM_BIT_WIDTH] |= 1 << (i & MASK);
    }
    pub fn get(&self, i: usize) -> bool {
        self.bits[i >> ELEM_BIT_WIDTH] & (1 << (i & MASK)) != 0
    }
    pub fn clear(&mut self, i: usize) {
        self.bits[i >> ELEM_BIT_WIDTH] &= !(1 << (i & MASK));
    }
}
pub fn find_primes(n: usize) -> Vec<usize> {
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

pub fn segmented_sieve(n: usize) -> Vec<usize> {
    let limit = (n as f64).sqrt().floor() as usize + 1;
    let mut primes = Vec::with_capacity((n as f64 / (n as f64).log(3.0)) as usize);
    assert_eq!(primes.capacity(), (n as f64 / (n as f64).log(3.0)) as usize);

    primes.push(2);

    let mut sieve = BitArray::zeroed((limit + 1) / 2);
    let sqrt_limit = ((limit as f64).sqrt().ceil() as usize) & !1;

    for num in (3..=sqrt_limit).step_by(2) {
        if !sieve.get(num >> 1) {
            primes.push(num);
            for multiple in (num * num..=limit).step_by(num << 1) {
                sieve.set(multiple >> 1);
            }
        }
    }
    primes.extend(
        ((sqrt_limit + 1)..=limit)
            .step_by(2)
            .filter(|&i| !sieve.get(i >> 1)),
    );

    let first_primes_count = primes.len();

    let mut low = limit;
    let mut high = limit << 1;
    dbg!((limit + 1) / 2);
    while low < n {
        sieve.zero();
        if high >= n {
            high = n;
        }
        for &prime in &primes[1..first_primes_count] {
            let mut lo_lim = if low % prime == 0 {
                low
            } else {
                low - low % prime + prime
            };
            if (lo_lim / prime) & 1 == 0 {
                lo_lim += prime;
            }
            for i in (lo_lim..high).step_by(prime << 1) {
                sieve.set((i - low) >> 1);
            }
        }
        primes.extend(
            (low | 1..high)
                .step_by(2)
                .filter(|&i| !sieve.get((i - low) >> 1)),
        );
        low += limit;
        high += limit;
    }
    primes
}
use std::time::Instant;

pub fn main() {
    dbg!(segmented_sieve(256));
    const N: usize = 1e6 as usize;
    let start = Instant::now();
    dbg!(segmented_sieve(N).len());
    let end = start.elapsed();
    println!("{:?}", end);
    let start = Instant::now();
    dbg!(find_primes(N).len());
    let end = start.elapsed();
    println!("{:?}", end);
    let start = Instant::now();
    dbg!(sundaram_sieve(N).len());
    let end = start.elapsed();
    println!("{:?}", end);
}
