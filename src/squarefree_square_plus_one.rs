/*
let C(n) be the number of squarefree integers of the form x^2+1 such that 1 <= x <= n.
For example, C(10) = 9 and C(1000) = 895.
Find C(123567101113).
*/
const N: usize = 123_567_101_113;
pub struct BitArray {
    bits: Box<[u64]>,
}
impl BitArray {
    pub fn zeroed(n: usize) -> BitArray {
        BitArray {
            bits: vec![0; (n >> 6) + 1].into_boxed_slice(),
        }
    }
    pub fn zero(&mut self) {
        self.bits.iter_mut().for_each(|x| *x = 0);
    }
    pub fn set(&mut self, i: usize) {
        self.bits[i >> 6] |= 1 << (i & 63);
    }
    pub fn get(&self, i: usize) -> bool {
        self.bits[i >> 6] & (1 << (i & 63)) != 0
    }
    pub fn clear(&mut self, i: usize) {
        self.bits[i >> 6] &= !(1 << (i & 63));
    }
}
pub fn find_primes(n: usize) -> Vec<usize> {
    let mut primes = vec![2];
    let mut sieve = BitArray::zeroed(n / 2);
    let sqrt_n = ((n as f64).sqrt().ceil() as usize) & !1;
    for num in (3..sqrt_n).step_by(2) {
        if !sieve.get((num - 1) >> 1) {
            primes.push(num);
            for multiple in (num * num..=n).step_by(num << 1) {
                sieve.set((multiple - 1) >> 1);
            }
        }
    }
    //dbg!(&primes);
    /* for i in (sqrt_n..n).step_by(2) {
        if !sieve.get(i >> 1) {
            primes.push(i + 1);
        }
    } */
    primes.extend((sqrt_n..n).step_by(2).filter_map(|i| {
        if sieve.get(i >> 1) {
            None
        } else {
            Some(i + 1)
        }
    }));
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
    let mut ret = vec![2];
    ret.extend((0..k).filter_map(|i| if sieve.get(i) { None } else { Some(2 * i + 3) }));
    ret
}

pub fn segmented_sieve(n: usize) -> Vec<usize> {
    let limit = (n as f64).sqrt().floor() as usize + 1;
    //dbg!(limit);
    let first_primes = find_primes(limit);
    //dbg!(&first_primes);
    let mut primes = Vec::with_capacity(limit << 7);
    primes.extend(first_primes.iter().map(|&i| i));
    let mut low = limit;
    let mut high = limit << 1;
    let mut mark = BitArray::zeroed(limit + 1);
    while low < n {
        if high >= n {
            high = n;
        }
        for &prime in first_primes.iter() {
            let mut lo_lim = (low / prime) * prime;
            if lo_lim < low {
                lo_lim += prime;
            }
            for j in (lo_lim..high).step_by(prime) {
                mark.set(j - low);
            }
        }
        /* let len = (low..high)
            .filter_map(|i| if mark.get(i - low) { None } else { Some(i) })
            .count();
        primes.reserve(len); */
        primes.extend((low..high).filter_map(|i| if mark.get(i - low) { None } else { Some(i) }));
        mark.zero();
        low += limit;
        high += limit;
    }
    primes
}
fn count_squarefrees(n: usize) -> usize {
    let mut square_frees = 0;
    let primes = find_primes((n as f64).sqrt().ceil() as usize);
    for i in primes {
        if i * i + 1 <= n {}
    }
    square_frees
}
use std::time::Instant;

pub fn main() {
    const N: usize = 100_000_000;
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
