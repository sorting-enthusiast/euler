/*
let C(n) be the number of squarefree integers of the form x^2+1 such that 1 <= x <= n.
For example, C(10) = 9 and C(1000) = 895.
Find C(123567101113).
*/
const N: usize = 123_567_101_113;
struct BitArray {
    bits: Box<[u64]>,
}
impl BitArray {
    fn zeroed(n: usize) -> BitArray {
        BitArray {
            bits: vec![0; (n >> 6) + 1].into_boxed_slice(),
        }
    }
    fn set(&mut self, i: usize) {
        self.bits[i >> 6] |= 1 << (i & 63);
    }
    fn get(&self, i: usize) -> bool {
        self.bits[i >> 6] & (1 << (i & 63)) != 0
    }
    fn clear(&mut self, i: usize) {
        self.bits[i >> 6] &= !(1 << (i & 63));
    }
}
fn find_primes(n: usize) -> Vec<usize> {
    let mut primes = Vec::new();
    let mut sieve = BitArray::zeroed(n);
    for i in 2..n {
        if !sieve.get(i - 1) {
            primes.push(i);
            for j in 1..=(n / i) {
                sieve.set(i * j - 1);
            }
        }
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
pub fn main() {
    //let range = 1..N;
    //println!("{}",(N as f64).sqrt().ceil() as usize);
    //let primes = find_primes((N as f64).sqrt().ceil() as usize);
    //dbg!(&primes);
    //dbg!(primes.len());
    for i in 1..=5 {
        println!("{}", i * i + 1);
    }
}
