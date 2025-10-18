use itertools::Itertools;

use crate::utils::{gaussian_ints::GaussianI64, prime_sieves::sift};

const fn powmod(mut x: u64, mut exp: u64, modulo: u64) -> u64 {
    if exp == 0 {
        return 1;
    }
    let mut r = 1;
    while exp > 1 {
        if exp & 1 == 1 {
            r = (r * x) % modulo;
        }
        x = (x * x) % modulo;
        exp >>= 1;
    }
    (r * x) % modulo
}
fn sum_of_squares(p: u64) -> (i64, i64) {
    assert_eq!(
        p & 3,
        1,
        "cannot find divisor, prime is not equivalent to 1 mod 4"
    );
    let r = (2..p)
        .map(|i| powmod(i, (p - 1) >> 2, p))
        .find(|r| p - 1 == (r * r) % p)
        .expect("divisor was not found, something is very wrong");
    let mut a = p;
    let mut b = r;

    // Euclidean algorithm
    while b * b > p {
        let tmp = a % b;
        a = b;
        b = tmp;
    }

    let tmp = (p - b * b).isqrt() as i64;
    let b = b as i64;
    if tmp > b { (b, tmp) } else { (tmp, b) }
}

struct Context {
    pub stack: Vec<GaussianI64>,
}
impl Context {
    fn new() -> Self {
        Self { stack: vec![] }
    }
    fn helper(&self, i: usize, prod: GaussianI64) -> u64 {
        if i >= self.stack.len() {
            let GaussianI64 { re: a, im: b } = prod;
            return a.unsigned_abs().min(b.unsigned_abs());
        }
        let p = self.stack[i];
        self.helper(i + 1, prod) + self.helper(i + 1, p.conj() * (prod / p))
    }
    fn dfs(&mut self, prod: GaussianI64, primes: &[GaussianI64]) -> u64 {
        let mut sum = self.helper(0, prod) >> 1;
        for (i, &p) in primes.iter().enumerate() {
            self.stack.push(p);
            sum += self.dfs(p * prod, &primes[i + 1..]);
            self.stack.pop();
        }
        sum
    }
}
// no way is this 70%
pub fn main() {
    let start = std::time::Instant::now();

    let primes = sift(150)
        .into_iter()
        .filter(|&p| p & 3 == 1)
        .map(sum_of_squares)
        .map(GaussianI64::from)
        .collect_vec();
    let res = Context::new().dfs(GaussianI64::from((1, 0)), &primes);
    let end = start.elapsed();
    println!("res = {res}, took {end:?}");
}
// slower, uses way more memory, example of approach taken by others on the forum
fn alt() {
    let start = std::time::Instant::now();

    let primes = sift(150)
        .into_iter()
        .filter_map(|p| (p & 3 == 1).then(|| sum_of_squares(p)))
        .collect_vec();
    let mut solutions = vec![primes[0]];
    let mut buffer = vec![];
    for (a, b) in primes.into_iter().skip(1) {
        for &(u, v) in &solutions {
            buffer.push((u, v));
            buffer.push((u * a + v * b, u * b - v * a));
            buffer.push((u * a - v * b, u * b + v * a));
        }
        buffer.push((a, b));
        core::mem::swap(&mut solutions, &mut buffer);
        buffer.clear();
    }
    let res = solutions
        .into_iter()
        .map(|(a, b)| a.unsigned_abs().min(b.unsigned_abs()))
        .sum::<u64>();
    let end = start.elapsed();
    println!("res = {res}, took {end:?}");
}
