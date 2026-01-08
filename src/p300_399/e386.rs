use crate::utils::primes::{
    count_signatures::{count_signature, iroot},
    primecount::lucy_strengthreduce_alt,
    wheel_sieve,
};
const N: usize = 1e8 as _;
const MAX_LEN: usize = N.ilog2() as usize + 1;
fn polymulmod(a: &[usize; MAX_LEN], b: &[usize; MAX_LEN], c: &mut [usize; MAX_LEN]) {
    c.fill(0);
    for i in 0..MAX_LEN {
        for j in 0..MAX_LEN - i {
            c[i + j] += a[i] * b[j];
        }
    }
}
// all numbers up to 10^8 have at most 8 distinct prime factors

pub fn main() {
    let primes = wheel_sieve(30);
    let mut prod = 1;
    for (i, p) in primes.iter().enumerate() {
        prod *= p;
        println!("{i}: {prod}");
    }
    dbg!(iroot(N, 10));
    let mut a = [0; MAX_LEN];
    a[0] = 1;
    a[1] = 1;
    //let b = a;
    let mut c = [0; MAX_LEN];
    polymulmod(&a, &a, &mut c);
    polymulmod(&c, &c, &mut a);
    polymulmod(&a, &a, &mut c);
    dbg!(c.iter().max());
    dbg!(c);
}
