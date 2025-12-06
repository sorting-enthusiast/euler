// all multiples of 10 have 0 as a root
// all rational roots p/q must have p be a divisor of the lsd and q a divisor of the msd
// this limits the set of possible roots to -9..=-1
// when the gcd of the first and last digits is 1, the only possible root is -1

pub fn main() {
    dbg!(eval(&convert(1012), -1));
    dbg!(eval(&convert(2024), -2));
    dbg!(eval(&convert(1012 * 6), -1));
}
fn convert(mut n: u64) -> [u8; 16] {
    let mut res = [0; 16];
    for i in 0..16 {
        res[i] = (n % 10) as u8;
        n /= 10;
    }
    res
}
fn eval(p: &[u8; 16], x: i64) -> i64 {
    //dbg!(p);
    let mut res = 0;
    let mut pow = 1;
    for &c in p {
        res += i64::from(c) * pow;
        pow *= x;
    }
    res
}
