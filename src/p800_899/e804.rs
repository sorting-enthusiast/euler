// counting repr's of integers in Q(sqrt(-163))
// simpler to just use lattice point counting in an ellipse :(
pub fn main() {
    const chi: [i64; 163] = {
        let mut ret = [0; _];
        let mut i = 1;
        while i < 163 {
            ret[i] = legendre_symbol(i) as i64;
            i += 1;
        }
        ret
    };
    const chi_sums: [i64; 163] = {
        let mut ret = [0; _];
        let mut i = 1;
        while i < 163 {
            ret[i] = ret[i - 1] + chi[i];
            i += 1;
        }
        ret
    };
    const N: usize = 1e16 as _;
    const SQRT_N: usize = N.isqrt();
    let mut sum = 0;
    for d in 1..=SQRT_N {
        sum += chi[d % 163] * (N / d) as i64;
        sum += chi_sums[(N / d) % 163];
    }
    sum -= SQRT_N as i64 * chi_sums[SQRT_N % 163];
    dbg!(2 * sum);
}
const fn legendre_symbol(a: usize) -> i32 {
    let mut base = a % 163;
    let mut exp = 81;
    let mut ls = 1;
    while exp > 0 {
        if exp & 1 == 1 {
            ls = (ls * base) % 163;
        }
        base = (base * base) % 163;
        exp >>= 1;
    }
    match ls {
        0 => 0,
        1 => 1,
        _ => -1,
    }
}
