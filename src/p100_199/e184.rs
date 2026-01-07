fn r2(t: i64) -> i64 {
    let mut res = 0;
    let tsqrt = t.isqrt();
    for k in 1..=tsqrt {
        let tk = t / k;
        res += [0, 1, 1, 0][(tk & 3) as usize];
        res += [0, 1, 0, -1][(k & 3) as usize] * tk;
    }
    res -= [0, 1, 1, 0][(tsqrt & 3) as usize] * tsqrt;
    (res << 2) | 1
}
pub fn main() {
    dbg!(r2(9 - 1));
}
