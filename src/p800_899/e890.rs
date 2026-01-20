const N: usize = 7usize.pow(7);
const MOD: usize = 1e9 as usize + 7;
// gf F(x) = \prod_{k\ge 0} {\frac1{1-x^{2^k}}}
// (1-x)F(x) = F(x^2)
// n even: [x^n]F(x) = [x^{n-1}]F(x) + [x^{\frac{n}2}]F(x)
// n odd: [x^n]F(x) = [x^{n-1}]F(x)
// maybe compute using the following algorithm: https://math.stackexchange.com/a/4944090
pub fn main() {
    let mut gf = vec![0; N + 1];
    gf[0] = 1;
    /* for i in 0..=N.ilog2() {
        let step = 1 << i;
        for a in step..=N {
            gf[a] += gf[a - step];
            if gf[a] >= MOD {
                gf[a] -= MOD;
            }
        }
    } */
    /* for i in 1..=N {
        gf[i] = gf[i - 1];
        if i & 1 == 0 {
            gf[i] += gf[i / 2];
            if gf[i] >= MOD {
                gf[i] -= MOD;
            }
        }
    } */
    for i in 1..=N / 2 {
        gf[i] = gf[i - 1];
        gf[i] += gf[i / 2];
        if gf[i] >= MOD {
            gf[i] -= MOD;
        }
    }
    dbg!(gf[N / 2], N);
}
const fn powmod(mut x: usize, mut exp: usize) -> usize {
    if exp == 0 {
        return 1;
    }
    let mut r = 1;
    x %= MOD;
    while exp > 1 {
        if exp & 1 == 1 {
            r = (r * x) % MOD;
        }
        x = (x * x) % MOD;
        exp >>= 1;
    }
    (r * x) % MOD
}
const fn modinv(x: usize) -> usize {
    powmod(x, MOD - 2)
}
fn solve() {
    /* const LG2: usize = 2181;
    let mut cur = 0;
    let mut poly = [vec![1], vec![]];
    let mut eval = vec![];
    for i in 0..LG2 {
        let nv = poly[cur].len();
    } */
}
