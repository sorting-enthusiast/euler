use std::time::Instant;

const MODULO: u64 = 20_092_010;
const LEN: usize = 2000;

// x^2000 = x + 1
fn polymulmod(poly1: &[u64; LEN], poly2: &[u64; LEN], res: &mut [u64; LEN]) {
    res.fill(0);
    for (i, &c1) in poly1.iter().enumerate().filter(|&(_, &c)| c != 0) {
        for (j, &c2) in poly2[..LEN - i].iter().enumerate() {
            res[i + j] += c1 * c2;
        }
        for (j, &c2) in poly2[LEN - i..].iter().enumerate() {
            let c = c1 * c2;
            res[j] += c;
            res[j + 1] += c;
        }
    }
    for c in res.iter_mut() {
        *c %= MODULO;
    }
}
pub fn main() {
    // Notice that the matrix M that describes the recurrence relation satisfies M^2000 = M + I (noticed by trial and error, proven by Cayley-Hamilton theorem)
    // Hence, M^(10^18 - 2000 + 1) = (M^2000)^(5 * 10^15 - 1) * M = M * (M + I)^(5 * 10^15 - 1)
    // Moreover, for every power of M in 1..1999 the sum of the bottom row is 2, while for 0 (the identity) the sum is 1
    // Finding the coefficients of the polynomial x * (x + 1)^(5 * 10^15 - 1) modulo x^2000 - x - 1 over the natural numbers modulo 20092010,
    // and calculating their weighted sum solves the question
    let start = Instant::now();

    let mut exp = 1e18 as usize / LEN - 1;

    let mut x_store = [0; LEN];
    let mut r_store = [0; LEN];
    let mut tmp_store = [0; LEN];
    let (mut x, mut r, mut tmp) = (&mut x_store, &mut r_store, &mut tmp_store);
    x[0] = 1;
    x[1] = 1;
    r[0] = 1;
    while exp > 1 {
        if exp & 1 == 1 {
            polymulmod(r, x, tmp);
            std::mem::swap(&mut tmp, &mut r);
            //r = (r * x) % MODULO;
        }
        polymulmod(x, x, tmp);
        std::mem::swap(&mut tmp, &mut x);
        //x = (x * x) % MODULO;
        exp >>= 1;
    }
    polymulmod(r, x, tmp);

    let sum = (2 * (*tmp).into_iter().sum::<u64>() + tmp[LEN - 1]) % MODULO;
    println!("{:?}", start.elapsed());
    println!("{sum}");
    solve_ext();
}
// https://arxiv.org/pdf/2008.08822
fn solve_ext() {
    // use num_bigint::BigUint;
    fn conv(a: &[u64; 2001], b: &[u64; 2001], c: &mut [u64; 4001]) {
        c.fill(0);
        for (i, &ai) in a.iter().enumerate().filter(|&(_, &v)| v != 0) {
            for (j, &bj) in b.iter().enumerate()
            /* .filter(|&(_, &v)| v != 0) */
            {
                c[i + j] += ai * bj;
            }
        }
        for ci in c.iter_mut() {
            *ci %= MODULO;
        }
    }

    let start = Instant::now();
    let mut p = [1; 2001]; // technically never goes above degree 2000, but nvm
    let mut q = [0; 2001];
    let mut tmp = [0; 4001];

    p[1999] = 0;
    p[2000] = 0;

    q[0] = 1;
    q[1999] = MODULO - 1;
    q[2000] = MODULO - 1;
    let mut n = 1e18 as usize; //BigUint::from(10u32).pow(9 << 13);
    //dbg!(n.bits());
    while n != 0
    /* BigUint::ZERO */
    {
        let lsb = n & 1; // usize::from(n.bit(0));
        let mut tmp2 = q;
        for e in tmp2[1..].iter_mut().step_by(2) {
            *e = MODULO - *e;
        }
        conv(&tmp2, &p, &mut tmp);
        for i in 0..2000 {
            p[i] = tmp[lsb..][i << 1];
        }
        conv(&q, &tmp2, &mut tmp);
        for i in 0..=2000 {
            q[i] = tmp[i << 1];
        }
        /* dbg!(p.iter().filter(|&&v| v != 0).count());
        dbg!(q.iter().filter(|&&v| v != 0).count()); */
        n >>= 1;
    }
    let res = p[0];
    println!("res = {res}, took {:?}", start.elapsed());
}
