use std::time::Instant;
// x^COEFFS = x + 1
fn polymulmod<const MODULO: u64, const COEFFS: usize>(
    poly1: &[u64; COEFFS],
    poly2: &[u64; COEFFS],
    res: &mut [u64; COEFFS],
) {
    res.fill(0);
    for (i, &c1) in poly1.iter().enumerate().filter(|&(_, &c)| c != 0) {
        for (j, &c2) in poly2[..COEFFS - i].iter().enumerate() {
            res[i + j] += c1 * c2;
        }
        for (j, &c2) in poly2[COEFFS - i..].iter().enumerate() {
            let c = c1 * c2;
            res[j] += c;
            res[j + 1] += c;
        }
    }
    res.iter_mut().for_each(|c| *c %= MODULO);
}
pub fn main() {
    const MODULO: u64 = 20_092_010;
    const LEN: usize = 2000;

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
            polymulmod::<MODULO, LEN>(r, x, tmp);
            std::mem::swap(&mut tmp, &mut r);
            //r = (r * x) % MODULO;
        }
        polymulmod::<MODULO, LEN>(x, x, tmp);
        std::mem::swap(&mut tmp, &mut x);
        //x = (x * x) % MODULO;
        exp >>= 1;
    }
    polymulmod::<MODULO, LEN>(r, x, tmp);

    let sum = (2 * (*tmp).into_iter().sum::<u64>() + tmp[LEN - 1]) % MODULO;
    println!("{:?}", start.elapsed());
    println!("{sum}");
}
