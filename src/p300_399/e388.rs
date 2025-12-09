use crate::utils::multiplicative_function_summation::mertens;
const fn sum_squares(x: i128) -> i128 {
    (x * (x + 1) * (2 * x + 1)) / 6
}
const fn sum_n(x: i128) -> i128 {
    if x & 1 == 0 {
        (x / 2) * (x + 1)
    } else {
        ((x + 1) / 2) * x
    }
}
const N: i128 = 1e12 as _;
const SQRT_N: i128 = N.isqrt();

// count coprime pairs * 3 for x,y,z planes, using the fact that the count is 2*totient sum - 1
// add contribution of all coprime triplets a,b,c, using generalized mobius inversion and dirichlet hyperbola method on
// n^3 = sum_{d<=n} C(n/d)
pub fn main() {
    let start = std::time::Instant::now();
    let m = mertens(N as _);

    let len = m.arr.len();

    let mut r = i128::from(m[N as i64]) + N.pow(3) - i128::from(m[SQRT_N as i64]) * SQRT_N.pow(3);
    let mut tot = i128::from(m[N as i64]) + sum_n(N) - sum_n(SQRT_N) * i128::from(m[SQRT_N as i64]);
    for i in 2..=SQRT_N as usize {
        r += i128::from(m.arr[i - 1] - m.arr[i - 2]) * (N / i as i128).pow(3);
        r += (3 * i * i - 3 * i + 1) as i128 * i128::from(m.arr[len - i]);

        tot += i as i128 * i128::from(m.arr[len - i]);
        tot += i128::from(m.arr[i - 1] - m.arr[i - 2]) * sum_n(N / i as i128);
    }
    let res = r + 3 + 3 * (2 * tot - 1);
    println!("res = {res}, took {:?}", start.elapsed());
    alt();
}
// alternate solution
fn alt() {
    let start = std::time::Instant::now();
    let m = mertens(N as _);

    let len = m.arr.len();

    let mut j2 = i128::from(m[N as i64]) + sum_squares(N)
        - sum_squares(SQRT_N) * i128::from(m[SQRT_N as i64]);
    let mut j1 = i128::from(m[N as i64]) + sum_n(N) - sum_n(SQRT_N) * i128::from(m[SQRT_N as i64]);
    for i in 2..=SQRT_N as usize {
        let mu_i = i128::from(m.arr[i - 1] - m.arr[i - 2]);

        j2 += (i * i) as i128 * i128::from(m.arr[len - i]);
        j2 += mu_i * sum_squares(N / i as i128);

        j1 += i as i128 * i128::from(m.arr[len - i]);
        j1 += mu_i * sum_n(N / i as i128);
    }
    let res = 1 + 3 * (j2 + j1);
    println!("res = {res}, took {:?}", start.elapsed());
}
