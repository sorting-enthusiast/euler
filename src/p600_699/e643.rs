use crate::utils::multiplicative_function_summation::totient_sum;
const N: i64 = 1e11 as _;
const MOD: i64 = 1e9 as i64 + 7;

pub fn main() {
    let tot = totient_sum::<MOD>(N);
    let mut res = 0;
    let mut n = N >> 1;
    while n > 0 {
        res += tot[n] - 1;
        if res >= MOD {
            res -= MOD;
        }
        n >>= 1;
    }
    dbg!(res);
}
