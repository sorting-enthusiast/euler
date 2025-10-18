use crate::utils::multiplicative_function_summation::{sum_n_i64, totient_sum};

// sum of convolution of N and totient
// dirichlet hyperbola method
pub fn main() {
    const N: i64 = 1e11 as _;
    const MOD: i64 = 998_244_353;
    let start = std::time::Instant::now();

    let sqrtn = N.isqrt();
    let totient_sums = totient_sum::<MOD>(N);
    let mut sum = (totient_sums[N] + sum_n_i64::<MOD>(N)) % MOD;
    for i in 2..=sqrtn {
        sum += (i * totient_sums[N / i]) % MOD;
        sum += ((totient_sums.arr[i as usize - 1] - totient_sums.arr[i as usize - 2]) % MOD
            * sum_n_i64::<MOD>(N / i))
            % MOD;
        sum %= MOD;
    }
    sum -= (totient_sums[sqrtn] * sum_n_i64::<MOD>(sqrtn)) % MOD;
    if sum < 0 {
        sum += MOD;
    }
    let end = start.elapsed();
    println!("res = {sum}, took {end:?}");
}
