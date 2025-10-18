use crate::utils::{
    FIArray::FIArrayI64,
    multiplicative_function_summation::{sum_n_i64, totient_sieve},
};
const N: i64 = 99_999_999_019;
const MOD: i64 = 999_999_017;
const fn sum_squares(x: i64) -> i64 {
    let x = x % MOD;
    (((x * (x + 1)) % MOD * (2 * x + 1)) % MOD * const { modinv(6) }) % MOD
}

fn sum_f(x: i64) -> FIArrayI64 {
    let y = if x > 1023 {
        (1e8 as usize).min((x as f64).powf(2. / 3.) as usize >> 1)
    } else {
        x as usize
    };
    let mut small = totient_sieve(y + 1);
    for i in 2..=y {
        small[i] *= i as i64 % MOD;
        small[i] %= MOD;
        small[i] += small[i - 1];
        small[i] %= MOD;
    }

    let mut ret = FIArrayI64::new(x);

    for (i, v) in FIArrayI64::keys(x).enumerate() {
        if v as usize <= y {
            ret.arr[i] = small[v as usize];
            continue;
        }
        let mut f_v = sum_squares(v) - sum_n_i64::<MOD>(v);
        let vsqrt = v.isqrt();
        for i in 2..=vsqrt {
            f_v -= (i * ret[v / i]) % MOD;
            f_v -= ((ret.arr[i as usize - 1] - ret.arr[i as usize - 2]) * sum_n_i64::<MOD>(v / i))
                % MOD;
            f_v %= MOD;
        }
        f_v += (ret[vsqrt] * sum_n_i64::<MOD>(vsqrt)) % MOD;
        f_v %= MOD;
        if f_v < 0 {
            f_v += MOD;
        }
        ret.arr[i] = f_v;
    }

    ret
}

const fn powmod(mut x: i64, mut exp: i64) -> i64 {
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
const fn modinv(x: i64) -> i64 {
    powmod(x, MOD - 2)
}
// average lcm is (1 + sum d|n d totient(d)) / 2, found proof on mathexchange:
// https://math.stackexchange.com/questions/761670/how-to-find-this-lcm-sum-function-textlcm1-n-textlcm2-n-cdots-t
// nested dirichlet hyperbola:
// first on u * (N phi)
// then on N phi: (N phi) = (N (mu * N)) = (N^2) * (N mu) =>
// N * (N phi) = N^2, use this to compute prefix sum of N phi in O(n^(2/3)) time
pub fn main() {
    let start = std::time::Instant::now();
    let prefix_sums = sum_f(N);
    let mut sum = N % MOD + prefix_sums[N];
    sum -= MOD;
    let sqrt = N.isqrt();
    for i in 2..=sqrt {
        sum += ((prefix_sums.arr[i as usize - 1] - prefix_sums.arr[i as usize - 2]) * (N / i)
            % MOD)
            % MOD;
        sum += prefix_sums[N / i];
        sum %= MOD;
    }
    sum -= (sqrt * prefix_sums[sqrt]) % MOD;
    sum %= MOD;
    if sum < 0 {
        sum += MOD;
    }
    sum += N % MOD;
    sum %= MOD;
    sum *= const { modinv(2) };
    sum %= MOD;
    let end = start.elapsed();
    println!("sum = {sum}, took {end:?}");
}
