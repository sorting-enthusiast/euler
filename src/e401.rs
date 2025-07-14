use std::time::Instant;

const MOD: i128 = 1e9 as i128;
const N: i128 = 1e15 as i128;

const fn sum_squares(x: i128) -> i128 {
    ((x * (x + 1) * (2 * x + 1)) / 6) % MOD
}

// dirichlet hyperbola method, probably overkill
pub fn main() {
    let start = Instant::now();

    let mut sum = 0i128;
    let sqrtn = N.isqrt();

    for i in 1..=sqrtn {
        let ni = (N / i) % MOD; // prevent overflow of i128
        sum += sum_squares(ni);
        if sum >= MOD {
            sum -= MOD;
        }
        sum += (i * i * ni) % MOD;
        if sum >= MOD {
            sum -= MOD;
        }
    }
    sum -= (sqrtn * sum_squares(sqrtn)) % MOD;
    if sum < 0 {
        sum += MOD;
    }
    let end = start.elapsed();
    println!("sum = {sum}, took {end:?}");
}
