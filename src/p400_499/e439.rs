use crate::utils::{FIArray::FIArrayI64, multiplicative_function_summation::dirichlet_mulmod_i64};

const MOD: i64 = 1e9 as _;
const N: i64 = 1e11 as _;

// D(n)^2 = sum over 1<=d<=n d * S(n/d) = S(n) + sum over 1<d<=n d * S(n/d) =>
// S(n) = D(n)^2 - sum over 1<d<=n d * S(n/d)
pub fn main() {
    let start = std::time::Instant::now();
    let id = FIArrayI64::id::<MOD>(N);
    let mut s = FIArrayI64::unit(N);
    let d = dirichlet_mulmod_i64::<MOD>(&id, &s, N as _);
    s.arr.fill(0);

    for (j, v) in FIArrayI64::keys(N).enumerate() {
        let mut sum = d[v];
        sum *= sum;
        sum %= MOD;
        let mut l = 2;
        while l <= v {
            let div = v / l;
            let r = v / div;
            sum += MOD - (s[div] * (id[r] + MOD - id[l - 1]) % MOD) % MOD;
            if sum >= MOD {
                sum -= MOD;
            }
            l = r + 1;
        }
        s.arr[j] = sum;
    }
    let end = start.elapsed();
    println!("{}, {end:?}", s[N]);
}
