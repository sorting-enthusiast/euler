use crate::utils::FIArray::FIArray;
const N_2: usize = N * N;
const N: usize = 1e8 as _;
const MOD: usize = 1e9 as usize + 7;
// https://oeis.org/A046895, can optimize from O(N) to O(N^2/3) using convex-hull based summation
pub fn main() {
    let start = std::time::Instant::now();
    let unit = FIArray::unit(N_2);
    let id = FIArray::id::<MOD>(N_2);
    let len = id.arr.len();
    let mut sdn =
        id.arr[len - 1] + unit.arr[len - 1] % MOD + MOD - (unit.arr[N - 1] * id.arr[N - 1]) % MOD;
    sdn %= MOD;

    for i in 2..=N {
        sdn += ((unit.arr[i - 1] - unit.arr[i - 2]) * id.arr[len - i]) % MOD;
        sdn += (id.arr[i - 1] + MOD - id.arr[i - 2]) % MOD * unit.arr[len - i] % MOD;
        sdn %= MOD;
    }
    let mut sdn4 = id.arr[len - 4] + unit.arr[len - 4] % MOD + MOD
        - (unit.arr[N / 2 - 1] * id.arr[N / 2 - 1]) % MOD;
    sdn4 %= MOD;
    for i in 2..=N / 2 {
        sdn4 += ((unit.arr[i - 1] - unit.arr[i - 2]) * id[N_2 / (4 * i)]) % MOD;
        sdn4 += ((id.arr[i - 1] + MOD - id.arr[i - 2]) % MOD * unit[N_2 / (4 * i)] % MOD) % MOD;
        sdn4 %= MOD;
    }
    let mut res = 1 + 8 * sdn + 32 * MOD - 32 * sdn4;
    res %= MOD;
    println!("res = {res}, took {:?}", start.elapsed());
}
