use crate::{
    utils::fenwick::FenwickTreeU32Mod, utils::multiplicative_function_summation::totient_sieve,
};
const MOD: u32 = 1e8 as _;
const N: usize = 2e7 as _;
pub fn main() {
    let start = std::time::Instant::now();
    let phi = totient_sieve(N + 1);
    let mut bit = FenwickTreeU32Mod::<MOD>::new(N + 1, 0);
    let mut f_n = 0;
    for n in (6..=N).rev() {
        let k = phi[n] as usize;
        f_n = 1 + bit.sum(n - 1) + MOD - bit.sum(k);
        f_n %= MOD;
        bit.add(k, f_n);
    }
    println!("res = {f_n}, took {:?}", start.elapsed());
}
