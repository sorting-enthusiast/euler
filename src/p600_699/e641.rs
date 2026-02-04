use crate::utils::primes::wheel_sieve;
#[must_use]
pub const fn iroot<const k: u128>(x: u128) -> u128 {
    let mut rt: u128 = 1 << (1 + x.ilog2().div_ceil(k as _));
    let mut x_div_rtk1 = x / rt.pow(k as u32 - 1);
    while rt > x_div_rtk1 {
        rt = (rt * (k - 1) + x_div_rtk1) / k;
        x_div_rtk1 = x / rt.pow(k as u32 - 1);
    }
    rt
}
// can be optimized a ton, but is fast enough for me
pub fn main() {
    // only 1 and 5 are invertible mod 6,
    // n'th die is 1 iff it is a product of an even number of primes with multiplicity = 4 mod 6, and any number of primes with multiplicity = 0 mod 6
    const N: u128 = 1_000_000_000u128.pow(4);
    const LIM: u64 = (iroot::<4>(N) / 2) as u64;
    fn dfs(lim: u128, parity: bool, primes: &[u64]) -> usize {
        let mut res = parity.into();
        let lim_1 = iroot::<6>(lim);
        let lim_5 = iroot::<4>(lim);

        for (i, &p) in primes.iter().enumerate() {
            let p = u128::from(p);
            if p > lim_5 {
                break;
            }
            let p6 = p.pow(6);
            let mut new_lim = lim / p.pow(4);
            loop {
                res += dfs(new_lim, !parity, &primes[i + 1..]);
                if new_lim < p6 {
                    break;
                }
                new_lim /= p6;
            }
            if p > lim_1 {
                continue;
            }
            new_lim = lim / p6;
            loop {
                res += dfs(new_lim, parity, &primes[i + 1..]);
                if new_lim < p6 {
                    break;
                }
                new_lim /= p6;
            }
        }
        res
    }
    let start = std::time::Instant::now();
    let primes = wheel_sieve(LIM);
    let res = dfs(N, true, &primes);
    println!("res = {res}, took {:?}", start.elapsed());
}
