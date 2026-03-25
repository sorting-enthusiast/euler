use crate::utils::primes::primality::is_prime;

const N: usize = 1e12 as _;
pub fn main() {
    let mut primes = vec![];
    for i in 0..=N.ilog2() {
        let p2 = 2usize.pow(i);
        for j in 0..=(N / p2).ilog(3) {
            let p3 = 3usize.pow(j);
            for k in 0..=(N / p2 / p3).ilog(5) {
                let p5 = 5usize.pow(k);
                let p = (p2 * p3 * p5) + 1;
                if p <= N && is_prime(p as _) {
                    primes.push(p);
                }
            }
        }
    }
    primes.sort_unstable();
    dbg!(dfs(1, N, &primes[3..]));
}
fn dfs(acc: usize, lim: usize, primes: &[usize]) -> u32 {
    let mut res = acc as u32 * sum_5_smooth(lim);
    for (i, &p) in primes.iter().enumerate() {
        if p > lim {
            break;
        }
        res += dfs(acc * p, lim / p, &primes[i + 1..]);
    }
    res
}
fn sum_5_smooth(lim: usize) -> u32 {
    let mut res = 0;
    for i in 0..=lim.ilog2() {
        let p2 = 2usize.pow(i);
        for j in 0..=(lim / p2).ilog(3) {
            let p3 = 3usize.pow(j);
            for k in 0..=(lim / p2 / p3).ilog(5) {
                let p5 = 5usize.pow(k);
                res += (p2 * p3 * p5) as u32;
            }
        }
    }
    res
}
