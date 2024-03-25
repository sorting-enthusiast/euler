use crate::{math, sieve_of_pritchard::sift};
use std::time::Instant;
static mut V: Vec<(u64, u64)> = vec![];
const N: u64 = 1e18 as u64;

pub fn main() {
    //302:
    //achilles number: every prime factor has multiplicity > 1 and gcd(multiplicities)=1
    //strong achilles number: both num and totient(num) are achilles numbers
    //note: the multiplicity of the largest prime factor has to be at least 3 for num to be a strong achilles number
    let start = Instant::now();
    let primes = sift((N as f64 / 4.0).cbrt() as u64);
    let ans = dfs(1, 1, 0, 1, 0, false, &primes);
    let end = start.elapsed();
    dbg!(end);
    dbg!(ans);
}
fn check_achilles(k: u64, primes: &[u64]) -> bool {
    let mut s = 0;
    let mut x = k;
    unsafe {
        for &(first, second) in V.iter() {
            let mut ret = second;
            while x % first == 0 {
                x /= first;
                ret += 1;
            }
            if ret < 2 {
                return false;
            }
            s = math::gcd(s, ret);
        }
        for &prime in primes {
            if prime * prime > x {
                break;
            }
            let mut multiplicity = 0;
            while x % prime == 0 {
                multiplicity += 1;
                x /= prime;
            }
            if multiplicity == 1 {
                return false;
            }
            s = math::gcd(s, multiplicity);
        }
    }
    s == 1 && x == 1
}
fn dfs(depth: u64, acc: u64, ls: u64, phi: u64, la: u64, z: bool, primes: &[u64]) -> u64 {
    let mut ans = 0;
    if ls == 1 && la > 2 && z {
        ans += check_achilles(phi, primes) as u64;
    }
    if depth as usize > primes.len() {
        return ans;
    }
    let prime = primes[depth as usize - 1];
    let mut prime_power = prime * prime;
    let mut power = 2;
    if acc > N / prime_power {
        return ans;
    }
    while acc <= N / prime_power {
        unsafe {
            V.push((prime, power - 1));
        }
        ans += dfs(
            depth + 1,
            acc * prime_power,
            math::gcd(ls, power),
            phi * (prime - 1),
            power,
            true,
            primes,
        );
        unsafe {
            V.pop();
        }
        if prime > N / prime_power {
            break;
        }
        prime_power *= prime;
        power += 1;
    }
    ans + dfs(depth + 1, acc, ls, phi, la, false, primes)
}
