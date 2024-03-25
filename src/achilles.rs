use crate::sieve_of_pritchard::sift;
use std::time::Instant;
static mut V: Vec<(u32, u32)> = vec![];
const N: u64 = 1e18 as u64;
pub fn gcd(mut u: u64, mut v: u64) -> u64 {
    if u == 0 || v == 0 {
        return u | v;
    }
    let shift = (u | v).trailing_zeros();
    u >>= u.trailing_zeros();
    loop {
        v >>= v.trailing_zeros();
        v -= u;
        let m = (v as i64 >> 63) as u64;
        u += v & m;
        v = (v + m) ^ m;
        if v == 0 {
            break;
        }
    }
    u << shift
}
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
            let mut ret = second as u64;
            while x % first as u64 == 0 {
                x /= first as u64;
                ret += 1;
            }
            if ret < 2 {
                return false;
            }
            s = gcd(s, ret);
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
            s = gcd(s, multiplicity);
        }
    }
    s == 1 && x == 1
}
fn dfs(depth: u32, w: u64, ls: u32, phi: u64, la: u32, z: bool, primes: &[u64]) -> u64 {
    let mut ans = 0;
    if ls == 1 && la >= 3 && z {
        ans += check_achilles(phi, primes) as u64;
    }
    if depth as usize > primes.len() {
        return ans;
    }
    let s = primes[depth as usize - 1];
    let mut tmp = s * s;
    let mut p = 2;
    if w > N / tmp {
        return ans;
    }
    while w <= N / tmp {
        unsafe {
            V.push((s as u32, p - 1));
        }
        ans += dfs(
            depth + 1,
            w * tmp,
            gcd(ls as u64, p as u64) as u32,
            phi * (s - 1),
            p,
            true,
            primes,
        );
        unsafe {
            V.pop();
        }
        if s > N / tmp {
            break;
        }
        tmp *= s;
        p += 1;
    }
    ans + dfs(depth + 1, w, ls, phi, la, false, primes)
}
