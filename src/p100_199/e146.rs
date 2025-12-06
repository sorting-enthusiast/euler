use crate::utils::primes::primality::is_prime;

const N: u64 = 150e6 as _;
pub fn main() {
    let start = std::time::Instant::now();
    let mut sum = 0;
    'outer: for k in (0..).step_by(210) {
        for r in [10, 80, 130, 200] {
            let n = k + r;
            if n > N {
                break 'outer;
            }
            let nn = n * n;
            if !(is_prime(nn + 13)
                && is_prime(nn + 1)
                && is_prime(nn + 3)
                && is_prime(nn + 7)
                && is_prime(nn + 9)
                && is_prime(nn + 27))
            {
                continue;
            }
            let primes = (nn + 1..=nn + 27).filter(|&n| is_prime(n)).count();
            if primes == 6 {
                sum += n;
            }
        }
    }
    println!("res = {sum}, took {:?}", start.elapsed());
}
