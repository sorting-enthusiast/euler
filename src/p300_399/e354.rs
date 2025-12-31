use crate::utils::primes::wheel_sieve;

const N: u64 = 5e11 as _;
const SQRT_N: u64 = N.isqrt();
const MAX_KEY: u64 = N / (3 * 7u64.pow(2) * 13u64.pow(2) * 19);
const MAX_PRIME: u64 = N / (3 * 7u64.pow(2) * 13u64.pow(2));

pub fn main() {
    let start = std::time::Instant::now();
    let primes = wheel_sieve(MAX_PRIME);
    let (splitting, inert_ramified): (Vec<u64>, Vec<u64>) =
        primes.iter().partition(|&&p| p % 3 == 1);
    dbg!(&splitting[..10]);

    let mut res = 0;
    let mut smooth = vec![1u64; MAX_KEY as usize];
    for &p in &inert_ramified {
        if p * p > MAX_KEY {
            break;
        }
        for v in p..=MAX_KEY {
            smooth[v as usize - 1] += smooth[(v / p) as usize - 1];
        }
    }

    // 2 2 1
    for (i, &p1) in splitting.iter().enumerate() {
        let lim = N / (p1 * p1);
        if lim < 3 {
            break;
        }
        for &p2 in &splitting[i + 1..] {
            let lim = lim / (p2 * p2);
            if lim < 3 {
                break;
            }
            for &p3 in &splitting {
                if p1 == p3 || p2 == p3 {
                    continue;
                }
                let lim = lim / p3;
                if lim < 3 {
                    break;
                }
                res += smooth[(lim as usize) / 3 - 1];
            }
        }
    }

    // 7 2
    for &p1 in &splitting {
        let lim = N / p1.pow(7);
        if lim < 3 {
            break;
        }
        for &p2 in &splitting {
            if p1 == p2 {
                continue;
            }
            let lim = lim / (p2 * p2);
            if lim < 3 {
                break;
            }
            res += smooth[(lim as usize) / 3 - 1];
        }
    }
    // 12 1
    for &p1 in &splitting {
        let lim = N / p1.pow(12);
        if lim < 3 {
            break;
        }
        for &p2 in &splitting {
            if p1 == p2 {
                continue;
            }
            let lim = lim / p2;
            if lim < 3 {
                break;
            }
            res += smooth[(lim as usize) / 3 - 1];
        }
    }
    println!("res = {res}, took {:?}", start.elapsed());
}
