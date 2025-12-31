use crate::utils::primes::wheel_sieve;

const N: u64 = 5e11 as _;
const SQRT_N: u64 = N.isqrt();
const BOUND_INT: u64 = N / 3;
const BOUND_SQRT: u64 = ((N as u128).pow(2) / 3).isqrt() as u64;
const MAX_BOUND: u64 = BOUND_SQRT;
const MAX_KEY: u64 = MAX_BOUND / (7u64.pow(2) * 13u64.pow(2) * 19);
const MAX_PRIME: u64 = MAX_BOUND / (7u64.pow(2) * 13u64.pow(2));

fn dfs(lim: u64, primes: &[u64]) -> u64 {
    let mut sum = 1;
    for (i, &p) in primes.iter().enumerate() {
        if p > lim {
            break;
        }
        let mut new_lim = lim;
        loop {
            new_lim /= p;
            sum += dfs(new_lim, &primes[i + 1..]);
            if p > new_lim {
                break;
            }
        }
    }
    sum
}

pub fn main() {
    let start = std::time::Instant::now();
    let primes = wheel_sieve(MAX_PRIME);
    let (splitting, inert_ramified): (Vec<u64>, Vec<u64>) =
        primes.iter().partition(|&&p| p % 3 == 1);
    let mut res = 0;
    let mut cache = vec![None; MAX_KEY as usize];
    // 2 2 1
    for (i, &p1) in splitting.iter().enumerate() {
        let lim = BOUND_INT / (p1 * p1);
        if lim == 0 {
            break;
        }
        for &p2 in &splitting[i + 1..] {
            let lim = lim / (p2 * p2);
            if lim == 0 {
                break;
            }
            for &p3 in &splitting {
                if p1 == p3 || p2 == p3 {
                    continue;
                }
                let lim = lim / p3;
                if lim == 0 {
                    break;
                }
                res += *cache[lim as usize - 1].get_or_insert_with(|| dfs(lim, &inert_ramified)); //smooth[lim as usize - 1];
            }
        }
    }

    // 7 2
    for &p1 in &splitting {
        let lim = BOUND_INT / p1.pow(7);
        if lim == 0 {
            break;
        }
        for &p2 in &splitting {
            if p1 == p2 {
                continue;
            }
            let lim = lim / (p2 * p2);
            if lim == 0 {
                break;
            }
            res += *cache[lim as usize - 1].get_or_insert_with(|| dfs(lim, &inert_ramified)); //smooth[lim as usize - 1];
        }
    }

    // 12 1
    for &p1 in &splitting {
        let lim = BOUND_INT / p1.pow(12);
        if lim == 0 {
            break;
        }
        for &p2 in &splitting {
            if p1 == p2 {
                continue;
            }
            let lim = lim / p2;
            if lim == 0 {
                break;
            }
            res += *cache[lim as usize - 1].get_or_insert_with(|| dfs(lim, &inert_ramified)); //smooth[lim as usize - 1];
        }
    }

    // 2 2 1
    for (i, &p1) in splitting.iter().enumerate() {
        let lim = BOUND_SQRT / (p1 * p1);
        if lim == 0 {
            break;
        }
        for &p2 in &splitting[i + 1..] {
            let lim = lim / (p2 * p2);
            if lim == 0 {
                break;
            }
            for &p3 in &splitting {
                if p1 == p3 || p2 == p3 {
                    continue;
                }
                let lim = lim / p3;
                if lim == 0 {
                    break;
                }
                res += *cache[lim as usize - 1].get_or_insert_with(|| dfs(lim, &inert_ramified)); //smooth[lim as usize - 1];
            }
        }
    }

    // 7 2
    for &p1 in &splitting {
        let lim = BOUND_SQRT / p1.pow(7);
        if lim == 0 {
            break;
        }
        for &p2 in &splitting {
            if p1 == p2 {
                continue;
            }
            let lim = lim / (p2 * p2);
            if lim == 0 {
                break;
            }
            res += *cache[lim as usize - 1].get_or_insert_with(|| dfs(lim, &inert_ramified)); //smooth[lim as usize - 1];
        }
    }

    // 12 1
    for &p1 in &splitting {
        let lim = BOUND_SQRT / p1.pow(12);
        if lim == 0 {
            break;
        }
        for &p2 in &splitting {
            if p1 == p2 {
                continue;
            }
            let lim = lim / p2;
            if lim == 0 {
                break;
            }
            res += *cache[lim as usize - 1].get_or_insert_with(|| dfs(lim, &inert_ramified)); //smooth[lim as usize - 1];
        }
    }

    println!("res = {res}, took {:?}", start.elapsed());
    dbg!(cache.iter().flatten().count());
}
