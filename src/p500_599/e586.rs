use crate::utils::primes::wheel_sieve;

const N: u64 = 1e15 as _;
const SQRT_N: u64 = N.isqrt();
const MAX_KEY: u64 = N / (11u64.pow(4) * 19 * 29 * 31 * 41);
const MAX_PRIME: u64 = N / (11u64.pow(4) * 19 * 29 * 31);
const K: usize = 40;

pub fn main() {
    let start = std::time::Instant::now();
    let primes = wheel_sieve(MAX_PRIME);
    let (splitting, inert_ramified): (Vec<u64>, Vec<u64>) =
        primes.iter().partition(|&&p| p % 5 == 1 || p % 5 == 4);
    let mut res = 0;
    let mut smooth = vec![1u64; MAX_KEY as usize];
    for &p in &inert_ramified {
        if p == 5 {
            continue;
        }
        if p * p > MAX_KEY {
            break;
        }
        for v in p * p..=MAX_KEY {
            smooth[v as usize - 1] += smooth[(v / (p * p)) as usize - 1];
        }
    }
    for v in 5..=MAX_KEY {
        smooth[v as usize - 1] += smooth[(v / 5) as usize - 1];
    }

    // possible multiplicative partitions of 81:

    // 2 2 2 2
    for (i, &p1) in splitting.iter().enumerate() {
        if p1 * p1 > N {
            break;
        }
        for (j, &p2) in splitting[i + 1..].iter().enumerate() {
            if p2 * p2 > N / p1.pow(2) {
                break;
            }
            for (k, &p3) in splitting[i + j + 2..].iter().enumerate() {
                if p3.pow(2) > N / (p2 * p1).pow(2) {
                    break;
                }
                for &p4 in &splitting[i + j + k + 3..] {
                    let val = N / (p3 * p2 * p1).pow(2);
                    if p4.pow(2) > val {
                        break;
                    }
                    res += smooth[(val / p4.pow(2)) as usize - 1];
                }
            }
        }
    }

    // 8 2 2
    for &p1 in &splitting {
        if p1.pow(8) > N {
            break;
        }
        for (i, &p2) in splitting.iter().enumerate() {
            if p1 == p2 {
                continue;
            }
            if p2.pow(2) > N / p1.pow(8) {
                break;
            }
            for &p3 in &splitting[i + 1..] {
                if p1 == p3 {
                    continue;
                }
                let val = N / (p2.pow(2) * p1.pow(8));
                if p3.pow(2) > val {
                    break;
                }
                res += smooth[(val / p3.pow(2)) as usize - 1];
            }
        }
    }

    // possible multiplicative partitions of 80:

    // 9 3 1
    for &p1 in &splitting {
        if p1.pow(9) > N {
            break;
        }
        for &p2 in &splitting {
            if p1 == p2 {
                continue;
            }
            if p2.pow(3) > N / p1.pow(9) {
                break;
            }
            for &p3 in &splitting {
                if p1 == p3 || p2 == p3 {
                    continue;
                }
                let val = N / (p1.pow(9) * p2.pow(3));
                if p3 > val {
                    break;
                }
                res += smooth[(val / p3) as usize - 1];
            }
        }
    }

    // 7 4 1
    for &p1 in &splitting {
        if p1.pow(7) > N {
            break;
        }
        for &p2 in &splitting {
            if p1 == p2 {
                continue;
            }
            if p2.pow(4) > N / p1.pow(7) {
                break;
            }
            for &p3 in &splitting {
                if p1 == p3 || p2 == p3 {
                    continue;
                }
                let val = N / (p1.pow(7) * p2.pow(4));
                if p3 > val {
                    break;
                }
                res += smooth[(val / p3) as usize - 1];
            }
        }
    }

    // 9 1 1 1
    for &p1 in &splitting {
        if p1.pow(9) > N {
            break;
        }
        for (i, &p2) in splitting.iter().enumerate() {
            if p1 == p2 {
                continue;
            }
            if p2 > N / p1.pow(9) {
                break;
            }
            for (j, &p3) in splitting[i + 1..].iter().enumerate() {
                if p1 == p3 {
                    continue;
                }
                if p2 * p3 > N / p1.pow(9) {
                    break;
                }
                for &p4 in &splitting[i + j + 2..] {
                    if p1 == p4 {
                        continue;
                    }
                    let val = N / (p1.pow(9) * p2 * p3);
                    if p4 > val {
                        break;
                    }
                    res += smooth[(val / p4) as usize - 1];
                }
            }
        }
    }

    // 4 3 1 1
    for &p1 in &splitting {
        if p1.pow(4) > N {
            break;
        }
        for &p2 in &splitting {
            if p1 == p2 {
                continue;
            }
            if p2.pow(3) > N / p1.pow(4) {
                break;
            }
            for (i, &p3) in splitting.iter().enumerate() {
                if p1 == p3 || p2 == p3 {
                    continue;
                }
                if p3 > N / (p2.pow(3) * p1.pow(4)) {
                    break;
                }
                for &p4 in &splitting[i + 1..] {
                    if p1 == p4 || p2 == p4 {
                        continue;
                    }
                    let val = N / (p3 * p2.pow(3) * p1.pow(4));
                    if p4 > val {
                        break;
                    }
                    res += smooth[(val / p4) as usize - 1];
                }
            }
        }
    }

    // 4 3 3
    for &p1 in &splitting {
        if p1.pow(4) > N {
            break;
        }
        for (i, &p2) in splitting.iter().enumerate() {
            if p1 == p2 {
                continue;
            }
            if p2.pow(3) > N / p1.pow(4) {
                break;
            }
            for &p3 in &splitting[i + 1..] {
                if p1 == p3 {
                    continue;
                }
                let val = N / (p2.pow(3) * p1.pow(4));
                if p3.pow(3) > val {
                    break;
                }
                res += smooth[(val / p3.pow(3)) as usize - 1];
            }
        }
    }
    dbg!(start.elapsed());
    // 4 1 1 1 1
    // oddly slow, can probably be optimised
    for &p1 in &splitting {
        let lim = N / p1.pow(4);
        if lim == 0 {
            break;
        }
        for (i, &p2) in splitting.iter().enumerate() {
            if p1 == p2 {
                continue;
            }
            let lim = lim / p2;
            if lim == 0 {
                break;
            }
            for (j, &p3) in splitting[i + 1..].iter().enumerate() {
                if p1 == p3 {
                    continue;
                }
                let lim = lim / p3;
                if lim == 0 {
                    break;
                }
                for (k, &p4) in splitting[i + j + 2..].iter().enumerate() {
                    if p1 == p4 {
                        continue;
                    }
                    let lim = lim / p4;
                    if lim == 0 {
                        break;
                    }
                    for &p5 in &splitting[i + j + k + 3..] {
                        if p1 == p5 {
                            continue;
                        }
                        let lim = lim / p5;
                        if lim == 0 {
                            break;
                        }
                        res += smooth[lim as usize - 1];
                    }
                }
            }
        }
    }

    println!("res = {res}, took {:?}", start.elapsed());
}
