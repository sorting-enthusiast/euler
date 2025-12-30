use itertools::Itertools;

use crate::utils::{FIArray::FIArrayU64, primes::wheel_sieve};

const N: u64 = 1e15 as _;
const SQRT_N: u64 = N.isqrt();

const K: usize = 40;

pub fn main() {
    let primes = wheel_sieve(SQRT_N);
    let (splitting, inert_ramified): (Vec<u64>, Vec<u64>) =
        primes.iter().partition(|&&p| p % 5 == 1 || p % 5 == 4);
    dbg!(&splitting[..10]);
    let mut res = 0;
    let mut smooth = FIArrayU64::eps(N);
    let keys = FIArrayU64::keys(N).collect_vec().into_boxed_slice();
    dbg!(N / 10_253_399_761_u64);
    for &p in &inert_ramified {
        if p == 5 {
            continue;
        }
        if p > SQRT_N {
            break;
        }
        if p * p > N / 10_253_399_761_u64 {
            break;
        }
        let start = keys.partition_point(|&v| v < p * p);
        for (i, &v) in keys.iter().enumerate().skip(start) {
            if v > N / 10_253_399_761_u64 {
                break;
            }
            smooth.arr[i] += smooth[v / (p * p)];
        }
    }
    for (i, &v) in keys.iter().enumerate() {
        if v > N / 10_253_399_761_u64 {
            break;
        }
        smooth.arr[i] += smooth[v / 5];
    }
    println!("smooth counting done");
    let mut max = 0;
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
                res += smooth[val / p3];
                if val / p3 > max {
                    max = val / p3;
                }
            }
        }
    }
    println!("9 3 1");

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
                res += smooth[val / p3];
                if val / p3 > max {
                    max = val / p3;
                }
            }
        }
    }
    println!("7 4 1");

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
                    res += smooth[val / p4];
                    if val / p4 > max {
                        max = val / p4;
                    }
                }
            }
        }
    }
    println!("9 1 1 1");

    // 4 1 1 1 1
    for &p1 in &splitting {
        if p1.pow(4) > N {
            break;
        }
        for (i, &p2) in splitting.iter().enumerate() {
            if p1 == p2 {
                continue;
            }
            if p2 > N / p1.pow(4) {
                break;
            }
            for (j, &p3) in splitting[i + 1..].iter().enumerate() {
                if p1 == p3 {
                    continue;
                }
                if p3 > N / (p2 * p1.pow(4)) {
                    break;
                }
                for (k, &p4) in splitting[i + j + 2..].iter().enumerate() {
                    if p1 == p4 {
                        continue;
                    }
                    if p4 > N / (p3 * p2 * p1.pow(4)) {
                        break;
                    }
                    for &p5 in &splitting[i + j + k + 3..] {
                        if p1 == p5 {
                            continue;
                        }
                        let val = N / (p4 * p3 * p2 * p1.pow(4));
                        if p5 > val {
                            break;
                        }
                        res += smooth[val / p5];
                        if val / p5 > max {
                            max = val / p5;
                        }
                    }
                }
            }
        }
    }
    println!("4 1 1 1 1");

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
                    res += smooth[val / p4];
                    if val / p4 > max {
                        max = val / p4;
                    }
                }
            }
        }
    }
    println!("4 3 1 1");

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
                res += smooth[val / p3.pow(3)];
                if val / p3.pow(3) > max {
                    max = val / p3.pow(3);
                }
            }
        }
    }
    println!("4 3 3");

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
                    if p4 > val {
                        break;
                    }
                    res += smooth[val / p4];
                    if val / p4.pow(2) > max {
                        max = val / p4.pow(2);
                    }
                }
            }
        }
    }
    println!("2 2 2 2");

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
                res += smooth[val / p3.pow(2)];
                if val / p3.pow(2) > max {
                    max = val / p3.pow(2);
                }
            }
        }
    }
    println!("8 2 2");
    dbg!(max);
    dbg!(res);
}
