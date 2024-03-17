use crate::gcd;
use crate::sieve_of_pritchard::sift;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::time::Instant;

pub fn main() {
    //302:
    //achilles number: every prime factor has multiplicity > 1 and gcd(multiplicities)=1
    //strong achilles number: both num and totient(num) are achilles numbers
    //note: the multiplicity of the largest prime factor has to be at least 3 for num to be a strong achilles number
    const N: u64 = 1e18 as u64;

    let achilles = gen_achilles(N);
    println!("found achilles nums");
    //dbg!(achilles.len());

    let start = Instant::now();
    let count2 = achilles.iter().fold(0u64, |acc, (_, tot_x)| {
        acc + achilles.contains_key(tot_x) as u64
    });
    let end = start.elapsed();
    println!("{:?}", end);

    dbg!(count2);

    let start = Instant::now();
    let count1 = achilles
        .par_iter()
        .fold_with(0u64, |acc, (_, tot_x)| {
            acc + achilles.contains_key(tot_x) as u64
        })
        .sum::<u64>();
    let end = start.elapsed();
    println!("{:?}", end);

    dbg!(count1);
}

fn gen_achilles(limit: u64) -> HashMap<u64, u64> {
    let primes = sift((limit as f64).sqrt().sqrt() as u64);
    println!("found primes");
    let mut totients = vec![1; ((limit / 8) as f64).sqrt() as usize].into_boxed_slice();
    for i in 2..((limit / 8) as f64).sqrt() as u64 {
        for &prime in &primes {
            if prime * prime > i {
                totients[i as usize - 1] = i - 1;
                break;
            }
            if i % prime == 0 {
                /* totient = if num / prime % prime == 0 {
                    totients[(num / prime) as usize - 1] * prime
                } else {
                    totients[prime as usize - 1] * totients[(num / prime) as usize - 1]
                }; */
                let d = if i / prime % prime == 0 { prime } else { 1 };
                totients[i as usize - 1] =
                    totients[prime as usize - 1] * totients[(i / prime) as usize - 1] * d
                        / totients[d as usize - 1];
                break;
            }
        }
    }
    drop(primes);
    println!("precomputed totients");
    //let pps = perfect_powers(limit);
    //println!("found perfect powers");
    let mut achilles = HashMap::new();
    for b in 2..((limit / 4) as f64).cbrt() as u128 {
        let b_cubed = b * b * b;
        let tot_b = b * b * totients[b as usize - 1] as u128;
        for a in 2..((limit / 8) as f64).sqrt() as u128 {
            if a == b {
                continue;
            }
            let p = b_cubed * a * a;
            if p > limit as u128 {
                break;
            }
            //if !pps.contains(&(p as u64)) {
            let tot_a = a * totients[a as usize - 1] as u128;
            let d = gcd(b as u64, a as u64);
            let tot_p = tot_a * tot_b * d as u128 / totients[d as usize - 1] as u128;
            achilles.insert(p as u64, tot_p as u64);
            //}
        }
    }
    for i in (2..(limit as f64).sqrt() as u128).rev() {
        if i % 1_000_000 == 0 {
            println!("{i}");
        }
        let mut p = i * i;
        while p <= limit as u128 {
            achilles.remove(&(p as u64));
            p *= i;
        }
    }
    achilles
}
fn perfect_powers(limit: u64) -> HashSet<u64> {
    let mut pps = HashSet::new();
    for i in (2..(limit as f64).sqrt() as u128).rev() {
        if i % 1_000_000 == 0 {
            println!("{i}");
        }
        let mut p = i * i;
        while p <= limit as u128 {
            pps.insert(p as u64);
            p *= i;
        }
    }
    pps
}
