use crate::sieve_of_pritchard::sift;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::time::Instant;

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
    const N: u64 = 1e18 as u64;

    let start = Instant::now();
    let achilles = gen_achilles(N);
    let end = start.elapsed();
    println!("{:?}", end);

    println!("found achilles nums");

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
//consider precomputing prime decomposition of all numbers in 2..((limit / 8) as f64).sqrt(), and using that/gcd to tell if num is a perfect power
fn gen_achilles(limit: u64) -> HashMap<u64, u64> {
    let primes = sift(((limit / 8) as f64).sqrt().sqrt() as u64);
    println!("found primes");
    let mut totients = vec![1; ((limit / 8) as f64).sqrt() as usize].into_boxed_slice();
    /* let mut factorization =
    vec![HashMap::new(); ((limit / 8) as f64).sqrt() as usize].into_boxed_slice(); */
    for i in 2..=((limit / 8) as f64).sqrt() as u64 {
        for &prime in &primes {
            if prime * prime > i {
                totients[i as usize - 1] = i as u32 - 1;
                break;
            }
            if i % prime == 0 {
                /* unsafe {
                    (*factorization.as_mut_ptr().add(i as usize - 1))
                        .clone_from(&factorization[(i / prime) as usize - 1]);
                }
                //factorization[i as usize - 1] = factorization[(i / prime) as usize - 1].clone();
                let multiplicity = *factorization[i as usize - 1]
                    .entry(prime as u32)
                    .and_modify(|x| *x += 1)
                    .or_insert(1);
                let d = if multiplicity > 1 { prime } else { 1 }; */

                let d = if i / prime % prime == 0 { prime } else { 1 };
                totients[i as usize - 1] = (totients[prime as usize - 1] as u64
                    * totients[(i / prime) as usize - 1] as u64
                    * d
                    / totients[d as usize - 1] as u64)
                    as u32;
                break;
            }
        }
        /* if factorization[i as usize - 1].is_empty() {
            factorization[i as usize - 1].insert(i as u32, 1u8);
            totients[i as usize - 1] = i as u32 - 1;
        } */
    }
    drop(primes);
    println!("precomputed totients");
    //let pps = perfect_powers(limit);
    //println!("found perfect powers");
    let mut achilles = HashMap::new();
    for b in 2..((limit / 4) as f64).cbrt() as u128 {
        let b_cubed = b * b * b;
        let tot_b = b * b * totients[b as usize - 1] as u128;
        for a in 2..=((limit / b_cubed as u64) as f64).sqrt() as u128 {
            if a == b {
                continue;
            }
            let p = b_cubed * a * a;
            /* let mut iter = factorization[a as usize - 1].iter();
            let mut g = {
                let (prime, &multiplicity) = iter.next().unwrap();
                if let Some(&multiplicity1) = factorization[b as usize - 1].get(prime) {
                    2 * multiplicity as u64 + 3 * multiplicity1 as u64
                } else {
                    2 * multiplicity as u64
                }
            };
            for (prime, &multiplicity) in iter {
                if g == 1 {
                    break;
                }
                if let Some(&multiplicity2) = factorization[b as usize - 1].get(prime) {
                    g = gcd(g, 2 * multiplicity as u64 + 3 * multiplicity2 as u64);
                } else {
                    g = gcd(g, 2 * multiplicity as u64);
                }
            }
            for (_, &multiplicity) in factorization[b as usize - 1]
                .iter()
                .filter(|(prime, _)| !factorization[a as usize - 1].contains_key(prime))
            {
                if g == 1 {
                    break;
                }
                g = gcd(g, 3 * multiplicity as u64);
            }
            if g != 1 {
                continue;
            } */
            let tot_a = a * totients[a as usize - 1] as u128;
            let d = gcd(b as u64, a as u64);
            let tot_p = tot_a * tot_b * d as u128 / totients[d as usize - 1] as u128;
            achilles.insert(p as u64, tot_p as u64);
        }
    }
    let mut count = 0;
    for i in 2..(limit as f64).sqrt() as u128 {
        let mut p = i * i;
        while p <= limit as u128 {
            count += achilles.remove(&(p as u64)).is_some() as u64;
            p *= i;
        }
    }
    dbg!(count);
    dbg!(achilles.len());
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
