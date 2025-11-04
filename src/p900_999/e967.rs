use std::collections::HashSet;

use itertools::Itertools;

const N: usize = 1e18 as _;
const B: usize = 120;

const PRIMES: [usize; 30] = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
    101, 103, 107, 109, 113,
];

pub fn main() {
    fn get_keys(keys: &mut HashSet<usize>, lim: usize, primes: &[usize]) {
        keys.insert(lim);
        for (i, &p) in primes.iter().enumerate() {
            if p > lim {
                break;
            }
            get_keys(keys, lim / p, &primes[i + 1..]);
        }
    }
    let start = std::time::Instant::now();
    let mut keys = HashSet::new();
    get_keys(&mut keys, N, &PRIMES);
    println!("Collected {} keys, took {:?}", keys.len(), start.elapsed());

    let mut keys = keys.into_iter().collect_vec();
    keys.sort_unstable();
    let keys = keys;

    let get_index = |v| keys.partition_point(|&e| e < v);

    let mut c = [keys.clone(), vec![0; keys.len()], vec![0; keys.len()]];
    for p in PRIMES {
        if p > B {
            break;
        }
        for (i, &v) in keys.iter().enumerate().rev() {
            if v < p {
                break;
            }
            let j = get_index(v / p);
            for k in 0..3 {
                c[k][i] -= c[k][j];
                c[k][i] += c[(k + 3 - (p % 3)) % 3][j];
            }
        }
    }
    let res = c[0][get_index(N)];
    println!("res = {res}, took {:?}", start.elapsed());
}
