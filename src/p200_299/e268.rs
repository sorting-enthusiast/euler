use std::collections::HashSet;

use itertools::Itertools;

const N: usize = 1e16 as _;

const PRIMES: [usize; 25] = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
];

// same idea as 967
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
    let keys = keys.into_boxed_slice();

    let get_index = |v| keys.partition_point(|&e| e < v);

    let mut c = [
        keys.clone(),
        vec![0; keys.len()].into_boxed_slice(),
        vec![0; keys.len()].into_boxed_slice(),
        vec![0; keys.len()].into_boxed_slice(),
    ];
    for p in PRIMES {
        for (i, &v) in keys.iter().enumerate().rev() {
            if v < p {
                break;
            }
            let j = get_index(v / p);
            c[0][i] -= c[0][j];

            c[1][i] -= c[1][j];
            c[1][i] += c[0][j];

            c[2][i] -= c[2][j];
            c[2][i] += c[1][j];

            c[3][i] -= c[3][j];
            c[3][i] += c[2][j];
        }
    }
    let res = N
        - c[0][keys.len() - 1]
        - c[1][keys.len() - 1]
        - c[2][keys.len() - 1]
        - c[3][keys.len() - 1];
    println!("res = {res}, took {:?}", start.elapsed());
}
