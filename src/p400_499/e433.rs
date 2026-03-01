//use itertools::Itertools;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

pub fn main() {
    const N: u32 = 5e6 as _;
    const fn gcd(a: u32, b: u32) -> usize {
        let (mut x, mut y) = (b, a % b);
        let mut cnt = 1;
        while y != 0 {
            (x, y) = (y, x % y);
            cnt += 1;
        }
        cnt
    }
    /* let mut acc = 0;
    for a in 1..=N {
        for b in 1..=N {
            acc += gcd(a, b);
        }
    } */
    let ret = (1..=N)
        .into_par_iter()
        .fold(
            || 0,
            |acc, a| {
                acc + (1..=N)
                    .into_par_iter()
                    .fold(|| 0, |acc, b| acc + gcd(a, b))
                    .sum::<usize>()
            },
        )
        .sum::<usize>();
    dbg!(ret);
}
