use crate::utils::binary_search::binsearch;
use std::{ops::RangeInclusive, time::Instant};

const GUESS: u64 = 7_000_000;
const L: u64 = 10_000_000_000_000;
const R: u64 = 100_000_000_000_000;

const fn rounded_sqrt_iters(n: u64) -> u8 {
    let mut iter_count = 1;

    let mut x = GUESS;
    let mut x_new = n.div_ceil(x).midpoint(x);

    while x != x_new {
        x = x_new;
        x_new = n.div_ceil(x).midpoint(x);
        iter_count += 1;
    }
    iter_count
}

fn iteration(r: RangeInclusive<u64>) -> usize {
    let mut acc = 0;
    let &(mut start) = r.start();
    let &end = r.end();
    let mut n = r.count();
    let mut initial_val = rounded_sqrt_iters(start);
    let final_val = rounded_sqrt_iters(end);
    if initial_val == final_val {
        return initial_val as usize * n;
    }
    let mut f = binsearch(
        n,
        |i| rounded_sqrt_iters(i as u64 + start),
        |iter| iter == initial_val,
    );
    acc += f * initial_val as usize;
    loop {
        n -= f;
        start += f as u64;
        initial_val = rounded_sqrt_iters(start);
        if initial_val == final_val {
            return acc + n * initial_val as usize;
        }
        f = binsearch(
            n,
            |i| rounded_sqrt_iters(i as u64 + start),
            |iter| iter == initial_val,
        );
        acc += f * initial_val as usize;
        if n == f {
            break;
        }
    }
    acc
}

fn solve() {
    let start = Instant::now();
    let acc = (3_162_279..10_000_000)
        .map(|root| {
            ((f64::from(root) - 0.5).powi(2).ceil() as u64)
                ..=((f64::from(root) + 0.5).powi(2).floor() as u64)
        })
        .map(iteration)
        .sum::<usize>()
        + iteration(L..=10_000_005_311_562)
        + iteration(99_999_990_000_001..=(R - 1));
    let res = acc as f64 / (R - L) as f64;
    let end = start.elapsed();
    println!("res = {res:.10}, took {end:?}");
}
// adapted from forum post
fn alt() {
    let start = Instant::now();
    let mut acc = 0;
    let mut stack = vec![(L, R - 1, GUESS, 0)];
    while let Some((a, b, x, n)) = stack.pop() {
        let (qmin, qmax) = (a.div_ceil(x), b.div_ceil(x));
        for q in qmin..=qmax {
            let (na, nb) = (a.max((q - 1) * x + 1), b.min(q * x));
            let x_next = (x + q) >> 1;
            if x_next == x {
                acc += (nb - na + 1) * (n + 1);
            } else {
                stack.push((na, nb, x_next, n + 1));
            }
        }
    }
    let res = acc as f64 / (R - L) as f64;
    let end = start.elapsed();
    println!("res = {res:.10}, took {end:?}");
}
pub fn main() {
    solve();
    alt();
}
