#![warn(clippy::pedantic)]
#![allow(long_running_const_eval)]
#![allow(dead_code)]
#![allow(clippy::too_many_lines)]
#![allow(clippy::cast_possible_truncation)]
#![allow(clippy::inline_always)]
#![allow(clippy::cast_sign_loss)]
#![allow(clippy::cast_precision_loss)]
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]
#![allow(clippy::large_stack_arrays)]
use chrono::Local;

use crate::utils::{
    FIArray::{DirichletFenwickU128, FIArray, FIArrayU128},
    fenwick::FenwickTreeU128,
    multiplicative_function_summation::{
        inverse_pseudo_euler_transform, mertens, mertens_slow, pseudo_euler_transform,
    },
    primes::primecount::{lucy_fenwick, mertens_min25},
};

pub mod p0_99;
pub mod p100_199;
pub mod p200_299;
pub mod p300_399;
pub mod p400_499;
pub mod p500_599;
pub mod p600_699;
pub mod p700_799;
pub mod p800_899;
pub mod p900_999;
pub mod utils;

// digital root of n is just n mod 9 if n mod 9 != 0, otherwise 9
const fn is_target_little_endian() -> bool {
    u16::from_ne_bytes([1, 0]) == 1
}
// TODO: understand convex hull based lattice point counting, optimize dirichlet mul
pub fn main() {
    const { assert!(is_target_little_endian()) }; // some code relies on this
    println!("Started running at: {} ", Local::now().time());
    //p500_599::e580::main();
    /* assert_eq!(
        inverse_pseudo_euler_transform(pseudo_euler_transform(&lucy_fenwick(1e10 as _))),
        lucy_fenwick(1e10 as _)
    ); */
    let pi = lucy_fenwick(1000 << (2 * 27));
    for i in 23..=27 {
        print!("{i}:{},", pi[1000 << (2 * i)]);
    }
    println!();
    /* let sums = lucy_fenwick_simple(1000u128 << 47); //inverse_pseudo_euler_transform_u128(FIArrayU128::id::<0>(1000u128 << 43));
    dbg!(sums.x);
    for i in 44..=47 {
        print!("{i}:{}, ", sums[1000u128 << i]);
    }
    println!(); */

    p800_899::e890::main();
    p300_399::e362::main();
    utils::primes::primecount::main();
    const N: i64 = 1e12 as _;
    let start = std::time::Instant::now();
    let s1 = mertens(N);
    let end = start.elapsed();
    dbg!(end, s1[N]);
    let start = std::time::Instant::now();
    let s1 = mertens_min25(N);
    let end = start.elapsed();
    dbg!(end, s1[N]);
    //utils::primes::prime_sieves::main();
    println!("Finished running at: {} ", Local::now().time());
}
const fn iroot<const k: u128>(x: u128) -> u128 {
    let mut rt = 1u128 << (1 + x.ilog2().div_ceil(k as _));
    let mut x_div_rtk1 = x / rt.pow(k as u32 - 1);
    while rt > x_div_rtk1 {
        rt = (rt * (k - 1) + x_div_rtk1) / k;
        x_div_rtk1 = x / rt.pow(k as u32 - 1);
    }
    rt
}
fn inverse_pseudo_euler_transform_u128(a: FIArrayU128) -> FIArrayU128 {
    let n = a.x;
    let len = a.arr.len();
    let mut r = FIArrayU128::new(n);
    let mut a_bit = DirichletFenwickU128::from(a);
    let x = iroot::<3>(n) + 1;
    for i in 2..x {
        let cur = a_bit.bit.sum(i as usize - 1) - 1;
        if cur == 0 {
            continue;
        }
        r.arr[i as usize - 1] = cur;
        a_bit.sparse_mul_at_most_one(i, cur);
    }
    let a = FIArrayU128::from(a_bit);
    for i in x as usize..=len {
        r.arr[i - 1] = a.arr[i - 1] - a.arr[i - 2];
    }
    for i in x..=a.isqrt {
        let vi = r.arr[i as usize - 1];
        if vi == 0 {
            continue;
        }
        let n_i = n / i;
        for z in 1..=n_i / i {
            let v = vi * (a[n_i / z] - a.arr[i as usize - 2]);
            r.arr[len - z as usize] -= v;
            if z > 1 {
                r.arr[len - z as usize + 1] += v;
            }
        }
    }
    for i in 1..len {
        r.arr[i] += r.arr[i - 1];
    }
    r
}
fn lucy_fenwick_simple(x: u128) -> FIArrayU128 {
    let mut s = FIArrayU128::new(x);
    let xsqrt = s.isqrt;
    let len = s.arr.len();

    for (i, v) in FIArrayU128::keys(x).enumerate() {
        s.arr[i] = v.div_ceil(2).pow(2) - 1 + 2;
    }
    s.arr[0] = 0;
    for i in (2..len).rev() {
        let j = i & (i + 1);
        if j != 0 {
            s.arr[i] -= s.arr[j - 1];
        }
    }
    let mut s_fenwick = FenwickTreeU128(s.arr);

    let get_index = |v| -> usize {
        if v == 0 {
            unsafe { core::hint::unreachable_unchecked() };
        } else if v <= xsqrt {
            v as usize - 1
        } else {
            len - (x / v) as usize
        }
    };
    let mut sp = 0;

    let cutoff = xsqrt
        .isqrt()
        .max(2 * iroot::<3>((xsqrt / x.ilog2() as u128).pow(2)))
        | 1; // iroot::<3>(x) | 1;
    //dbg!(cutoff, iroot::<3>(x) | 1);
    for p in /* (2 <= cutoff)
        .then_some(2)
        .into_iter()
        .chain( */
        (3..=cutoff).step_by(2)
    //)
    {
        let sp1 = s_fenwick.sum(p as usize - 1);
        if sp1 == sp {
            continue;
        }

        let lim = x / p;
        let mut j = 1;
        //assert_eq!(get_index(lim), len - p);
        let mut cur = s_fenwick.sum(get_index(lim));
        while (j + 1) <= lim / (j + 1) {
            let next = s_fenwick.sum(get_index(lim / (j + 1)));
            if next != cur {
                s_fenwick.sub(len - j as usize, p * (cur - next));
                cur = next;
            }
            j += 1;
        }
        for i in (p..=lim / j).rev() {
            let next = s_fenwick.sum(i as usize - 2);
            if next != cur {
                s_fenwick.sub(get_index(p * i), p * (cur - next));
                cur = next;
            }
        }
        sp = sp1;
    }
    s.arr = s_fenwick.flatten();
    for p in /* (2 > cutoff)
        .then_some(2)
        .into_iter()
        .chain( */
        (cutoff + 2..=xsqrt).step_by(2)
    //)
    {
        let sp1 = s.arr[p as usize - 1];
        if sp1 == sp {
            continue;
        }
        let mut ip = 0;
        for i in 1..=(x / p) / p {
            ip += p;
            s.arr[len - i as usize] -= p
                * (if ip <= xsqrt {
                    s.arr[len - ip as usize]
                } else {
                    s.arr[(x / ip) as usize - 1]
                } - sp);
        }
        sp = sp1;
    }
    s
}
