use itertools::Itertools;

use crate::{
    divisor_sums, inverse_pseudo_euler_transform_fraction_i64, p300_399::e362::{mult, mult_sparse}, pseudo_euler_transform_fraction_i64, utils::{
        FIArray::{FIArray, FIArrayI64},
        math::iroot,
        primes::{
            primecount::{lucy_fenwick, lucy_fenwick_simple},
            wheel_sieve,
        },
    }
};

use paste::paste;
macro_rules! DynamicPrefixSum_impl_for {
    ($($type:ty),+) => { $(
        paste!{
            #[derive(Debug, Clone, PartialEq, Eq)]
            pub struct [<DynamicPrefixSum $type:camel>](pub Box<[$type]>, pub usize); // 1: len of flat prefix
            impl [<DynamicPrefixSum $type:camel>] {
                #[must_use]
                pub fn new(len: usize, init: $type) -> Self {
                    let mut ret = Self(vec![init; len].into_boxed_slice(), 1);
                    if init != 0 {
                        ret.construct();
                    }
                    ret
                }
                pub fn construct(&mut self) {
                    let len = self.0.len();
                    for i in 1..len {
                        let r = i + (i & (!i + 1));
                        if r <= len {
                            self.0[r - 1] += self.0[i - 1];
                        }
                    }
                }
                #[must_use]
                pub fn new_with(len: usize, init: impl FnMut(usize) -> $type) -> Self {
                    let mut ret = Self((0..len).map(init).collect_vec().into_boxed_slice(), 1);
                    ret.construct();
                    ret
                }
                #[must_use]
                pub fn sum(&self, i: usize) -> $type {
                    let k = self.1;
                    let mut i = i + 1;
                    let mut sum = 0;
                    while i > k {
                        sum += self.0[i - 1];
                        i &= i - 1;
                    }
                    if i > 0 {
                        sum += self.0[i - 1];
                    }
                    sum
                }
                pub fn add(&mut self, mut i: usize, v: $type) {
                    while i < self.0.len() {
                        self.0[i] += v;
                        i |= i + 1;
                    }
                }
                pub fn sub(&mut self, mut i: usize, v: $type) {
                    while i < self.0.len() {
                        self.0[i] -= v;
                        i |= i + 1;
                    }
                }

                pub fn inc(&mut self, mut i: usize) {
                    while i < self.0.len() {
                        self.0[i] += 1;
                        i |= i + 1;
                    }
                }
                pub fn dec(&mut self, mut i: usize) {
                    while i < self.0.len() {
                        self.0[i] -= 1;
                        i |= i + 1;
                    }
                }
                pub fn extend_flattened_prefix(&mut self, new_len: usize) {
                    let k = self.1;
                    for i in k..new_len {
                        let j = i & (i + 1);
                        if j != 0 {
                            self.0[i] += self.0[j - 1];
                        }
                    }
                    self.1 = new_len;
                }
                pub fn shrink_flattened_prefix(&mut self, new_len: usize) {
                    let k = self.1;
                    for i in (new_len..k).rev() {
                        let j = i & (i + 1);
                        if j != 0 {
                            self.0[i] -= self.0[j - 1];
                        }
                    }
                    self.1 = new_len;
                }

                #[must_use]
                pub fn flatten(self) -> Box<[$type]> {
                    let mut ret = self.0;
                    let k = self.1;
                    for i in k..ret.len() {
                        let j = i & (i + 1);
                        if j != 0 {
                            ret[i] += ret[j - 1];
                        }
                    }
                    ret
                }
            }

        }
    )+ };
}
DynamicPrefixSum_impl_for!(i64, u64, usize, isize, u128, i128);
pub type DynamicPrefixSum = DynamicPrefixSumUsize;
#[must_use]
pub fn testing(x: usize) -> FIArray {
    const LUT: [usize; 30] = [
        0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 8,
    ];

    let mut s = FIArray::new(x);
    let xsqrt = s.isqrt;
    let len = s.arr.len();

    /* for (i, v) in FIArray::keys(x).enumerate() {
        s.arr[i] = /* v - 1;  */(v + 1) >> 1;
    } */

    for (i, v) in FIArray::keys(x).enumerate() {
        //s.arr[i] = (v + 1) >> 1;
        s.arr[i] = ((v / 30) << 3) + LUT[v % 30] - 1 + 3;
    }
    s[1] = 0;
    s[2] = 1;
    s[3] = 2;
    s[4] = 2;
    s[5] = 3;
    //s.arr[0] = 0;

    let mut s_fenwick = DynamicPrefixSumUsize(s.arr, len);
    s_fenwick.shrink_flattened_prefix(1);

    let get_index = |v| -> usize {
        if v == 0 {
            unsafe { core::hint::unreachable_unchecked() };
        } else if v <= xsqrt {
            v - 1
        } else {
            len - (x / v)
        }
    };
    let mut sp = 3; //0;
    let cutoff = xsqrt.isqrt().max(iroot::<3>(x / x.ilog2() as usize)) | 1; /* xsqrt
    .isqrt()
    .max(2 * iroot::<3>((xsqrt / x.ilog2() as usize).pow(2)))
    | 1;*/
    //iroot::<3>(x) | 1;
    dbg!(cutoff, iroot::<3>(x) | 1);
    for p in /* (2 <= cutoff)
        .then_some(2)
        .into_iter()
        .chain( */
        (7..=cutoff).step_by(2)
    //)
    {
        let sp1 = s_fenwick.sum(p - 1);
        if sp1 == sp {
            continue;
        }
        s_fenwick.extend_flattened_prefix(1 + get_index(p * p));

        let lim = x / p;
        let mut j = 1;
        //assert_eq!(get_index(lim), len - p);
        let mut cur = s_fenwick.sum(get_index(lim));
        while (j + 1) <= lim / (j + 1) {
            let next = s_fenwick.sum(get_index(lim / (j + 1)));
            if next != cur {
                s_fenwick.sub(len - j, cur - next);
                cur = next;
            }
            j += 1;
        }
        for i in (p..=lim / j).rev() {
            let next = s_fenwick.sum(i - 2);
            if next != cur {
                s_fenwick.sub(get_index(p * i), cur - next);
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
        let sp1 = s.arr[p - 1];
        if sp1 == sp {
            continue;
        }
        let mut ip = 0;
        for i in 1..=(x / p) / p {
            ip += p;
            s.arr[len - i] -= if ip <= xsqrt {
                s.arr[len - ip]
            } else {
                s.arr[(x / ip) - 1]
            } - sp;
        }
        sp = sp1;
    }
    s
}
// O\left(n^\frac23 \log^{-\frac23}n \right)?
#[must_use]
pub fn testing_basic(x: usize) -> FIArray {
    let mut s = FIArray::new(x);
    let xsqrt = s.isqrt;
    let len = s.arr.len();

    for (i, v) in FIArray::keys(x).enumerate() {
        s.arr[i] = v - 1;
    }
    let mut s_fenwick = DynamicPrefixSumUsize(s.arr, len);
    s_fenwick.shrink_flattened_prefix(1);
    let get_index = |v| -> usize {
        if v == 0 {
            unsafe { core::hint::unreachable_unchecked() };
        } else if v <= xsqrt {
            v - 1
        } else {
            len - (x / v)
        }
    };
    let mut sp = 0;
    let cutoff = xsqrt.isqrt().max(iroot::<3>(x / x.ilog2() as usize));
    dbg!(cutoff, iroot::<3>(x));
    for p in 2..=cutoff {
        let sp1 = s_fenwick.sum(p - 1);
        if sp1 == sp {
            continue;
        }
        s_fenwick.extend_flattened_prefix(1 + get_index(p * p));

        let lim = x / p;
        let mut j = 1;
        let mut cur = s_fenwick.sum(get_index(lim));
        while (j + 1) <= lim / (j + 1) {
            let next = s_fenwick.sum(get_index(lim / (j + 1)));
            if next != cur {
                s_fenwick.sub(len - j, cur - next);
                cur = next;
            }
            j += 1;
        }
        for i in (p..=lim / j).rev() {
            let next = s_fenwick.sum(i - 2);
            if next != cur {
                s_fenwick.sub(get_index(p * i), cur - next);
                cur = next;
            }
        }

        sp = sp1;
    }
    s.arr = s_fenwick.flatten();
    for p in cutoff + 1..=xsqrt {
        let sp1 = s.arr[p - 1];
        if sp1 == sp {
            continue;
        }
        let mut ip = 0;
        for i in 1..=(x / p) / p {
            ip += p;
            s.arr[len - i] -= if ip <= xsqrt {
                s.arr[len - ip]
            } else {
                s.arr[(x / ip) - 1]
            } - sp;
        }
        sp = sp1;
    }
    s
}

/// The cutoff is optimized for the sparse case, losing a factor of O(\log^\frac16 n) for dense inputs compared to the optimum,
/// while gaining a factor of O(\log^\frac13 n) on sparse inputs
#[must_use]
pub fn inverse_pseudo_euler_transform_lucy_i64(mut a: FIArrayI64) -> FIArrayI64 {
    let x = a.x;
    let xsqrt = a.isqrt;
    let len = a.arr.len();
    let mut a_fenwick = DynamicPrefixSumI64(a.arr, len);
    a_fenwick.shrink_flattened_prefix(1);
    a_fenwick.dec(0);
    let get_index = |v| -> usize {
        if v == 0 {
            unsafe { core::hint::unreachable_unchecked() };
        } else if v <= xsqrt {
            v - 1
        } else {
            len - (x / v)
        }
    };
    let mut sp = 0;
    let cutoff = xsqrt.isqrt().max(iroot::<3>(x / x.ilog2() as usize));
    for p in 2..=cutoff {
        let sp1 = a_fenwick.sum(p - 1);
        if sp1 == sp {
            continue;
        }

        a_fenwick.extend_flattened_prefix(1 + get_index(p * p));

        let w = sp1 - sp;
        let lim = x / p;
        let mut j = 1;
        let mut cur = a_fenwick.sum(get_index(lim));
        while (j + 1) <= lim / (j + 1) {
            let next = a_fenwick.sum(get_index(lim / (j + 1)));
            if next != cur {
                a_fenwick.sub(len - j, w * (cur - next));
                cur = next;
            }
            j += 1;
        }
        for i in (p..=lim / j).rev() {
            let next = a_fenwick.sum(i - 2);
            if next != cur {
                a_fenwick.sub(get_index(p * i), w * (cur - next));
                cur = next;
            }
        }

        sp = sp1;
    }
    a.arr = a_fenwick.flatten();
    for p in cutoff + 1..=xsqrt {
        let sp1 = a.arr[p - 1];
        if sp1 == sp {
            continue;
        }
        let w = sp1 - sp;
        let mut ip = 0;
        for i in 1..=(x / p) / p {
            ip += p;
            a.arr[len - i] -= w
                * (if ip <= xsqrt {
                    a.arr[len - ip]
                } else {
                    a.arr[(x / ip) - 1]
                } - sp);
        }
        sp = sp1;
    }
    a
}
#[must_use]
pub fn pseudo_euler_transform_lucy_i64(mut a: FIArrayI64) -> FIArrayI64 {
    let x = a.x;
    let xsqrt = a.isqrt;
    let len = a.arr.len();
    let cutoff = xsqrt.isqrt().max(iroot::<3>(x / x.ilog2() as usize));

    let mut sp = a.arr[xsqrt - 1];

    for p in (cutoff + 1..=xsqrt).rev() {
        let sp1 = a.arr[p - 2];
        if sp1 == sp {
            continue;
        }
        let w = sp - sp1;
        //let mut ip = 0;
        for i in (1..=(x / p) / p).rev() {
            a.arr[len - i] += w
                * (if i * p <= xsqrt {
                    a.arr[len - i * p]
                } else {
                    a.arr[(x / (i * p)) - 1]
                } - sp1);
        }
        sp = sp1;
    }

    let mut a_fenwick = DynamicPrefixSumI64(a.arr, len);
    let get_index = |v| -> usize {
        if v == 0 {
            unsafe { core::hint::unreachable_unchecked() };
        } else if v <= xsqrt {
            v - 1
        } else {
            len - (x / v)
        }
    };
    for p in (2..=cutoff).rev() {
        let sp1 = a_fenwick.sum(p - 2);
        if sp1 == sp {
            continue;
        }

        a_fenwick.shrink_flattened_prefix(1 + get_index(p * p));

        let w = sp - sp1;
        let lim = x / p;
        let mut prev = sp1;
        let mut i = p;
        while i <= lim / i {
            let cur = a_fenwick.sum(i - 1);
            if cur != prev {
                a_fenwick.add(get_index(i * p), w * (cur - prev));
                prev = cur;
            }
            i += 1;
        }
        for j in (1..=lim / i).rev() {
            let cur = a_fenwick.sum(get_index(lim / j));
            if cur != prev {
                a_fenwick.add(len - j, w * (cur - prev));
                prev = cur;
            }
        }
        sp = sp1;
    }
    a_fenwick.shrink_flattened_prefix(1);
    if a_fenwick.sum(0) == 0 {
        a_fenwick.inc(0);
    }
    a.arr = a_fenwick.flatten();
    a
}

pub fn main() {
    const N: usize = 1e14 as _;
    println!("Entering {} main", file!());

    let start = std::time::Instant::now();
    let s1 = lucy_fenwick_simple(N);
    println!("{} | {:?}", s1[N], start.elapsed());
    let start = std::time::Instant::now();
    let s1 = lucy_fenwick(N);
    println!("{} | {:?}", s1[N], start.elapsed());
    let start = std::time::Instant::now();
    let s1 = testing(N);
    println!("{} | {:?}", s1[N], start.elapsed());
    let start = std::time::Instant::now();
    let s1 = testing_basic(N);
    println!("{} | {:?}", s1[N], start.elapsed());
    
    let start = std::time::Instant::now();
    let s2 = log_zeta_fast_partial_flatten(N);
    println!("{} | {:?}", s2[N], start.elapsed());
    //assert_eq!(s1, s2,);
    
    let start = std::time::Instant::now();
    let s2 = log_zeta_2_partial_flatten(N);
    println!("{} | {:?}", s2[N], start.elapsed());
    
    let start = std::time::Instant::now();
    let s2 = log_zeta_partial_flatten(N);
    println!("{} | {:?}", s2[N], start.elapsed());
    /* let start = std::time::Instant::now();
    let s = legendre_test(N);
    println!("{} | {:?}", s, start.elapsed()); // 2^40: 1.6238139s, 2^49: 114.6576979s

    let start = std::time::Instant::now();
    let s = lucy(N)[N];
    println!("{} | {:?}", s, start.elapsed()); */

    /* let mut chi4 = FIArrayI64::new(N);
    for (i, v) in FIArrayI64::keys(N).enumerate() {
        chi4.arr[i] = [0, 1, 1, 0][v & 3];
    }
    assert_eq!(
        inverse_pseudo_euler_transform_lucy_i64(chi4.clone()),
        inverse_pseudo_euler_transform_fraction_i64(chi4.clone())
    );
    assert_eq!(
        pseudo_euler_transform_lucy_i64(inverse_pseudo_euler_transform_lucy_i64(chi4.clone())),
        pseudo_euler_transform_fraction_i64(inverse_pseudo_euler_transform_lucy_i64(chi4))
    ); */
    //dbg!(testing_basic(1000));
    assert_eq!(
        lucy_pet_slow(/* dbg! */ (lucy_ipet_slow(FIArray::unit(1000)))),
        FIArray::unit(1000)
    );
    println!("Exiting {} main", file!());
}
// O(n^\frac34 \log ^{-1} n) with most of the time spent on cheap O(1) sum calls, much faster than without flattening
fn legendre_test(x: usize) -> usize {
    let s = FIArray::unit(x);
    let xsqrt = s.isqrt;
    let len = s.arr.len();
    let mut s_fenwick = DynamicPrefixSumUsize(s.arr, len);
    s_fenwick.shrink_flattened_prefix(1);

    let get_index = |v| -> usize {
        if v == 0 {
            unsafe { core::hint::unreachable_unchecked() };
        } else if v <= xsqrt {
            v - 1
        } else {
            len - (x / v)
        }
    };
    let primes = wheel_sieve(s.isqrt as u64);
    for &p in &primes {
        let p = p as usize;
        //s.sparse_mul_at_most_one(p as _, 1);
        s_fenwick.extend_flattened_prefix(p);
        let lim = x / p;
        let mut j = 1;
        let mut cur = s_fenwick.sum(get_index(lim));
        while (j + 1) <= lim / (j + 1) {
            let next = s_fenwick.sum(get_index(lim / (j + 1)));
            if next != cur {
                s_fenwick.sub(len - j, cur - next);
                cur = next;
            }
            j += 1;
        }
        for i in (2..=lim / j).rev() {
            let next = s_fenwick.sum(i - 2);
            if next != cur {
                s_fenwick.sub(get_index(p * i), cur - next);
                cur = next;
            }
        }
        if cur != 0 {
            s_fenwick.sub(p - 1, cur);
        }
    }
    s_fenwick.sum(len - 1) + primes.len() - 1
}

#[must_use]
pub fn lucy_ipet_slow(mut s: FIArray) -> FIArray {
    s.adjacent_difference();
    s.arr[0] = 0;
    s.partial_sum();
    let keys = FIArray::keys(s.x).collect_vec();
    let mut sp1 = 0;
    for p in 2..=s.isqrt {
        let sp = s.arr[p - 1];
        if sp == sp1 {
            continue;
        }
        let fp = sp - sp1;
        for (i, &v) in keys.iter().enumerate().rev() {
            if v < p * p {
                break;
            }
            s.arr[i] -= fp * (s[v / p] - sp1);
        }
        sp1 = sp;
    }
    s
}
#[must_use]
pub fn lucy_pet_slow(mut s: FIArray) -> FIArray {
    let keys = FIArray::keys(s.x).collect_vec();
    let mut sp = s.arr[s.isqrt - 1];
    for p in (2..=s.isqrt).rev() {
        let sp1 = s.arr[p - 2];
        if sp == sp1 {
            continue;
        }
        let fp = sp - sp1;
        for (i, &v) in keys.iter().enumerate().skip(s.get_index(p * p)) {
            s.arr[i] += fp * (s[v / p] - sp1);
        }
        sp = sp1;
    }
    s.adjacent_difference();
    s.arr[0] = 1;
    s.partial_sum();
    s
}
#[must_use]
pub fn log_zeta_2_partial_flatten(n: usize) -> FIArray {
    const fn inv_odd(mut k: usize) -> usize {
        let mut exp = (1u64 << 63) - 1;

        let mut r: usize = 1;
        while exp > 1 {
            r = r.overflowing_mul(k).0;
            k = k.overflowing_mul(k).0;
            exp >>= 1;
        }
        r.overflowing_mul(k).0
    }

    const INVS: [usize; 5] = [0, 4, 2, inv_odd(3) << 2, 1];
    let mut zeta = FIArray::unit(n);
    let rt = zeta.isqrt;
    let len = zeta.arr.len();

    let x = iroot::<5>(n) + 1;
    // remove contributions of small primes
    /* for p in 2..x {
        if zeta.bit.sum(p - 1) == 1 {
            //not prime
            continue;
        }
        ret.arr[p - 1] = 1;
        zeta.sparse_mul_at_most_one(p, 1);
    }
    zeta.bit.dec(0);
     */
    let mut zeta_fenwick = DynamicPrefixSumUsize(zeta.arr, len);
    zeta_fenwick.shrink_flattened_prefix(1);
    zeta_fenwick.dec(0);
    let get_index = |v| -> usize {
        if v == 0 {
            unsafe { core::hint::unreachable_unchecked() };
        } else if v <= rt {
            v - 1
        } else {
            len - (n / v)
        }
    };
    let mut sp = 0;
    for p in 2..x {
        let sp1 = zeta_fenwick.sum(p - 1);
        if sp1 == sp {
            continue;
        }
        zeta_fenwick.extend_flattened_prefix(1 + get_index(p * p));

        let lim = n / p;
        let mut j = 1;
        let mut cur = zeta_fenwick.sum(get_index(lim));
        while (j + 1) <= lim / (j + 1) {
            let next = zeta_fenwick.sum(get_index(lim / (j + 1)));
            if next != cur {
                zeta_fenwick.sub(len - j, cur - next);
                cur = next;
            }
            j += 1;
        }
        for i in (p..=lim / j).rev() {
            let next = zeta_fenwick.sum(i - 2);
            if next != cur {
                zeta_fenwick.sub(get_index(p * i), cur - next);
                cur = next;
            }
        }

        sp = sp1;
    }
    zeta.arr = zeta_fenwick.flatten();

    //let zeta = FIArray::from(zeta);
    let mut ret = FIArray::new(n);
    zeta.adjacent_difference();
    ret.arr[..x - 1].copy_from_slice(&zeta.arr[..x - 1]);
    zeta.arr[..x - 1].fill(0);
    zeta.partial_sum();

    // zeta now equals zeta_t - 1
    // compute log(zeta_t) using log(x + 1) = x - x^2 / 2 + x^3 / 3 - x^4 / 4
    // in order to not have to deal with rational numbers, we compute 12 * log(zeta_t)
    // and adjust later
    for i in x..=len {
        ret.arr[i - 1] = zeta.arr[i - 1] * INVS[1];
    }
    let zeta_2 = mult(&zeta, &zeta);
    let mut x_pow = x * x;

    for i in ret.get_index(x_pow)..=len {
        ret.arr[i - 1] -= zeta_2.arr[i - 1] * INVS[2];
    }
    let zeta_3 = mult(&zeta, &zeta_2);
    x_pow *= x;
    for i in ret.get_index(x_pow)..=len {
        ret.arr[i - 1] += zeta_3.arr[i - 1] * INVS[3];
    }
    let zeta_4 = mult_sparse(&zeta, &zeta_3);
    x_pow *= x;
    for i in ret.get_index(x_pow)..=len {
        ret.arr[i - 1] -= zeta_4.arr[i - 1] * INVS[4];
    }
    // correction phase: get rid of contributions of prime powers
    for i in (x..len).rev() {
        ret.arr[i] -= ret.arr[i - 1];
    }

    for x in x..=rt {
        let v = ret.arr[x - 1] / INVS[1];
        if v == 0 {
            continue;
        }
        let mut e = 1;
        let mut px = x;
        while px <= n / x {
            e += 1;
            px *= x;

            ret[px] -= INVS[e];
        }
    }
    for i in 1..x - 1 {
        ret.arr[i] += ret.arr[i - 1];
    }
    for i in x - 1..len {
        ret.arr[i] /= INVS[1];
        ret.arr[i] += ret.arr[i - 1];
    }
    ret
}
#[must_use]
pub fn log_zeta_partial_flatten(n: usize) -> FIArray {
    const fn inv_odd(mut k: usize) -> usize {
        let mut exp = (1u64 << 63) - 1;

        let mut r: usize = 1;
        while exp > 1 {
            r = r.overflowing_mul(k).0;
            k = k.overflowing_mul(k).0;
            exp >>= 1;
        }
        r.overflowing_mul(k).0
    }

    const INVS: [usize; 4] = [0, 2, 1, inv_odd(3) << 1];
    let mut zeta = FIArray::unit(n);
    let rt = zeta.isqrt;
    let len = zeta.arr.len();
    let x = iroot::<4>(n) + 1;
    // remove contributions of small primes
    /* let mut ret = FIArray::new(n);
    for p in 2..x {
        if zeta.bit.sum(p - 1) == 1 {
            //not prime
            continue;
        }
        ret.arr[p - 1] = 1;
        zeta.sparse_mul_at_most_one(p, 1);
    }
    zeta.bit.dec(0);
    let mut zeta = FIArray::from(zeta); */
    let mut zeta_fenwick = DynamicPrefixSumUsize(zeta.arr, len);
    zeta_fenwick.shrink_flattened_prefix(1);
    zeta_fenwick.dec(0);
    let get_index = |v| -> usize {
        if v == 0 {
            unsafe { core::hint::unreachable_unchecked() };
        } else if v <= rt {
            v - 1
        } else {
            len - (n / v)
        }
    };
    let mut sp = 0;
    for p in 2..x {
        let sp1 = zeta_fenwick.sum(p - 1);
        if sp1 == sp {
            continue;
        }
        zeta_fenwick.extend_flattened_prefix(1 + get_index(p * p));

        let lim = n / p;
        let mut j = 1;
        let mut cur = zeta_fenwick.sum(get_index(lim));
        while (j + 1) <= lim / (j + 1) {
            let next = zeta_fenwick.sum(get_index(lim / (j + 1)));
            if next != cur {
                zeta_fenwick.sub(len - j, cur - next);
                cur = next;
            }
            j += 1;
        }
        for i in (p..=lim / j).rev() {
            let next = zeta_fenwick.sum(i - 2);
            if next != cur {
                zeta_fenwick.sub(get_index(p * i), cur - next);
                cur = next;
            }
        }

        sp = sp1;
    }
    zeta.arr = zeta_fenwick.flatten();

    let mut ret = FIArray::new(n);
    zeta.adjacent_difference();
    ret.arr[..x - 1].copy_from_slice(&zeta.arr[..x - 1]);
    zeta.arr[..x - 1].fill(0);
    zeta.partial_sum();

    // zeta now equals zeta_t - 1
    // compute log(zeta_t) using log(x + 1) = x^3 / 3 - x^2 / 2 + x
    // in order to not have to deal with rational numbers, we compute 6 * log(zeta_t)
    // and adjust later

    for i in x..=len {
        ret.arr[i - 1] = zeta.arr[i - 1] * INVS[1];
    }

    let pa = {
        let mut vec = vec![];
        for i in x..=rt {
            if zeta.arr[i - 1] != zeta.arr[i - 2] {
                vec.push(i);
            }
        }
        vec.push(rt + 1);
        vec
    };
    let va = &pa[..pa.len() - 1];
    let mut pow_zeta = //mult(&zeta, &zeta);
    {
        let mut res = FIArray::new(n);
        let mut r = va.len();
        let mut l = 0;
        for &x in va {
            res[x * x] += 1;

            l += 1;
            let Nx = n / x;
            while r > l && Nx / pa[r - 1] < r - 1 {
                r -= 1;
            }
            if r < l {
                r = l;
            }
            let X = rt / x;
            let Nx = n / x;

            let mut i = l;
            /* while i != r {
                res[x * pa[i]] += 2;
                i += 1;
            } */
            while i != r && pa[i] <= X {
                res.arr[x * pa[i] - 1] += 2 ;
                i += 1;
            }
            while i != r {
                res.arr[len - Nx / pa[i]] += 2 ;
                i += 1;
            }
            
            if r != 0 && pa[r] <= Nx {
                res[x * pa[r]] -= zeta.arr[pa[r - 1] - 1] * 2;
            }
        }
        res.partial_sum();
        r = va.len();
        l = 0;
        for &x in va {
            l += 1;
            let Nx = n / x;
            while r > l && Nx / pa[r - 1] < r - 1 {
                r -= 1;
            }
            if r < l {
                r = l;
            }
            /* for i in (1..=Nx / pa[r]).rev() {
                res[n / i] += 2 * zeta[Nx / i];
            } */
            let mut i = Nx / pa[r];
            let X = rt / x;
            while i > X {
                res.arr[len - i] += 2 * zeta.arr[Nx / i - 1];
                i -= 1;
            }
            while i > 0 {
                res.arr[len - i] += 2 * zeta.arr[len - x * i];
                i -= 1;
            }
        }
        res
    };
    /* let ind = pow_zeta.get_index(x.pow(2));
    assert!(pow_zeta.arr[..ind].iter().all(|&e| e == 0));
    dbg!(ind, x.pow(2), rt - 1, len, len - ind); */

    for i in rt + 1..=len {
        ret.arr[i - 1] -= pow_zeta.arr[i - 1] * INVS[2];
    }
    {
        //pow_zeta = mult_sparse(&zeta, &pow_zeta);

        zeta.arr[rt..].fill(0);
        for &i in va {
            for j in 1..=rt / i {
                //zeta[n / j] += pow_zeta[n / (i * j)];
                zeta.arr[len - j] += pow_zeta.arr[len - i * j];
            }
        }
        //zeta.arr[..rt].fill(0);
        core::mem::swap(&mut pow_zeta, &mut zeta);
    }
    let ind = pow_zeta.get_index(x.pow(3));
    //dbg!(ind, len, len - ind);
    //assert!(pow_zeta.arr[..ind].iter().all(|&e| e == 0)); */
    for i in ind + 1..=len {
        ret.arr[i - 1] += pow_zeta.arr[i - 1] * INVS[3];
    }

    // correction phase: get rid of contributions of prime powers
    for i in (x..len).rev() {
        ret.arr[i] -= ret.arr[i - 1];
    }

    for x in x..=rt {
        let v = ret.arr[x - 1] / INVS[1];
        if v == 0 {
            continue;
        }
        let mut e = 1;
        let mut px = x;
        while px <= n / x {
            e += 1;
            px *= x;

            ret[px] -= INVS[e];
        }
    }
    for i in 1..x - 1 {
        ret.arr[i] += ret.arr[i - 1];
    }
    for i in x - 1..len {
        ret.arr[i] /= INVS[1];
        ret.arr[i] += ret.arr[i - 1];
    }
    ret
}

fn log_zeta_fast_partial_flatten(n: usize) -> FIArray {
    fn mult_correction(d: &FIArray, primes: &[usize]) -> FIArray {
        struct Correction(FIArray, usize);
        impl Correction {
            fn fill(&mut self, primes: &[usize], lim: usize, x: usize, y: usize) {
                self.0[x] += y;
                self.1 += 1;
                for (i, &p) in primes.iter().enumerate() {
                    if p > lim / p {
                        break;
                    }
                    let mut pp = p * p;
                    let mut new_lim = lim / pp;
                    for e in 2.. {
                        let hp = 1 << (e - 2);
                        if hp != 0 {
                            self.fill(&primes[i + 1..], new_lim, x * pp, y * hp);
                        }
                        if p > new_lim {
                            break;
                        }
                        pp *= p;
                        new_lim /= p;
                    }
                }
            }
        }
        let mut correction = Correction(FIArray::new(d.x), 0);
        correction.fill(primes, d.x, 1, 1);
        for i in 1..correction.0.arr.len() {
            correction.0.arr[i] += correction.0.arr[i - 1];
        }
        dbg!(correction.1);
        mult_sparse(d, &correction.0)
    }
    const fn inv_odd(mut k: usize) -> usize {
        let mut exp = (1u64 << 63) - 1;

        let mut r: usize = 1;
        while exp > 1 {
            r = r.overflowing_mul(k).0;
            k = k.overflowing_mul(k).0;
            exp >>= 1;
        }
        r.overflowing_mul(k).0
    }

    const INVS: [usize; 4] = [0, 2, 1, inv_odd(3) << 1];
    let start = std::time::Instant::now();
    let mut zeta = FIArray::unit(n);
    let rt = zeta.isqrt;
    let len = zeta.arr.len();

    let mut primes = vec![];
    let x = iroot::<4>(n) + 1;
    let mut zeta_fenwick = DynamicPrefixSumUsize(zeta.arr, len);
    zeta_fenwick.shrink_flattened_prefix(1);
    zeta_fenwick.dec(0);
    let get_index = |v| -> usize {
        if v == 0 {
            unsafe { core::hint::unreachable_unchecked() };
        } else if v <= rt {
            v - 1
        } else {
            len - (n / v)
        }
    };
    let mut sp = 0;
    for p in 2..x {
        let sp1 = zeta_fenwick.sum(p - 1);
        if sp1 == sp {
            continue;
        }
        primes.push(p);
        zeta_fenwick.extend_flattened_prefix(1 + get_index(p * p));

        let lim = n / p;
        let mut j = 1;
        let mut cur = zeta_fenwick.sum(get_index(lim));
        while (j + 1) <= lim / (j + 1) {
            let next = zeta_fenwick.sum(get_index(lim / (j + 1)));
            if next != cur {
                zeta_fenwick.sub(len - j, cur - next);
                cur = next;
            }
            j += 1;
        }
        for i in (p..=lim / j).rev() {
            let next = zeta_fenwick.sum(i - 2);
            if next != cur {
                zeta_fenwick.sub(get_index(p * i), cur - next);
                cur = next;
            }
        }

        sp = sp1;
    }
    zeta.arr = zeta_fenwick.flatten();
    dbg!(start.elapsed());
    let mut zeta_2 = divisor_sums(n);
    dbg!(start.elapsed());
    zeta_2 = mult_correction(&zeta_2, &primes);
    dbg!(start.elapsed());

    let mut zeta_2_fenwick = DynamicPrefixSumUsize(zeta_2.arr, len);
    zeta_2_fenwick.shrink_flattened_prefix(1);
    zeta_2_fenwick.dec(0);
    for &p in &primes {
        zeta_2_fenwick.extend_flattened_prefix(1 + get_index(p * p));

        let lim = n / p;
        let mut j = 1;
        let mut cur = zeta_2_fenwick.sum(get_index(lim));
        while (j + 1) <= lim / (j + 1) {
            let next = zeta_2_fenwick.sum(get_index(lim / (j + 1)));
            if next != cur {
                zeta_2_fenwick.sub(len - j, 2 * (cur - next));
                cur = next;
            }
            j += 1;
        }
        for i in (p..=lim / j).rev() {
            let next = zeta_2_fenwick.sum(i - 2);
            if next != cur {
                zeta_2_fenwick.sub(get_index(p * i), 2 * (cur - next));
                cur = next;
            }
        }
    }
    zeta_2.arr = zeta_2_fenwick.flatten();
    dbg!(start.elapsed());

    let mut ret = FIArray::new(n);
    zeta.adjacent_difference();
    ret.arr[..x - 1].copy_from_slice(&zeta.arr[..x - 1]);
    zeta.arr[..x - 1].fill(0);
    zeta.partial_sum();

    zeta_2.adjacent_difference();
    zeta_2.arr[..x - 1].fill(0);
    zeta_2.partial_sum();
    for i in 0..len {
        zeta_2.arr[i] -= 2 * zeta.arr[i];
    }
    dbg!(start.elapsed());

    // zeta now equals zeta_t - 1, and zeta_2 (zeta_t - 1)^2
    // compute log(zeta_t) using log(x + 1) = x^3 / 3 - x^2 / 2 + x
    // in order to not have to deal with rational numbers, we compute 6 * log(zeta_t)
    // and adjust later

    for i in x..=len {
        ret.arr[i - 1] = zeta.arr[i - 1] * INVS[1];
    }

    let pa = {
        let mut vec = vec![];
        for i in x..=rt {
            if zeta.arr[i - 1] != zeta.arr[i - 2] {
                vec.push(i);
            }
        }
        vec.push(rt + 1);
        vec
    };
    let va = &pa[..pa.len() - 1];

    /* let ind = zeta_2.get_index(x.pow(2));
    assert!(pow_zeta.arr[..ind].iter().all(|&e| e == 0));
    dbg!(ind, x.pow(2), rt - 1, len, len - ind); */

    for i in rt + 1..=len {
        ret.arr[i - 1] -= zeta_2.arr[i - 1] * INVS[2];
    }

    //pow_zeta = mult_sparse(&zeta, &pow_zeta);

    zeta.arr[rt..].fill(0);
    for &i in va {
        for j in 1..=rt / i {
            //zeta[n / j] += pow_zeta[n / (i * j)];
            zeta.arr[len - j] += zeta_2.arr[len - i * j];
        }
    }
    //zeta.arr[..rt].fill(0);
    let zeta_3 = zeta;

    let ind = zeta_3.get_index(x.pow(3));
    //dbg!(ind, len, len - ind);
    //assert!(pow_zeta.arr[..ind].iter().all(|&e| e == 0)); */
    for i in ind + 1..=len {
        ret.arr[i - 1] += zeta_3.arr[i - 1] * INVS[3];
    }
    dbg!(start.elapsed());
    // correction phase: get rid of contributions of prime powers
    for i in (x..len).rev() {
        ret.arr[i] -= ret.arr[i - 1];
    }

    for x in x..=rt {
        let v = ret.arr[x - 1] / INVS[1];
        if v == 0 {
            continue;
        }
        let mut e = 1;
        let mut px = x;
        while px <= n / x {
            e += 1;
            px *= x;

            ret[px] -= INVS[e];
        }
    }
    for i in 1..x - 1 {
        ret.arr[i] += ret.arr[i - 1];
    }
    for i in x - 1..len {
        ret.arr[i] /= INVS[1];
        ret.arr[i] += ret.arr[i - 1];
    }
    ret
}
