use itertools::Itertools;

use crate::{
    inverse_pseudo_euler_transform_fraction_i64, inverse_pseudo_euler_transform_i64,
    pseudo_euler_transform_fraction_i64,
    utils::{
        FIArray::{FIArray, FIArrayI64},
        math::iroot,
        primes::{
            primecount::{lucy, lucy_fenwick, lucy_fenwick_simple},
            wheel_sieve,
        },
    },
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
/*#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DynamicPrefixSum(pub Box<[usize]>, usize); // 1: len of flat prefix
impl DynamicPrefixSum {
    #[must_use]
    pub fn new(len: usize, init: usize) -> Self {
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
    pub fn new_with(len: usize, init: impl FnMut(usize) -> usize) -> Self {
        let mut ret = Self((0..len).map(init).collect_vec().into_boxed_slice(), 1);
        ret.construct();
        ret
    }
    #[must_use]
    pub fn sum(&self, i: usize) -> usize {
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
    pub fn add(&mut self, mut i: usize, v: usize) {
        while i < self.0.len() {
            self.0[i] += v;
            i |= i + 1;
        }
    }
    pub fn sub(&mut self, mut i: usize, v: usize) {
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
    pub fn flatten(self) -> Box<[usize]> {
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
*/
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

// TODO: implement unflatten and reverse step of min25 sieve for totally multiplicative functions, isomorphic changes to sparse_mul_unlimited
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
    let s2 = testing_basic(N);
    println!("{} | {:?}", s2[N], start.elapsed());
    //assert_eq!(s1, s2,);
    /* let start = std::time::Instant::now();
    let s = legendre_test(N);
    println!("{} | {:?}", s, start.elapsed()); // 2^40: 1.6238139s, 2^49: 114.6576979s

    let start = std::time::Instant::now();
    let s = lucy(N)[N];
    println!("{} | {:?}", s, start.elapsed()); */

    let mut chi4 = FIArrayI64::new(N);
    for (i, v) in FIArrayI64::keys(N).enumerate() {
        chi4.arr[i] = [0, 1, 1, 0][v & 3];
    }
    assert_eq!(
        inverse_pseudo_euler_transform_lucy_i64(chi4.clone()),
        inverse_pseudo_euler_transform_fraction_i64(chi4.clone())
    );
    assert_eq!(
        pseudo_euler_transform_lucy_i64(inverse_pseudo_euler_transform_lucy_i64(chi4.clone())),
        chi4
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
