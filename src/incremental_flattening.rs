use itertools::Itertools;

use crate::utils::{
    FIArray::FIArray,
    math::iroot,
    primes::primecount::{lucy_fenwick, lucy_fenwick_simple},
};

#[derive(Debug, Clone, PartialEq, Eq)]
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
        if i < k {
            return self.0[i];
        }
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
    pub fn partial_flatten(&mut self, new_k: usize) {
        let k = self.1;
        for i in k..new_k {
            let j = i & (i + 1);
            if j != 0 {
                self.0[i] += self.0[j - 1];
            }
        }
        self.1 = new_k;
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
    for i in (2..len).rev() {
        let j = i & (i + 1);
        if j != 0 {
            s.arr[i] -= s.arr[j - 1];
        }
    }
    let mut s_fenwick = DynamicPrefixSum(s.arr, 1);

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

    let cutoff = xsqrt
        .isqrt()
        .max(2 * iroot::<3>((xsqrt / x.ilog2() as usize).pow(2)))
        | 1; //iroot::<3>(x) | 1;
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
        s_fenwick.partial_flatten(1 + get_index(p * p));

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

pub fn main() {
    const N: usize = 1e17 as _;

    let start = std::time::Instant::now();
    let s1 = lucy_fenwick_simple(N);
    println!("{} | {:?}", s1[N], start.elapsed());
    let start = std::time::Instant::now();
    let s1 = lucy_fenwick(N);
    println!("{} | {:?}", s1[N], start.elapsed());
    let start = std::time::Instant::now();
    let s2 = testing(N);
    println!("{} | {:?}", s2[N], start.elapsed());
    //assert_eq!(s1, s2,);
}
