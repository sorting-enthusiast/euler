use std::time::Instant;

use itertools::Itertools;
use paste::paste;

macro_rules! FenwickTree_impl_for {
    ($($type:ty),+) => { $(
        paste!{
            #[derive(Debug, Clone, PartialEq, Eq)]
            pub struct FenwickTree {
                pub len: usize,
                pub data: Box<[$type]>,
            }

            impl FenwickTree {
                #[inline(always)]
                const fn hole(k: usize) -> usize {
                    k + (k >> 10)
                }

                #[must_use]
                pub fn new(len: usize, init: $type) -> Self {
                    let data_len = Self::hole(len);
                    let mut ret = Self {
                        len,
                        data: vec![init; data_len].into_boxed_slice(),
                    };
                    if init != 0 {
                        ret.construct();
                    }
                    ret
                }

                pub fn construct(&mut self) {
                    let len = self.len;
                    for i in 1..len {
                        let r = i + (i & (!i + 1));
                        if r <= len {
                            let val = self.data[Self::hole(i) - 1];
                            self.data[Self::hole(r) - 1] += val;
                        }
                    }
                }

                #[must_use]
                pub fn new_with(len: usize, mut init: impl FnMut(usize) -> $type) -> Self {
                    let data_len = Self::hole(len);
                    let mut data = vec![0; data_len].into_boxed_slice();
                    for i in 0..len {
                        data[Self::hole(i + 1) - 1] = init(i);
                    }
                    let mut ret = Self { len, data };
                    ret.construct();
                    ret
                }

                #[must_use]
                pub fn sum(&self, i: usize) -> $type {
                    let mut i = i + 1;
                    let mut sum = 0;
                    while i > 0 {
                        sum += self.data[Self::hole(i) - 1];
                        i &= i - 1;
                    }
                    sum
                }

                pub fn add(&mut self, mut i: usize, v: $type) {
                    while i < self.len {
                        self.data[Self::hole(i + 1) - 1] += v;
                        i |= i + 1;
                    }
                }

                pub fn sub(&mut self, mut i: usize, v: $type) {
                    while i < self.len {
                        self.data[Self::hole(i + 1) - 1] -= v;
                        i |= i + 1;
                    }
                }

                pub fn inc(&mut self, mut i: usize) {
                    while i < self.len {
                        self.data[Self::hole(i + 1) - 1] += 1;
                        i |= i + 1;
                    }
                }

                pub fn dec(&mut self, mut i: usize) {
                    while i < self.len {
                        self.data[Self::hole(i + 1) - 1] -= 1;
                        i |= i + 1;
                    }
                }

                #[must_use]
                pub fn flatten(self) -> Box<[$type]> {
                    let mut ret = vec![0; self.len].into_boxed_slice();
                    for i in 0..self.len {
                        ret[i] = self.data[Self::hole(i + 1) - 1];
                        let j = i & (i + 1);
                        if j != 0 {
                            ret[i] += ret[j - 1];
                        }
                    }
                    ret
                }
            }

            #[derive(Debug, Clone, PartialEq, Eq)]
            pub struct [<FenwickTree $type:camel Mod>]<const MOD: $type> {
                pub len: usize,
                pub data: Box<[$type]>,
            }

            impl<const MOD: $type> [<FenwickTree $type:camel Mod>]<MOD> {
                #[inline(always)]
                const fn hole(k: usize) -> usize {
                    k + (k >> 10)
                }

                #[must_use]
                pub fn new(len: usize, init: $type) -> Self {
                    let data_len = Self::hole(len);
                    let mut ret = Self {
                        len,
                        data: vec![init % MOD; data_len].into_boxed_slice(),
                    };
                    if init != 0 {
                        ret.construct();
                    }
                    ret
                }

                pub fn construct(&mut self) {
                    let len = self.len;
                    for i in 1..len {
                        let r = i + (i & (!i + 1));
                        if r <= len {
                            let val = self.data[Self::hole(i) - 1];
                            let target = Self::hole(r) - 1;
                            self.data[target] += val;
                            if self.data[target] >= MOD {
                                self.data[target] -= MOD;
                            }
                        }
                    }
                }

                #[must_use]
                pub fn new_with(len: usize, mut init: impl FnMut(usize) -> $type) -> Self {
                    let data_len = Self::hole(len);
                    let mut data = vec![0; data_len].into_boxed_slice();
                    for i in 0..len {
                        data[Self::hole(i + 1) - 1] = init(i) % MOD;
                    }
                    let mut ret = Self { len, data };
                    ret.construct();
                    ret
                }

                #[must_use]
                pub fn sum(&self, i: usize) -> $type {
                    let mut i = i + 1;
                    let mut sum = 0;
                    while i > 0 {
                        sum += self.data[Self::hole(i) - 1];
                        if sum >= MOD {
                            sum -= MOD;
                        }
                        i &= i - 1;
                    }
                    sum
                }

                pub fn add(&mut self, mut i: usize, v: $type) {
                    let v = v % MOD;
                    while i < self.len {
                        let idx = Self::hole(i + 1) - 1;
                        self.data[idx] += v;
                        if self.data[idx] >= MOD {
                            self.data[idx] -= MOD;
                        }
                        i |= i + 1;
                    }
                }

                pub fn sub(&mut self, mut i: usize, v: $type) {
                    while i < self.len {
                        let idx = Self::hole(i + 1) - 1;
                        self.data[idx] += MOD - v;
                        if self.data[idx] >= MOD {
                            self.data[idx] -= MOD;
                        }
                        i |= i + 1;
                    }
                }

                pub fn inc(&mut self, mut i: usize) {
                    while i < self.len {
                        let idx = Self::hole(i + 1) - 1;
                        self.data[idx] += 1;
                        if self.data[idx] >= MOD {
                            self.data[idx] -= MOD;
                        }
                        i |= i + 1;
                    }
                }

                pub fn dec(&mut self, mut i: usize) {
                    while i < self.len {
                        let idx = Self::hole(i + 1) - 1;
                        self.data[idx] += MOD - 1;
                        if self.data[idx] >= MOD {
                            self.data[idx] -= MOD;
                        }
                        i |= i + 1;
                    }
                }

                #[must_use]
                pub fn flatten(self) -> Box<[$type]> {
                    let mut ret = vec![0; self.len].into_boxed_slice();
                    for i in 0..self.len {
                        ret[i] = self.data[Self::hole(i + 1) - 1];
                        let j = i & (i + 1);
                        if j != 0 {
                            ret[i] += ret[j - 1];
                            if ret[i] >= MOD {
                                ret[i] -= MOD;
                            }
                        }
                    }
                    ret
                }
            }
        }
    )+ };
}

FenwickTree_impl_for!(usize);

use crate::{divisor_sums, icy_mult, p300_399::e362::{mult, mult_sparse}, utils::{
    FIArray::FIArray, fast_divisor_sums::divisor_summatory, math::iroot, primes::{log_zeta::{log_zeta, log_zeta_2}, primecount::lucy_fenwick}
}};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DirichletFenwickHole {
    pub x: usize,
    pub isqrt: usize,
    pub bit: FenwickTree,
}

impl DirichletFenwickHole {
    #[must_use]
    pub fn new(x: usize) -> Self {
        let isqrt = x.isqrt();
        let l = (isqrt << 1) - usize::from(isqrt == x / isqrt);
        Self {
            x,
            isqrt,
            bit: FenwickTree::new(l as usize, 0),
        }
    }
    #[must_use]
    pub fn zeta(x: usize) -> Self {
        let isqrt = x.isqrt();
        let arr = FIArray::keys(x)
            .map(|v| v as usize)
            .collect_vec()
            .into_boxed_slice();
        Self {
            x,
            isqrt,
            bit: FenwickTree::new_with(arr.len(), |i| if i == 0 { 1 } else { arr[i] - arr[i - 1] }),
        }
    }

    #[must_use]
    pub fn get_prefix(&self, v: usize) -> usize {
        self.bit.sum(self.get_index(v))
    }
    #[must_use]
    pub fn get_bucket_prefix(&self, i: usize) -> usize {
        self.bit.sum(i)
    }
    #[must_use]
    pub fn get_index(&self, v: usize) -> usize {
        if v == 0 {
            unsafe { core::hint::unreachable_unchecked() };
        } else if v <= self.isqrt {
            v - 1
        } else {
            let l = self.bit.len;
            l - (self.x / v)
        }
    }
    /// 1 / (1 - w x^-s)
    pub fn sparse_mul_unlimited(&mut self, x: usize, w: usize) {
        let lim = self.x / x;
        let mut prev = 0;
        let mut i = 1;
        while i <= lim / i {
            let cur = self.bit.sum(i - 1);
            if cur != prev {
                self.bit.add(self.get_index(i * x), w * (cur - prev));
                prev = cur;
            }
            i += 1;
        }
        for j in (1..=lim / i).rev() {
            let cur = self.get_prefix(lim / j);
            if cur != prev {
                self.bit.add(self.bit.len - j as usize, w * (cur - prev));
                prev = cur;
            }
        }
    }
    /// 1 - w x^-s
    pub fn sparse_mul_at_most_one(&mut self, x: usize, w: usize) {
        let lim = self.x / x;
        let len = self.bit.len;
        let mut j = 1;
        let mut cur = self.get_prefix(lim);
        while (j + 1) <= lim / (j + 1) {
            let next = self.get_prefix(lim / (j + 1));
            if next != cur {
                self.bit.sub(len - j, w * (cur - next));
                cur = next;
            }
            j += 1;
        }
        for i in (2..=lim / j).rev() {
            let next = self.bit.sum(i - 2);
            if next != cur {
                self.bit.sub(self.get_index(x * i), w * (cur - next));
                cur = next;
            }
        }
        if cur != 0 {
            self.bit.sub(self.get_index(x), w * cur);
        }
    }
}
impl std::convert::From<FIArray> for DirichletFenwickHole {
    fn from(mut value: FIArray) -> Self {
        value.adjacent_difference();
        Self {
            x: value.x,
            isqrt: value.isqrt,
            bit: FenwickTree::new_with(value.arr.len(), |i| value.arr[i]),
        }
    }
}
impl std::convert::From<DirichletFenwickHole> for FIArray {
    fn from(value: DirichletFenwickHole) -> Self {
        Self {
            x: value.x,
            isqrt: value.isqrt,
            arr: value.bit.flatten(),
        }
    }
}
#[must_use]
pub fn lucy_fenwick_hole(x: usize) -> FIArray {
    const LUT: [usize; 30] = [
        0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 8,
    ];
    let mut s = FIArray::new(x);
    let xsqrt = s.isqrt;
    let len = s.arr.len();

    for (i, v) in FIArray::keys(x).enumerate() {
        //s.arr[i] = (v + 1) >> 1;
        s.arr[i] = ((v / 30) << 3) + LUT[v % 30] - 1 + 3;
    }
    s[1] = 0;
    s[2] = 1;
    s[3] = 2;
    s[4] = 2;
    s[5] = 3;

    let mut s = DirichletFenwickHole::from(s);
    let mut sp = s.get_bucket_prefix(7 - 2);
    let mut cnt_query: usize = 1;
    let mut cnt_modify: usize = 0;
    let cutoff = xsqrt
        .isqrt()
        .max(2 * iroot::<3>((xsqrt / x.ilog2() as usize).pow(2)))
        | 1; // iroot::<3>(x) | 1;
    //dbg!(cutoff, iroot::<3>(x) | 1);
    for p in (7..=cutoff).step_by(2) {
        let sp1 = s.get_bucket_prefix(p - 1);
        cnt_query += 1;
        if sp1 == sp {
            continue;
        }

        let lim = x / p;
        let mut j = 1;
        let mut cur = s.get_prefix(lim);
        cnt_query += 1;
        while (j + 1) <= lim / (j + 1) {
            let next = s.get_prefix(lim / (j + 1));
            cnt_query += 1;
            if next != cur {
                s.bit.sub(len - j, cur - next);
                cnt_modify += 1;
                cur = next;
            }
            j += 1;
        }
        for i in (p..=lim / j).rev() {
            let next = s.bit.sum(i - 2);
            cnt_query += 1;
            if next != cur {
                s.bit.sub(s.get_index(p * i), cur - next);
                cnt_modify += 1;
                cur = next;
            }
        }
        sp = sp1;
    }
    dbg!(cnt_modify,cnt_query);
    let mut s = FIArray::from(s);
    for p in (cutoff + 2..=xsqrt).step_by(2) {
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
// 1e17
// res = 2623557157654233, took 805.1480913s
// res = 2623557157654233, took 681.3617369s
// res = 2623557157654233, took 854.8861609s
// res = 2623557157654233, took 742.6827623s
// res = 2623557157654233, took 1295.5591953s
// res = 2623557157654233, took 1212.7197351s
pub fn main() {
    const N: usize = 1e14 as _;
    let start = Instant::now();
    let count = log_zeta_2(N)[N]; // n^(2/3) / \log n
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");   

    let start = Instant::now();
    let count = log_zeta_2_hole(N)[N]; // n^(2/3) / \log n
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    let start = Instant::now();
    let count = log_zeta_3_hole(N)[N]; // n^(2/3) / \log n
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    let start = Instant::now();
    let count = log_zeta_3(N)[N]; // n^(2/3) / \log n
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    let start = Instant::now();
    let count = log_zeta_3_odd(N)[N]; // n^(2/3) / \log n
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    let start = Instant::now();
    let count = log_zeta(N as _)[N as _];
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");
    
    let start = Instant::now();
    let count = log_zeta_hole(N as _)[N as _];
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    
    let start = Instant::now();
    let count = log_zeta_fast(N as _)[N as _];
    let end = start.elapsed();
    println!("n^5/8: res = {count}, took {end:?}");

    
    let start = Instant::now();
    let count = log_zeta_fast_alt(N as _)[N as _];
    let end = start.elapsed();
    println!("n^5/8: res = {count}, took {end:?}");

    
    let start = Instant::now();
    let count = log_zeta_fast_alt_2(N as _)[N as _];
    let end = start.elapsed();
    println!("n^5/8: res = {count}, took {end:?}");

    let start = Instant::now();
    let count = lucy_fenwick(N as _)[N as _];
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");

    let start = Instant::now();
    let count = lucy_fenwick_hole(N as _)[N as _];
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");
    //assert_eq!(log_zeta_fast(N),log_zeta_2_hole(N));
    //assert_eq!(lucy_fenwick(N),lucy_fenwick_hole(N));
}
// TODO: finish
fn log_zeta_fast(n: usize) -> FIArray {
    const INVS: [usize; 4] = [0, 6, 3, 2];
    let start = std::time::Instant::now();
    let mut zeta = DirichletFenwickHole::zeta(n);
    let mut zeta_2 = DirichletFenwickHole::from(divisor_sums(n));
    dbg!(start.elapsed());
    let rt = zeta.isqrt;
    let len = zeta.bit.len;

    let mut ret = FIArray::new(n);

    let x = iroot::<4>(n) + 1;
    // remove contributions of small primes
    for p in 2..x {
        if zeta.bit.sum(p - 1) == 1 {
            //not prime
            continue;
        }
        ret.arr[p - 1] = 1;
        zeta.sparse_mul_at_most_one(p, 1);
        zeta_2.sparse_mul_at_most_one(p, 1);
        zeta_2.sparse_mul_at_most_one(p, 1);
    }
    zeta.bit.dec(0);
    let mut zeta = FIArray::from(zeta);
    let mut zeta_2 = FIArray::from(zeta_2);
    for i in 0..len {
        zeta_2.arr[i] -= 2 * zeta.arr[i] + 1;
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
        let v = ret.arr[x - 1] / 6;
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
        ret.arr[i] /= 6;
        ret.arr[i] += ret.arr[i - 1];
    }
    ret
}

fn log_zeta_fast_alt(n: usize) -> FIArray {
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
    const INVS: [usize; 4] = [0, 6, 3, 2];
    let start = std::time::Instant::now();
    let mut zeta = DirichletFenwickHole::zeta(n);
    let mut zeta_2 = DirichletFenwickHole::from(divisor_sums(n));
    dbg!(start.elapsed());
    let rt = zeta.isqrt;
    let len = zeta.bit.len;

    let mut ret = FIArray::new(n);

    let mut primes = vec![];
    let x = iroot::<4>(n) + 1;
    // remove contributions of small primes
    for p in 2..x {
        if zeta.bit.sum(p - 1) == 1 {
            //not prime
            continue;
        }
        primes.push(p);
        ret.arr[p - 1] = 1;
        zeta.sparse_mul_at_most_one(p, 1);
        zeta_2.sparse_mul_at_most_one(p, 2);
    }
    dbg!(x,primes.len());

    zeta.bit.dec(0);
    let mut zeta = FIArray::from(zeta);
    let mut zeta_2 = FIArray::from(zeta_2);
    dbg!(start.elapsed());

    zeta_2 = mult_correction(&zeta_2, &primes);
    for i in 0..len {
        zeta_2.arr[i] -= 2 * zeta.arr[i] + 1;
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
        let v = ret.arr[x - 1] / 6;
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
        ret.arr[i] /= 6;
        ret.arr[i] += ret.arr[i - 1];
    }
    ret
}

// 1e17: res = 2623557157654233, took 2534.5070687s
// 1e16: res = 279238341033925, took 550.2544082s
pub fn log_zeta_fast_alt_2(n: usize) -> FIArray {
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
    const INVS: [usize; 4] = [0, 6, 3, 2];
    let start = std::time::Instant::now();
    let mut zeta = DirichletFenwickHole::zeta(n);
    let mut zeta_2 = DirichletFenwickHole::from(divisor_sums(n));
    dbg!(start.elapsed());
    let rt = zeta.isqrt;
    let len = zeta.bit.len;

    let mut ret = FIArray::new(n);
    let mut cnt_query = 0usize;
    let mut cnt_modify = 0usize;
    let mut cnt_query_2 = 0usize;
    let mut cnt_modify_2 = 0usize;
    let mut primes = vec![];
    let x = iroot::<4>(n) + 1;
    // remove contributions of small primes
    for p in 2..x {
        cnt_query += 1;
        if zeta.bit.sum(p - 1) == 1 {
            //not prime
            continue;
        }
        primes.push(p);
        ret.arr[p - 1] = 1;

        //zeta.sparse_mul_at_most_one(p, 1);
        //zeta_2.sparse_mul_at_most_one(p, 2);
        // instead of executing zeta_2.sparse_mul_at_most_one(p, 1) twice to remove p from zeta_2,
        // execute zeta_2.sparse_mul_at_most_one(p, 2) and correct later
        {
            let lim = n / p;
            let mut j = 1;
            let ind = zeta.get_index(lim);
            let mut cur_zeta = zeta.get_bucket_prefix(ind);
            let mut cur_zeta_2 = zeta_2.get_bucket_prefix(ind);
            cnt_query += 1;
            cnt_query += 1;

            while (j + 1) <= lim / (j + 1) {
                let ind = zeta.get_index(lim / (j + 1));
                let next_zeta = zeta.get_bucket_prefix(ind);
                let next_zeta_2 = zeta_2.get_bucket_prefix(ind);
                cnt_query += 1;
                cnt_query_2 += 1;
                if next_zeta != cur_zeta {
                    zeta.bit.sub(len - j, cur_zeta - next_zeta);
                    cnt_modify += 1;
                    cur_zeta = next_zeta;
                }
                if next_zeta_2 != cur_zeta_2 {
                    zeta_2.bit.sub(len - j, 2 * (cur_zeta_2 - next_zeta_2));
                    cnt_modify_2 += 1;
                    cur_zeta_2 = next_zeta_2;
                }
                j += 1;
            }
            for i in (2..=lim / j).rev() {
                let next_zeta = zeta.get_bucket_prefix(i - 2);
                let next_zeta_2 = zeta_2.get_bucket_prefix(i - 2);
                cnt_query += 1;
                cnt_query_2 += 1;
                let ind = zeta.get_index(p * i);
                if next_zeta != cur_zeta {
                    zeta.bit.sub(ind, cur_zeta - next_zeta);
                    cnt_modify += 1;
                    cur_zeta = next_zeta;
                }
                if next_zeta_2 != cur_zeta_2 {
                    zeta_2.bit.sub(ind, 2 * (cur_zeta_2 - next_zeta_2));
                    cnt_modify_2 += 1;
                    cur_zeta_2 = next_zeta_2;
                }
            }
            if cur_zeta != 0 {
                zeta.bit.sub(p - 1, cur_zeta);
                cnt_modify += 1;
            }
            if cur_zeta_2 != 0 {
                zeta_2.bit.sub(p - 1, 2 * cur_zeta_2);
                cnt_modify_2 += 1;
            }
        }
    }
    zeta.bit.dec(0);
    cnt_modify += 1;
    dbg!(cnt_modify,cnt_query);
    dbg!(cnt_modify_2,cnt_query_2);

    let mut zeta = FIArray::from(zeta);
    let mut zeta_2 = FIArray::from(zeta_2);

    zeta_2 = mult_correction(&zeta_2, &primes);
    for i in 0..len {
        zeta_2.arr[i] -= 2 * zeta.arr[i] + 1;
    }

    // zeta now equals zeta_t - 1, and zeta_2 (zeta_t - 1)^2
    // compute log(zeta_t) using log(x + 1) = x^3 / 3 - x^2 / 2 + x
    // in order to not have to deal with rational numbers, we compute 6 * log(zeta_t)
    // and adjust later

    for i in x..=len {
        ret.arr[i - 1] = zeta.arr[i - 1] * INVS[1];
    }

    for i in rt + 1..=len {
        ret.arr[i - 1] -= zeta_2.arr[i - 1] * INVS[2];
    }

    let zeta_3 = {
        zeta.arr[rt..].fill(0);
        for i in x..=rt {
            if zeta.arr[i - 1] == zeta.arr[i - 2] {
                continue;
            }
            for j in 1..=rt / i {
                zeta.arr[len - j] += zeta_2.arr[len - i * j];
            }
        }
        //zeta.arr[..rt].fill(0);
        zeta
    };

    let ind = zeta_3.get_index(x.pow(3));
    for i in ind + 1..=len {
        ret.arr[i - 1] += zeta_3.arr[i - 1] * INVS[3];
    }

    // correction phase: get rid of contributions of prime powers
    for i in (x..len).rev() {
        ret.arr[i] -= ret.arr[i - 1];
    }

    for x in x..=rt {
        let v = ret.arr[x - 1] / 6;
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
        ret.arr[i] /= 6;
        ret.arr[i] += ret.arr[i - 1];
    }
    ret
}

#[must_use]
pub fn log_zeta_2_hole(n: usize) -> FIArray {
    const fn inv_odd(mut k: usize) -> usize{
        let mut exp = (1u64 << 63) - 1;
        
        let mut r: usize = 1;
        while exp > 1 {
            r = r.overflowing_mul(k).0;
            k = k.overflowing_mul(k).0;
            exp >>= 1;
        }
        r.overflowing_mul(k).0    
    }
    const INVS: [usize; 5] = [0, 4, 2, 4usize.overflowing_mul(inv_odd(3)).0, 1];
    let start = std::time::Instant::now();
    let mut zeta = DirichletFenwickHole::zeta(n);
    let rt = zeta.isqrt;
    let len = zeta.bit.len;

    let mut ret = FIArray::new(n);

    let x = iroot::<5>(n) + 1;
    // remove contributions of small primes
    for p in 2..x {
        if zeta.bit.sum(p - 1) == 1 {
            //not prime
            continue;
        }
        ret.arr[p - 1] = 1;
        zeta.sparse_mul_at_most_one(p, 1);
    }
    zeta.bit.dec(0);
    let zeta = FIArray::from(zeta);
    println!("Finished sieving: {:?}",start.elapsed());
    // zeta now equals zeta_t - 1
    // compute log(zeta_t) using log(x + 1) = x - x^2 / 2 + x^3 / 3 - x^4 / 4
    // in order to not have to deal with rational numbers, we compute 12 * log(zeta_t)
    // and adjust later

    for i in x..=len {
        ret.arr[i - 1] = zeta.arr[i - 1] * INVS[1];
    }
    println!("started first mul: {:?}",start.elapsed());
    let zeta_2 = mult(&zeta, &zeta);
    println!("finished first mul: {:?}",start.elapsed());

    let mut x_pow = x * x;

    for i in ret.get_index(x_pow)..=len {
        ret.arr[i - 1] -= zeta_2.arr[i - 1] * INVS[2];
    }    
    println!("started second mul: {:?}",start.elapsed());
    let zeta_3 = mult(&zeta, &zeta_2);
    println!("finished second mul: {:?}",start.elapsed());

    x_pow *= x;
    for i in ret.get_index(x_pow)..=len {
        ret.arr[i - 1] += zeta_3.arr[i - 1] * INVS[3];
    }
    println!("started third mul: {:?}",start.elapsed());
    let zeta_4 = mult_sparse(&zeta, &zeta_3);
    println!("finished third mul: {:?}",start.elapsed());
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
pub fn log_zeta_3_hole(n: usize) -> FIArray {
    const INVS: [usize; 7] = [0, 60, 30, 20, 15, 12, 10];
    let start = std::time::Instant::now();
    let mut zeta = DirichletFenwickHole::zeta(n);
    let rt = zeta.isqrt;
    let len = zeta.bit.len;

    let mut ret = FIArray::new(n);

    let x = iroot::<7>(n) + 1;
    // remove contributions of small primes
    for p in 2..x {
        if zeta.bit.sum(p - 1) == 1 {
            //not prime
            continue;
        }
        ret.arr[p - 1] = 1;
        zeta.sparse_mul_at_most_one(p, 1);
    }
    zeta.bit.dec(0);
    let zeta = FIArray::from(zeta);
    println!("Finished sieving: {:?}",start.elapsed());
    // zeta now equals zeta_t - 1
    // compute log(zeta_t) using log(x + 1) = x - x^2 / 2 + x^3 / 3 - x^4 / 4 + x^5 / 5 - x^6 / 6
    // in order to not have to deal with rational numbers, we compute 12 * log(zeta_t)
    // and adjust later

    for i in x..=len {
        ret.arr[i - 1] = zeta.arr[i - 1] * INVS[1];
    }
    println!("started first mul: {:?}",start.elapsed());
    let zeta_2 = mult(&zeta, &zeta);
    println!("finished first mul: {:?}",start.elapsed());

    let mut x_pow = x * x;

    for i in ret.get_index(x_pow)..=len {
        ret.arr[i - 1] -= zeta_2.arr[i - 1] * INVS[2];
    }    
    println!("started second mul: {:?}",start.elapsed());
    let zeta_3 = mult(&zeta, &zeta_2);
    println!("finished second mul: {:?}",start.elapsed());

    x_pow *= x;
    for i in ret.get_index(x_pow)..=len {
        ret.arr[i - 1] += zeta_3.arr[i - 1] * INVS[3];
    }
    println!("started third mul: {:?}",start.elapsed());
    let zeta_4 = mult(&zeta_2, &zeta_2);
    println!("finished third mul: {:?}",start.elapsed());
    x_pow *= x;
    for i in ret.get_index(x_pow)..=len {
        ret.arr[i - 1] -= zeta_4.arr[i - 1] * INVS[4];
    }
    println!("started fourth mul: {:?}",start.elapsed());
    let zeta_5 = mult_sparse(&zeta, &zeta_4);
    println!("finished fourth mul: {:?}",start.elapsed());
    x_pow *= x;
    for i in ret.get_index(x_pow)..=len {
        ret.arr[i - 1] += zeta_5.arr[i - 1] * INVS[5];
    }
    println!("started fifth mul: {:?}",start.elapsed());
    let zeta_6 = mult_sparse(&zeta_2, &zeta_4);
    println!("finished fifth mul: {:?}",start.elapsed());
    x_pow *= x;
    for i in ret.get_index(x_pow)..=len {
        ret.arr[i - 1] -= zeta_6.arr[i - 1] * INVS[6];
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
pub fn log_zeta_3(n: usize) -> FIArray {
    const INVS: [usize; 7] = [0, 60, 30, 20, 15, 12, 10];
    //let start = std::time::Instant::now();
    let mut zeta = FIArray::unit(n);
    let rt = zeta.isqrt;
    let len = zeta.arr.len();

    let mut ret = FIArray::new(n);

    let x = iroot::<7>(n) + 1;
    // remove contributions of small primes
    let keys = zeta.clone();
    for p in 2..x {
        if zeta.arr[p - 1] == 1 {
            //not prime
            continue;
        }
        ret.arr[p - 1] = 1;
        for (i,&v) in keys.arr.iter().enumerate().rev() {
            if v < p {
                break;
            }
            zeta.arr[i] -= zeta[v / p];
        }
    }
    for e in &mut zeta.arr {
        *e -= 1;
    }
    let zeta = zeta;
    //println!("Finished sieving: {:?}",start.elapsed());
    // zeta now equals zeta_t - 1
    // compute log(zeta_t) using log(x + 1) = x - x^2 / 2 + x^3 / 3 - x^4 / 4 + x^5 / 5 - x^6 / 6
    // in order to not have to deal with rational numbers, we compute 12 * log(zeta_t)
    // and adjust later

    for i in x..=len {
        ret.arr[i - 1] = zeta.arr[i - 1] * INVS[1];
    }
    //println!("started first mul: {:?}",start.elapsed());
    let zeta_2 = mult(&zeta, &zeta);
    //println!("finished first mul: {:?}",start.elapsed());

    let mut x_pow = x * x;

    for i in ret.get_index(x_pow)..=len {
        ret.arr[i - 1] -= zeta_2.arr[i - 1] * INVS[2];
    }    
    //println!("started second mul: {:?}",start.elapsed());
    let zeta_3 = mult(&zeta, &zeta_2);
    //println!("finished second mul: {:?}",start.elapsed());

    x_pow *= x;
    for i in ret.get_index(x_pow)..=len {
        ret.arr[i - 1] += zeta_3.arr[i - 1] * INVS[3];
    }
    //println!("started third mul: {:?}",start.elapsed());
    let zeta_4 = mult(&zeta_2, &zeta_2);
    //println!("finished third mul: {:?}",start.elapsed());
    x_pow *= x;
    for i in ret.get_index(x_pow)..=len {
        ret.arr[i - 1] -= zeta_4.arr[i - 1] * INVS[4];
    }
    //println!("started fourth mul: {:?}",start.elapsed());
    let zeta_5 = mult_sparse(&zeta, &zeta_4);
    //println!("finished fourth mul: {:?}",start.elapsed());
    x_pow *= x;
    for i in ret.get_index(x_pow)..=len {
        ret.arr[i - 1] += zeta_5.arr[i - 1] * INVS[5];
    }
    //println!("started fifth mul: {:?}",start.elapsed());
    let zeta_6 = mult_sparse(&zeta_2, &zeta_4);
    //println!("finished fifth mul: {:?}",start.elapsed());
    x_pow *= x;
    for i in ret.get_index(x_pow)..=len {
        ret.arr[i - 1] -= zeta_6.arr[i - 1] * INVS[6];
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
pub fn log_zeta_3_odd(n: usize) -> FIArray {
    const fn inv_odd(mut k: usize) -> usize{
        let mut exp = (1u64 << 63) - 1;
        
        let mut r: usize = 1;
        while exp > 1 {
            r = r.overflowing_mul(k).0;
            k = k.overflowing_mul(k).0;
            exp >>= 1;
        }
        r.overflowing_mul(k).0    
    }
    const INVS: [usize; 7] = [0, 4, 2, inv_odd(3) << 2, 1, inv_odd(5) << 2, inv_odd(3) << 1];
    //let start = std::time::Instant::now();
    let mut zeta = FIArray::unit(n);
    let rt = zeta.isqrt;
    let len = zeta.arr.len();

    let mut ret = FIArray::new(n);

    let x = iroot::<7>(n) + 1;
    // remove contributions of small primes
    let keys = zeta.clone();
    for p in 2..x {
        if zeta.arr[p - 1] == 1 {
            //not prime
            continue;
        }
        ret.arr[p - 1] = 1;
        for (i,&v) in keys.arr.iter().enumerate().rev() {
            if v < p {
                break;
            }
            zeta.arr[i] -= zeta[v / p];
        }
    }
    for e in &mut zeta.arr {
        *e -= 1;
    }
    let zeta = zeta;
    //println!("Finished sieving: {:?}",start.elapsed());
    // zeta now equals zeta_t - 1
    // compute log(zeta_t) using log(x + 1) = x - x^2 / 2 + x^3 / 3 - x^4 / 4 + x^5 / 5 - x^6 / 6
    // in order to not have to deal with rational numbers, we compute 4 * log(zeta_t)
    // and adjust later

    for i in x..=len {
        ret.arr[i - 1] = zeta.arr[i - 1] * INVS[1];
    }
    //println!("started first mul: {:?}",start.elapsed());
    let zeta_2 = mult(&zeta, &zeta);
    //println!("finished first mul: {:?}",start.elapsed());

    let mut x_pow = x * x;

    for i in ret.get_index(x_pow)..=len {
        ret.arr[i - 1] -= zeta_2.arr[i - 1] * INVS[2];
    }    
    //println!("started second mul: {:?}",start.elapsed());
    let zeta_3 = mult(&zeta, &zeta_2);
    //println!("finished second mul: {:?}",start.elapsed());

    x_pow *= x;
    for i in ret.get_index(x_pow)..=len {
        ret.arr[i - 1] += zeta_3.arr[i - 1] * INVS[3];
    }
    //println!("started third mul: {:?}",start.elapsed());
    let zeta_4 = mult(&zeta_2, &zeta_2);
    //println!("finished third mul: {:?}",start.elapsed());
    x_pow *= x;
    for i in ret.get_index(x_pow)..=len {
        ret.arr[i - 1] -= zeta_4.arr[i - 1] * INVS[4];
    }
    //println!("started fourth mul: {:?}",start.elapsed());
    let zeta_5 = mult_sparse(&zeta, &zeta_4);
    //println!("finished fourth mul: {:?}",start.elapsed());
    x_pow *= x;
    for i in ret.get_index(x_pow)..=len {
        ret.arr[i - 1] += zeta_5.arr[i - 1] * INVS[5];
    }
    //println!("started fifth mul: {:?}",start.elapsed());
    let zeta_6 = mult_sparse(&zeta_2, &zeta_4);
    //println!("finished fifth mul: {:?}",start.elapsed());
    x_pow *= x;
    for i in ret.get_index(x_pow)..=len {
        ret.arr[i - 1] -= zeta_6.arr[i - 1] * INVS[6];
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
pub fn log_zeta_hole(n: usize) -> FIArray {
    const fn inv_odd(mut k: usize) -> usize{
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
    let mut zeta = DirichletFenwickHole::zeta(n);
    let rt = zeta.isqrt;
    let len = zeta.bit.len;

    let mut ret = FIArray::new(n);

    let x = iroot::<4>(n) + 1;
    // remove contributions of small primes
    for p in 2..x {
        if zeta.bit.sum(p - 1) == 1 {
            //not prime
            continue;
        }
        ret.arr[p - 1] = 1;
        zeta.sparse_mul_at_most_one(p, 1);
    }
    zeta.bit.dec(0);
    let mut zeta = FIArray::from(zeta);
    dbg!(start.elapsed());

    // zeta now equals zeta_t - 1
    // compute log(zeta_t) using log(x + 1) = x^3 / 3 - x^2 / 2 + x
    // in order to not have to deal with rational numbers, we compute 2 * log(zeta_t)
    // and adjust later

    for i in x..=len {
        ret.arr[i - 1] = zeta.arr[i - 1] * INVS[1];
    }
    dbg!(start.elapsed());

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
        for i in 1..len {
            res.arr[i] += res.arr[i - 1];
        }
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
    dbg!(start.elapsed());
    /* let ind = pow_zeta.get_index(x.pow(2));
    assert!(pow_zeta.arr[..ind].iter().all(|&e| e == 0));
    dbg!(ind, x.pow(2), rt - 1, len, len - ind); */

    for i in rt + 1..=len {
        ret.arr[i - 1] -= pow_zeta.arr[i - 1] * INVS[2];
    }
    dbg!(start.elapsed());
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
    dbg!(start.elapsed());

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
