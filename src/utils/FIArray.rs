use super::{
    fenwick::{
        FenwickTreeI32, FenwickTreeI64, FenwickTreeI128, FenwickTreeIsize, FenwickTreeU32,
        FenwickTreeU64, FenwickTreeU128, FenwickTreeUsize,
    },
    multiplicative_function_summation::{
        sum_n_i32, sum_n_i64, sum_n_i128, sum_n_isize, sum_n_u32, sum_n_u64, sum_n_u128,
        sum_n_usize,
    },
};
use itertools::Itertools;
use paste::paste;
use std::ops::{Index, IndexMut};
macro_rules! FIArray_impl_for {
    ($($type:ty),+) => { $(
        paste!{
            #[derive(Debug, Clone, PartialEq, Eq)]
            pub struct [<FIArray $type:camel>] {
                pub x: $type,
                pub isqrt: $type,
                pub arr: Box<[$type]>,
            }

            impl [<FIArray $type:camel>] {
                #[must_use] pub fn new(x: $type) -> Self {
                    let isqrt = x.isqrt();
                    let l = (isqrt << 1) - $type::from(isqrt == x / isqrt);
                    Self {
                        x,
                        isqrt,
                        arr: vec![0; l as usize].into_boxed_slice(),
                    }
                }
                #[must_use] pub fn unit(x: $type) -> Self {
                    let isqrt = x.isqrt();
                    let arr = Self::keys(x).collect_vec().into_boxed_slice();
                    Self { x, isqrt, arr }
                }
                #[must_use] pub fn id<const MOD: $type>(x: $type) -> Self {
                    let isqrt = x.isqrt();
                    let arr = Self::keys(x).map(|v| [<sum_n_ $type>]::<MOD>(v)).collect_vec().into_boxed_slice();
                    Self { x, isqrt, arr }
                }
                #[must_use] pub fn eps(x: $type) -> Self {
                    let isqrt = x.isqrt();
                    let l = (isqrt << 1) - $type::from(isqrt == x / isqrt);
                    Self {
                        x,
                        isqrt,
                        arr: vec![1; l as usize].into_boxed_slice(),
                    }
                }
                #[must_use] pub fn keys(x: $type) -> impl DoubleEndedIterator<Item = $type> + use<> {
                    let isqrt = x.isqrt();
                    (1..=isqrt)
                        .chain((isqrt != x / isqrt).then_some(x / isqrt))
                        .chain((1..isqrt).rev().map(move |n| x / n))
                }
                #[must_use] pub fn large_keys(x: $type) -> impl DoubleEndedIterator<Item = $type> + use<> {
                    let isqrt = x.isqrt();
                    (1..isqrt).chain((isqrt != x / isqrt).then_some(isqrt))
                }
                #[must_use] pub fn get_index(&self, v: $type) -> usize {
                    if v <= 0 {
                        0
                    } else if v <= self.isqrt {
                        v as usize - 1
                    } else {
                        let l = self.arr.len();
                        l - (self.x / v) as usize
                    }
                }
            }

            impl Index<$type> for [<FIArray $type:camel>] {
                type Output = $type;
                fn index(&self, v: $type) -> &Self::Output {
                    if v <= 0 {
                        &0
                    } else if v <= self.isqrt {
                        &self.arr[v as usize - 1]
                    } else {
                        let l = self.arr.len();
                        &self.arr[l - (self.x / v) as usize]
                    }
                }
            }
            impl IndexMut<$type> for [<FIArray $type:camel>] {
                fn index_mut(&mut self, v: $type) -> &mut Self::Output {
                    if v <= 0 {
                        unsafe { core::hint::unreachable_unchecked() };
                    }
                    if v <= self.isqrt {
                        &mut self.arr[v as usize - 1]
                    } else {
                        let l = self.arr.len();
                        &mut self.arr[l - (self.x / v) as usize]
                    }
                }
            }
            #[derive(Debug, Clone, PartialEq, Eq)]
            pub struct [<DirichletFenwick $type:camel>] {
                pub x: $type,
                pub isqrt: $type,
                pub bit: [<FenwickTree $type:camel>],
            }

            impl [<DirichletFenwick $type:camel>] {
                #[must_use] pub fn new(x: $type) -> Self {
                    let isqrt = x.isqrt();
                    let l = (isqrt << 1) - $type::from(isqrt == x / isqrt);
                    Self {
                        x,
                        isqrt,
                        bit: [<FenwickTree $type:camel>]::new(l as usize, 0),
                    }
                }
                #[must_use] pub fn zeta(x: $type) -> Self {
                    let isqrt = x.isqrt();
                    let mut arr = Self::keys(x).collect_vec().into_boxed_slice();
                    for i in (2..arr.len()).rev() {
                        let j = i & (i + 1);
                        if j != 0 {
                            arr[i] -= arr[j - 1];
                        }
                    }
                    Self { x, isqrt, bit: [<FenwickTree $type:camel>] { 0: arr } }
                }
                #[must_use] pub fn eps(x: $type) -> Self {
                    let isqrt = x.isqrt();
                    let l = (isqrt << 1) - $type::from(isqrt == x / isqrt);
                    let mut arr = vec![0; l as usize].into_boxed_slice();
                    arr[0] = 1;
                    let mut f = [<FenwickTree $type:camel>] { 0: arr };
                    f.construct();
                    Self {
                        x,
                        isqrt,
                        bit: f,
                    }
                }
                #[must_use] pub fn keys(x: $type) -> impl DoubleEndedIterator<Item = $type> + use<> {
                    let isqrt = x.isqrt();
                    (1..=isqrt)
                        .chain((isqrt != x / isqrt).then_some(x / isqrt))
                        .chain((1..isqrt).rev().map(move |n| x / n))
                }
                #[must_use] pub fn get_prefix(&self, v: $type) -> $type {
                    self.bit.sum(self.get_index(v))
                }
                #[must_use] pub fn get_bucket_prefix(&self, i: usize) -> $type {
                    self.bit.sum(i)
                }
                #[must_use] pub fn get_index(&self, v: $type) -> usize {
                    if v <= 0 {
                        unsafe { core::hint::unreachable_unchecked() };
                    } else if v <= self.isqrt {
                        v as usize - 1
                    } else {
                        let l = self.bit.0.len();
                        l - (self.x / v) as usize
                    }
                }
                /// 1 / (1 - w x^-s)
                pub fn sparse_mul_unlimited(&mut self, x: $type, w: $type) {
                    let lim = self.x / x;
                    let mut prev = 0;
                    let mut i = 1;
                    while i <= lim / i {
                        let cur = self.bit.sum(i as usize - 1);
                        if cur != prev {
                            self.bit.add(self.get_index(i * x), w * (cur - prev));
                            prev = cur;
                        }
                        i += 1;
                    }
                    for j in (1..=lim / i).rev() {
                        let cur = self.get_prefix(lim / j);
                        if cur != prev {
                            self.bit.add(self.bit.0.len() - j as usize, w * (cur - prev));
                            prev = cur;
                        }
                    }
                }
                /// 1 - w x^-s
                pub fn sparse_mul_at_most_one(&mut self, x: $type, w: $type) {
                    let lim = self.x / x;
                    let len = self.bit.0.len();
                    let mut j = 1;
                    let mut cur = self.get_prefix(lim);
                    while (j + 1) <= lim / (j + 1) {
                        let next = self.get_prefix(lim / (j + 1));
                        if next != cur {
                            self.bit.sub(len - j as usize, w * (cur - next));
                            cur = next;
                        }
                        j += 1;
                    }
                    for i in (2..=lim / j).rev() {
                        let next = self.bit.sum(i as usize - 2);
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
            impl std::convert::From<[<FIArray $type:camel>]> for [<DirichletFenwick $type:camel>] {
                fn from(value: [<FIArray $type:camel>]) -> Self {
                    let mut bit = value.arr;
                    for i in (2..bit.len()).rev() {
                        let j = i & (i + 1);
                        if j != 0 {
                            bit[i] -= bit[j - 1];
                        }
                    }
                    Self {
                        x: value.x,
                        isqrt: value.isqrt,
                        bit: [<FenwickTree $type:camel>] { 0: bit}
                    }
                }
            }
            impl std::convert::From<[<DirichletFenwick $type:camel>]> for [<FIArray $type:camel>] {
                fn from(value: [<DirichletFenwick $type:camel>]) -> Self {
                    Self {
                        x: value.x,
                        isqrt: value.isqrt,
                        arr: value.bit.flatten()
                    }
                }
            }
        }
    )+ };
}
FIArray_impl_for!(u32, i32, u64, i64, usize, isize, u128, i128);
pub type FIArray = FIArrayUsize;
pub type DirichletFenwick = DirichletFenwickUsize;
