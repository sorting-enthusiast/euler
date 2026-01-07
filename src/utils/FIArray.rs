use super::multiplicative_function_summation::{
    sum_n_i32, sum_n_i64, sum_n_i128, sum_n_isize, sum_n_u32, sum_n_u64, sum_n_u128, sum_n_usize,
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
        }
    )+ };
}
FIArray_impl_for!(u32, i32, u64, i64, usize, isize, u128, i128);
pub type FIArray = FIArrayUsize;
