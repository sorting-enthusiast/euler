use paste::paste;
use std::ops::{Index, IndexMut};
macro_rules! FIArray_impl_for {
    ($($type:ty),+) => { $(
        paste!{
            #[derive(Debug, Clone, PartialEq, Eq)]
            pub struct [<FIArray $type:camel>] {
                x: $type,
                isqrt: $type,
                pub arr: Vec<$type>,
            }

            impl [<FIArray $type:camel>] {
                pub fn new(x: $type) -> Self {
                    let isqrt = x.isqrt();
                    let l = (isqrt << 1) - $type::from(isqrt == x / isqrt);
                    Self {
                        x,
                        isqrt,
                        arr: vec![0; l as usize],
                    }
                }
                pub fn at(&self, i: usize) -> $type {
                    //assert!(i < self.arr.len());
                    unsafe { core::hint::assert_unchecked(i < self.arr.len()) };

                    self.arr[i]
                }
                pub fn set(&mut self, i: usize, v: $type) {
                    //assert!(i < self.arr.len());
                    unsafe { core::hint::assert_unchecked(i < self.arr.len()) };

                    self.arr[i] = v;
                }
                pub fn keys(x: $type) -> impl DoubleEndedIterator<Item = $type> + use<> {
                    let isqrt = x.isqrt();
                    (1..=isqrt)
                        .chain((isqrt != x / isqrt).then_some(x / isqrt))
                        .chain((1..isqrt).rev().map(move |n| x / n))
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
