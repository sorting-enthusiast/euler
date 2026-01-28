use itertools::Itertools;
use paste::paste;
macro_rules! FenwickTree_impl_for {
    ($($type:ty),+) => { $(
        paste!{
            #[derive(Debug, Clone, PartialEq, Eq)]
            pub struct [<FenwickTree $type:camel>](pub Box<[$type]>);
            impl [<FenwickTree $type:camel>] {
                #[must_use]
                pub fn new(len: usize, init: $type) -> Self {
                    let mut ret = Self(vec![init; len].into_boxed_slice());
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
                    let mut ret = Self((0..len).map(init).collect_vec().into_boxed_slice());
                    ret.construct();
                    ret
                }
                #[must_use]
                pub fn sum(&self, i: usize) -> $type {
                    let mut i = i + 1;
                    let mut sum = 0;
                    while i > 0 {
                        sum += self.0[i - 1];
                        i &= i - 1;
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
                #[must_use]
                pub fn flatten(self) -> Box<[$type]> {
                    let mut ret = self.0;
                    for i in 2..ret.len() {
                        let j = i & (i + 1);
                        if j != 0 {
                            ret[i] += ret[j - 1];
                        }
                    }
                    ret
                }
            }
        #[derive(Debug, Clone, PartialEq, Eq)]
            pub struct [<FenwickTree $type:camel Mod>]<const MOD: $type>(pub Box<[$type]>);
            impl<const MOD: $type> [<FenwickTree $type:camel Mod>]<MOD> {
                #[must_use]
                pub fn new(len: usize, init: $type) -> Self {
                    let mut ret = Self(vec![init % MOD; len].into_boxed_slice());
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
                            if self.0[r - 1] >= MOD {
                                self.0[r - 1] -= MOD;
                            }
                        }
                    }
                }
                #[must_use]
                pub fn new_with(len: usize, init: impl FnMut(usize) -> $type) -> Self {
                    let mut ret = Self((0..len).map(init).map(|v| v % MOD).collect_vec().into_boxed_slice());
                    ret.construct();
                    ret
                }
                #[must_use]
                pub fn sum(&self, i: usize) -> $type {
                    let mut i = i + 1;
                    let mut sum = 0;
                    while i > 0 {
                        sum += self.0[i - 1];
                        if sum >= MOD {
                            sum -= MOD;
                        }
                        i &= i - 1;
                    }
                    sum
                }
                pub fn add(&mut self, mut i: usize, v: $type) {
                    let v = v % MOD;
                    while i < self.0.len() {
                        self.0[i] += v;
                        if self.0[i] >= MOD {
                            self.0[i] -= MOD;
                        }
                        i |= i + 1;
                    }
                }
                pub fn sub(&mut self, mut i: usize, v: $type) {
                    while i < self.0.len() {
                        self.0[i] += MOD - v;
                        if self.0[i] >= MOD {
                            self.0[i] -= MOD;
                        }
                        i |= i + 1;
                    }
                }

                pub fn inc(&mut self, mut i: usize) {
                    while i < self.0.len() {
                        self.0[i] += 1;
                        if self.0[i] >= MOD {
                            self.0[i] -= MOD;
                        }
                        i |= i + 1;
                    }
                }
                pub fn dec(&mut self, mut i: usize) {
                    while i < self.0.len() {
                        self.0[i] += MOD - 1;
                        if self.0[i] >= MOD {
                            self.0[i] -= MOD;
                        }
                        i |= i + 1;
                    }
                }
                #[must_use]
                pub fn flatten(self) -> Box<[$type]> {
                    let mut ret = self.0;
                    for i in 2..ret.len() {
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
FenwickTree_impl_for!(u32, i32, u64, i64, usize, isize, u128, i128);
pub type FenwickTree = FenwickTreeUsize;
