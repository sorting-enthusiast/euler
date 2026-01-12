use itertools::Itertools;
#[derive(Clone)]
pub struct FenwickTree(pub Box<[i64]>);

impl FenwickTree {
    #[must_use]
    pub fn new(len: usize, init: i64) -> Self {
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
    pub fn new_with(len: usize, init: impl FnMut(usize) -> i64) -> Self {
        let mut ret = Self((0..len).map(init).collect_vec().into_boxed_slice());
        ret.construct();
        ret
    }
    #[must_use]
    pub fn sum(&self, i: usize) -> i64 {
        let mut i = i + 1;
        let mut sum = 0;
        while i != 0 {
            sum += self.0[i - 1];
            i &= i - 1;
        }
        sum
    }
    pub fn add(&mut self, i: usize, v: i64) {
        let mut k = i + 1;
        while k <= self.0.len() {
            self.0[k - 1] += v;
            k += k & (!k + 1);
        }
    }
    #[must_use]
    pub fn flatten(self) -> Box<[i64]> {
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
#[derive(Clone)]
pub struct FenwickTreeUsize(pub Box<[usize]>);

impl FenwickTreeUsize {
    #[must_use]
    pub fn new(len: usize, init: usize) -> Self {
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
    pub fn new_with(len: usize, init: impl FnMut(usize) -> usize) -> Self {
        let mut ret = Self((0..len).map(init).collect_vec().into_boxed_slice());
        ret.construct();
        ret
    }
    #[must_use]
    pub fn sum(&self, i: usize) -> usize {
        let mut i = i + 1;
        let mut sum = 0;
        while i > 0 {
            sum += self.0[i - 1];
            i &= i - 1;
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
    #[must_use]
    pub fn flatten(self) -> Box<[usize]> {
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

#[derive(Clone)]
pub struct FenwickTreeU32(pub Vec<u32>);

impl FenwickTreeU32 {
    #[must_use]
    pub fn new(len: usize, init: u32) -> Self {
        let mut v = vec![init; len];
        if init != 0 {
            for i in 1..len {
                let r = i + (i & (!i + 1));
                if r <= len {
                    v[r - 1] += v[i - 1];
                }
            }
        }
        Self(v)
    }
    pub fn new_with(len: usize, init: impl FnMut(usize) -> u32) -> Self {
        let mut v = (0..len).map(init).collect_vec();
        for i in 1..len {
            let r = i + (i & (!i + 1));
            if r <= len {
                v[r - 1] += v[i - 1];
            }
        }
        Self(v)
    }
    #[must_use]
    pub fn sum(&self, i: usize) -> u32 {
        let mut i = i + 1;
        let mut sum = 0;
        while i != 0 {
            sum += self.0[i - 1];
            i &= i - 1;
        }
        sum
    }
    pub fn add(&mut self, mut i: usize, v: u32) {
        while i < self.0.len() {
            self.0[i] += v;
            i |= i + 1;
        }
    }
    pub fn sub(&mut self, mut i: usize, v: u32) {
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
    pub fn flatten(self) -> Vec<u32> {
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
