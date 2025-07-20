pub struct FenwickTree(Vec<i64>);

impl FenwickTree {
    pub fn new(len: usize, init: i64) -> Self {
        let mut v = vec![init; len];
        for i in 1..len {
            let r = i + (i & (!i + 1));
            if r <= len {
                v[r - 1] += v[i - 1];
            }
        }
        Self(v)
    }
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
}
