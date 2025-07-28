use crate::utils::prime_sieves::sift;

//Yields (n, h(n) mod m) where n are the O(sqrt x) powerful numbers
//up to x, and h is any multiplicative function.
pub struct PowerfulExt<H, const M: i64>
where
    H: Fn(i64, i64) -> i64,
{
    x: i64,
    h: H,
    ps: Vec<u64>,
    stack: Vec<(i64, i64, usize)>,
    pub max_len: usize,
}
impl<H, const M: i64> PowerfulExt<H, M>
where
    H: Fn(i64, i64) -> i64,
{
    pub fn new(x: i64, h: H) -> Self {
        let stack = vec![(1, 1, 0)];
        let ps = sift(x.isqrt() as u64 + 1);
        Self {
            x,
            h,
            ps,
            stack,
            max_len: 1,
        }
    }
}
impl<H, const M: i64> Iterator for PowerfulExt<H, M>
where
    H: Fn(i64, i64) -> i64,
{
    type Item = (i64, i64);
    fn next(&mut self) -> Option<Self::Item> {
        while let Some((n, hn, i)) = self.stack.pop() {
            if let Some(&p) = self.ps.get(i) {
                let xdivn = (self.x / n) as u64;
                if p * p > xdivn {
                    return Some((n, hn));
                }
                self.stack.push((n, hn, i + 1));
                let mut pp = p * p;
                let mut e = 2;
                while pp <= xdivn {
                    self.stack
                        .push((n * pp as i64, (hn * (self.h)(p as i64, e)) % M, i + 1));
                    if pp > xdivn / p {
                        break;
                    }
                    pp *= p;
                    e += 1;
                }
                if self.stack.len() > self.max_len {
                    self.max_len = self.stack.len();
                }
            } else {
                return Some((n, hn));
            }
        }
        None
    }
}

pub struct PowerfulExtAlt<'a, H, const M: i64>
where
    H: Fn(i64, i64) -> i64,
{
    x: i64,
    h: H,
    ps: &'a [u64],
    stack: Vec<(i64, i64, usize)>,
    pub max_len: usize,
}
impl<'a, H, const M: i64> PowerfulExtAlt<'a, H, M>
where
    H: Fn(i64, i64) -> i64,
{
    pub fn new(x: i64, h: H, ps: &'a [u64]) -> Self {
        let stack = vec![(1, 1, 0)];
        Self {
            x,
            h,
            ps,
            stack,
            max_len: 1,
        }
    }
}
impl<H, const M: i64> Iterator for PowerfulExtAlt<'_, H, M>
where
    H: Fn(i64, i64) -> i64,
{
    type Item = (i64, i64);
    fn next(&mut self) -> Option<Self::Item> {
        while let Some((n, hn, i)) = self.stack.pop() {
            if let Some(&p) = self.ps.get(i) {
                let xdivn = (self.x / n) as u64;
                if p * p > xdivn {
                    return Some((n, hn));
                }
                self.stack.push((n, hn, i + 1));
                let mut pp = p * p;
                let mut e = 2;
                while pp <= xdivn {
                    self.stack
                        .push((n * pp as i64, (hn * (self.h)(p as i64, e)) % M, i + 1));
                    if pp > xdivn / p {
                        break;
                    }
                    pp *= p;
                    e += 1;
                }
                if self.stack.len() > self.max_len {
                    self.max_len = self.stack.len();
                }
            } else {
                return Some((n, hn));
            }
        }
        None
    }
}
