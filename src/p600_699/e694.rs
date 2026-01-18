use crate::utils::{math::iroot, primes::wheel_sieve};
struct CubeFull {
    x: usize,
    ps: Vec<u64>,
    stack: Vec<(usize, usize)>,
}
impl CubeFull {
    fn new(x: usize) -> Self {
        let stack = vec![(1, 0)];
        let ps = wheel_sieve(iroot::<3>(x) as _);
        Self { x, ps, stack }
    }
}
impl Iterator for CubeFull {
    type Item = usize;
    fn next(&mut self) -> Option<Self::Item> {
        while let Some((n, i)) = self.stack.pop() {
            if let Some(&p) = self.ps.get(i) {
                let p = p as usize;
                let xdivn = self.x / n;
                if p * p * p > xdivn {
                    return Some(n);
                }
                self.stack.push((n, i + 1));
                let mut pp = p * p * p;
                while pp <= xdivn {
                    self.stack.push((n * pp, i + 1));
                    if pp > xdivn / p {
                        break;
                    }
                    pp *= p;
                }
            } else {
                return Some(n);
            }
        }
        None
    }
}
pub fn main() {
    const N: usize = 1e18 as _;
    dbg!(CubeFull::new(N).map(|c| N / c).sum::<usize>());
}
