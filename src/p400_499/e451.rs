use crate::utils::math::extended_gcd;

pub fn main() {
    const N: usize = 2e7 as _;
    let start = std::time::Instant::now();
    let mut divisors = vec![vec![1]; N + 1];
    for d in 2..=N.isqrt() {
        for m in (d..=N).step_by(d) {
            divisors[m].push(d);
        }
    }
    dbg!(start.elapsed());

    let mut sum = 0;
    for n in 3..=N {
        // find maximal x s.t. (x+1)(x-1)=0 mod n
        let mut max_root = 1;
        for &d in &divisors[n] {
            let n = n as i64;
            let d = d as i64;
            let nd = n / d;
            if d > nd {
                break;
            }
            let (g, mut x, _) = extended_gcd(d, nd);
            if 2 % g != 0 {
                continue;
            }
            // dx+(n/d)y = g
            // dx(-2/g) + (n/d)y(-2/g) = -2
            // (n/d)y(-2/g)=dx(-2/g) + 2

            x *= -2 / g;

            let mut root = (d * x + 1) % n;
            if root < 0 {
                root += n;
            }
            for _ in 0..g {
                let larger_root = if root < n - root { n - root } else { root };
                if larger_root != n - 1 && larger_root > max_root {
                    max_root = larger_root;
                }
                root += n / g;
                if root >= n {
                    root -= n;
                }
            }
        }

        sum += max_root;
    }
    dbg!(sum);
    dbg!(start.elapsed());
}
