use crate::utils::math::extended_gcd;

pub fn main() {
    const N: usize = 1e7 as _;
    let start = std::time::Instant::now();
    let mut divisors = vec![vec![1]; N + 1];
    for d in 2..=N.isqrt() {
        for m in (d..=N).step_by(d) {
            divisors[m].push(d);
        }
    }
    dbg!(start.elapsed());

    let mut sum = 0;
    for n in 2..=N {
        // find maximal x s.t. x(x-1)=0 mod n

        let mut max_root = 1;
        for &d in &divisors[n] {
            let n = n as i64;
            let d = d as i64;
            let nd = n / d;
            if d > nd {
                break;
            }
            let (g, x, _) = extended_gcd(d, nd);
            if g != 1 {
                continue;
            }
            // dx + (n/d)y = 1
            // dx - 1 = -y(n/d)
            let mut root = (d * x) % n;
            if root < 0 {
                root += n;
            }
            if root > max_root {
                max_root = root;
            }
            root = (n + 1 - root) % n;
            if root > max_root {
                max_root = root;
            }
        }

        sum += max_root;
    }
    dbg!(sum);
    dbg!(start.elapsed());
}
