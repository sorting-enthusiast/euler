const fn f(mut n: u64, d: u64) -> u64 {
    let mut ret = 0;
    while n > 9 {
        let k = n.ilog10();
        let pow10 = 10_u64.pow(k);
        if n == 10 * pow10 - 1 {
            return ret + (k as u64 + 1) * pow10;
        }
        let (a, b) = (n / pow10, n % pow10);
        ret += a * k as u64 * 10_u64.pow(k - 1);
        if a > d {
            ret += pow10;
        } else if a == d {
            ret += 1 + b;
        }
        n = b;
    }
    ret + (n >= d) as u64
}

fn sum(l: u64, r: u64, d: u64) -> u64 {
    let mid = l.midpoint(r);
    if l == mid {
        return if f(l, d) == l { l } else { 0 };
    }
    let mut ret = 0;
    // exploit the fact that both the integers and f(n, d) are monotonically increasing:
    // If l > f(mid, d), then for all n in [l, mid] n > f(mid, d) >= f(n, d), hence no need to check the left range.
    // If mid < f(l, d), then for all n in [l, mid] n < f(l, d) <= f(n, d), hence no need to check the left range.
    // The same reasoning can be applied to the right range
    if f(l, d) <= mid && l <= f(mid, d) {
        ret += sum(l, mid, d);
    }
    if f(mid, d) <= r && mid <= f(r, d) {
        ret += sum(mid, r, d);
    }
    ret
}
fn s(d: u64) -> u64 {
    sum(1, u64::MAX, d)
}
pub fn main() {
    let start = std::time::Instant::now();
    let s = (1..10).map(s).sum::<u64>();
    let end = start.elapsed();

    println!("{end:?}");
    println!("s = {s}");
}
