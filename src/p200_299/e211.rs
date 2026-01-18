pub fn main() {
    const N: usize = 6.4e7 as _;
    let start = std::time::Instant::now();
    let mut s2 = vec![1; N];
    for d in 2..N {
        let dd = d * d;
        for m in (d..N).step_by(d) {
            s2[m] += dd;
        }
    }
    dbg!(start.elapsed());
    let mut res = 1;
    for n in 2..N {
        if s2[n].isqrt().pow(2) == s2[n] {
            res += n;
        }
    }
    dbg!(start.elapsed(), res);
}
