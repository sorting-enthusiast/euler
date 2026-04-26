pub fn main() {
    const N: usize = 1e4 as _;
    let mut proper_divisors = vec![1; N];
    proper_divisors[0] = 0;
    proper_divisors[1] = 0;
    for d in 2..N {
        for m in (2 * d..N).step_by(d) {
            proper_divisors[m] += d;
        }
    }
    let mut sum = 0;
    for n in 2..N {
        let pdn = proper_divisors[n];
        if (n + 1..N).contains(&pdn) && proper_divisors[pdn] == n {
            sum += n + pdn;
        }
    }
    dbg!(sum);
}
