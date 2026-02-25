const MOD: u32 = 1e9 as u32 + 7;
const N: usize = 1e14 as _;
pub fn main() {
    let start = std::time::Instant::now();
    let res = (2..=10).fold(0, |acc, k| (acc + f(k)) % MOD);
    println!("res = {res}, took {:?}", start.elapsed());
}
// essentially the same idea as 890
fn f(k: usize) -> u32 {
    let mut n = N;
    let mut p = vec![1];
    let q = vec![1; k];
    let mut q0 = q.clone();
    let mut tmp = vec![];
    while n != 0 {
        tmp.resize(q0.len() + q.len() - 1, 0);
        conv(&q, &q0, &mut tmp);
        core::mem::swap(&mut q0, &mut tmp);

        tmp.resize(q0.len() + p.len() - 1, 0);
        conv(&p, &q0, &mut tmp);
        p.clear();
        p.extend(tmp[n % k..].iter().step_by(k));
        n /= k;
    }
    p[0]
}
fn conv(a: &[u32], b: &[u32], c: &mut [u32]) {
    c.fill(0);
    for (i, &ai) in a.iter().enumerate() {
        for (j, &bj) in b.iter().enumerate() {
            let (ai, bj) = (ai as u64, bj as u64);
            let cij = (ai * bj) % MOD as u64;
            c[i + j] += cij as u32;
            if c[i + j] >= MOD {
                c[i + j] -= MOD;
            }
        }
    }
}
