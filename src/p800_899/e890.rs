use num_bigint::BigUint;
const MOD: u32 = 1e9 as u32 + 7;
fn conv(a: &[u32], b: &[u32], c: &mut [u32]) {
    c.fill(0);
    for (i, &ai) in a.iter().enumerate() {
        for (j, &bj) in b.iter().enumerate() {
            let (ai, bj) = (u64::from(ai), u64::from(bj));
            let cij = (ai * bj) % u64::from(MOD);
            c[i + j] += cij as u32;
            if c[i + j] >= MOD {
                c[i + j] -= MOD;
            }
        }
    }
}

// nice problem
// gf F(x) = \prod_{k\ge 0} {\frac1{1-x^{2^k}}}
// (1-x)F(x) = F(x^2)
// similar to the algorithm for the n'th coefficient of a rational polynomial: https://arxiv.org/pdf/2008.08822
// computed using https://qiita.com/ryuhe1/items/185e1a283f13ac638a53#%E4%BA%8C%E5%86%AA%E5%88%86%E5%89%B2%E6%95%B0
pub fn main() {
    let start = std::time::Instant::now();
    let mut N = BigUint::from(7u32).pow(777);

    let mut p = vec![1];
    let mut q = vec![1];
    let mut tmp = vec![];
    while N != BigUint::ZERO {
        let lsb = usize::from(N.bit(0));
        q.push(unsafe { *q.last().unwrap_unchecked() });
        for i in (1..q.len() - 1).rev() {
            q[i] += q[i - 1];
            if q[i] >= MOD {
                q[i] -= MOD;
            }
        }
        tmp.resize(q.len() + p.len() - 1, 0);
        conv(&p, &q, &mut tmp);
        p.clear();
        p.extend(tmp[lsb..].iter().step_by(2));
        N >>= 1;
    }
    println!("res = {}, took {:?}", p[0], start.elapsed());
}
