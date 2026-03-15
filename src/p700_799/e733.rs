use crate::utils::fenwick::FenwickTreeUsizeMod;
const MOD: usize = 1e9 as usize + 7;
const N: usize = 1e6 as _;

pub fn main() {
    let start = std::time::Instant::now();
    let mut sums: [_; 3] =
        core::array::from_fn(|_| FenwickTreeUsizeMod::<MOD>::new(10_000_019 + 1, 0));
    let mut cnts: [_; 3] =
        core::array::from_fn(|_| FenwickTreeUsizeMod::<MOD>::new(10_000_019 + 1, 0));
    let mut ret = 0;
    let mut ai = 1;
    for _ in 0..N {
        ai *= 153;
        ai %= 10_000_019;
        ret += (cnts[2].sum(ai) * ai) % MOD + sums[2].sum(ai);
        ret %= MOD;

        let new_cnts = [1, cnts[0].sum(ai), cnts[1].sum(ai)];
        let new_sums = [
            ai,
            (new_cnts[1] * ai) % MOD + sums[0].sum(ai),
            (new_cnts[2] * ai) % MOD + sums[1].sum(ai),
        ];

        for j in 0..3 {
            cnts[j].add(ai, new_cnts[j]);
            sums[j].add(ai, new_sums[j]);
        }
    }
    println!("ret = {ret}, took {:?}", start.elapsed());
}
