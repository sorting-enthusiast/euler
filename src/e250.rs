use std::time::Instant;

pub fn main() {
    const MOD: i64 = 1e16 as i64;
    let start = Instant::now();
    let mut row1 = [0; 250];
    let mut row2 = [0; 250];
    row1[0] = 1;
    let (mut prev, mut cur) = (&mut row1, &mut row2);
    for i in 1..=250_250 {
        let residue = modexp::<250>(i, i) as usize;
        *cur = *prev;
        for r in 0..residue {
            cur[r] += prev[r + 250 - residue];
            if cur[r] >= MOD {
                cur[r] -= MOD;
            }
        }
        for r in residue..250 {
            cur[r] += prev[r - residue];
            if cur[r] >= MOD {
                cur[r] -= MOD;
            }
        }
        core::mem::swap(&mut prev, &mut cur);
    }
    let res = prev[0] - 1;
    println!("{:?}", start.elapsed());
    println!("{res}");
}
const fn modexp<const MODULO: u64>(mut x: u64, mut exp: u64) -> u64 {
    let mut r = 1;
    while exp > 1 {
        if exp & 1 == 1 {
            r = (r * x) % MODULO;
        }
        x = (x * x) % MODULO;
        exp >>= 1;
    }
    (r * x) % MODULO
}
