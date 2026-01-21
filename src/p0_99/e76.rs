pub fn main() {
    let mut gf = vec![0; 100 + 1];
    gf[0] = 1;
    for x in 1..=100 {
        for n in x..=100 {
            gf[n] += gf[n - x];
        }
    }
    dbg!(gf[5] - 1);
    dbg!(gf[100] - 1);
}
