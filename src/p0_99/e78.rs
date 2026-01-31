pub fn main() {
    const lim: usize = 6e4 as _;
    let mut gf = vec![0u64; lim + 1];
    gf[0] = 1;
    for x in 1..=lim {
        for n in x..=lim {
            gf[n] += gf[n - x];
            if gf[n] >= 1000000 {
                gf[n] -= 1000000;
            }
        }
        if gf[x] == 0 {
            println!("{x}");
            break;
        }
    }
    dbg!(gf[5]);
}
