pub fn main() {
    const V: [u64; 2] = [4, 2];
    let mut seq = vec![];
    let mut which = 0;
    let mut run = V[which];
    let mut len = 0;
    for i in 0..100 {
        for _ in 0..run {
            seq.push(V[which]);
        }
        len += run * V[which];
        which ^= 1;
        run = seq[i + 1];
        if i.trailing_zeros() >= 10 {
            dbg!(len);
        }
    }
    dbg!(seq[..48].iter().sum::<u64>());
    for e in &seq[..48] {
        print!("{e}, ");
    }
    println!();
}
