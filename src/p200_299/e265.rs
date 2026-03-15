// brute force each 5-bit subsequence
fn is_de_bruijn(b: u32) -> bool {
    let b = u64::from(b);
    let y = (b << 32) | b;
    let mut seen = 0u32;
    for i in 0..32 {
        let w = (y >> i) as u32 & 31;
        seen ^= 1 << w;
    }
    seen == u32::MAX
}
pub fn main() {
    let start = std::time::Instant::now();
    // restricts i to iterate over numbers with msb's 000001 and lsb 1
    let sum = ((1 << 26) | 1..1 << 27)
        .step_by(2)
        .filter(|&i| is_de_bruijn(i))
        .map(u64::from)
        .sum::<u64>();
    println!("res = {sum}, took {:?}", start.elapsed());
}
