pub fn main() {
    let mut sum = 0;
    let a003313 = std::fs::read_to_string(r"src\p100_199\A003313.txt").unwrap();
    for l in a003313.lines() {
        sum += l.split(' ').next_back().unwrap().parse::<u32>().unwrap();
    }
    dbg!(sum);
}
