pub fn main() {
    for i in 0..25 {
        println!("f({i}) = {}", aviv_optimized(i));
    }
    println!(
        "{},\n{}",
        aviv_optimized(10_000_000_000_000_000_000_000_000_u128),
        aviv(10_000_000_000_000_000_000_000_000_u128)
    );
}
pub const fn aviv_optimized(mut num: u128) -> u64 {
    let mut ans = 1;
    let mut dp = 0;
    while num != 0 {
        let tz = num.trailing_zeros() as u64;
        dp += tz * ans;
        ans += dp;
        num >>= tz + 1;
    }
    ans
}
pub const fn aviv(mut num: u128) -> u64 {
    let mut ans = 1;
    let mut dp = 0;
    while num != 0 {
        if num & 1 == 0 {
            dp += ans;
        } else {
            ans += dp;
        }
        num >>= 1;
    }
    ans
}
