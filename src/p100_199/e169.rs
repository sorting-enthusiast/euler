pub fn main() {
    println!(
        "{},\n{},\n{}",
        aviv_optimized(10_000_000_000_000_000_000_000_000_u128),
        aviv(10_000_000_000_000_000_000_000_000_u128),
        gf_based_opt(10_000_000_000_000_000_000_000_000_u128)
    );
}
const fn aviv_optimized(mut num: u128) -> u64 {
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
const fn aviv(mut num: u128) -> u64 {
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
// using standard method for computation of [x^n]F(x) where F(x)=\prod_{k\ge 0}\frac{P(x^{2^k})}{Q(x^{2^k})}
// in this case P(x) = 1 + x + x^2 and Q(x) = 1, so the code can be simplified manyfold
fn gf_based(mut n: u128) -> u64 {
    const LEN: usize = 4;
    let mut P0 = [0; LEN];
    let mut tmp = [0; LEN];
    P0[0] = 1;
    while n > 0 {
        for i in (2..LEN).rev() {
            P0[i] += P0[i - 1] + P0[i - 2];
        }
        P0[1] += P0[0];
        for (i, v) in ((n & 1) as usize..LEN).step_by(2).enumerate() {
            tmp[i] = P0[v];
        }
        P0.copy_from_slice(&tmp);
        tmp.fill(0);
        n >>= 1;
    }
    P0[0]
}
// optimization of above method, using the fact that only 2 coefficients are ever in use
// exactly what aviv's solution does
const fn gf_based_opt(mut n: u128) -> u64 {
    let [mut x0, mut x1] = [1, 0];
    while n > 0 {
        if n & 1 == 0 {
            x1 += x0;
        } else {
            x0 += x1;
        }
        n >>= 1;
    }
    x0
}
