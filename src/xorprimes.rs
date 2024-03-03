use crate::squarefree_square_plus_one::BitArray;
const fn xor_product(X: usize, Y: usize) -> usize {
    let mut res = 0;
    let mut x = if X < Y { X } else { Y };
    let y = if X > Y { X } else { Y };
    let mut shift_amt = x.trailing_zeros();
    x >>= shift_amt;
    while x > 0 {
        res ^= y << shift_amt;
        x >>= 1;
        let tmp = x.trailing_zeros();
        x >>= tmp;
        shift_amt += 1 + tmp;
    }
    res
}
const fn naive_xor_product(X: usize, Y: usize) -> usize {
    let mut res = 0;
    let mut x = if X < Y { X } else { Y };
    let y = if X > Y { X } else { Y };
    let mut shift_amt = 0;
    while x > 0 {
        if x & 1 == 1 {
            res ^= y << shift_amt;
        }
        x >>= 1;
        shift_amt += 1;
    }
    res
}
pub fn main() {
    const X: usize = 2374659487;
    const Y: usize = 25072475244257;
    dbg!(xor_product(X, Y));
    dbg!(naive_xor_product(X, Y));
    assert_eq!(xor_product(X, Y), naive_xor_product(X, Y));
}
