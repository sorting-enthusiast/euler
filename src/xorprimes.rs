use crate::squarefree_square_plus_one::BitArray;
const fn xor_product(X: u128, Y: u128) -> u128 {
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
fn is_xor_prime(x: u128) -> Result<(), (u128, u128)> {
    for i in 2..x - 1 {
        for j in i..x - 1 {
            if x == xor_product(i, j) {
                return Err((i, j));
            }
        }
    }
    Ok(())
}
fn sieve_of_eratosthenes_xor_primes(limit: usize) -> Option<usize> {
    let mut sieve = BitArray::zeroed((limit + 1) / 2);
    let new_lim = 1 << (usize::BITS - 1 - limit.leading_zeros());
    let mut counter = 1;
    for num in (3..=new_lim).step_by(2) {
        if !sieve.get(num >> 1) {
            counter += 1;
            if counter == 5_000_000 {
                return Some(num);
            }
            let _limit = new_lim >> ((usize::BITS - 1) ^ num.leading_zeros());

            for i in (num.._limit).step_by(2) {
                let multiple = xor_product(num as u128, i as u128);
                if multiple > new_lim as u128 {
                    continue;
                }
                sieve.set(multiple as usize >> 1);
            }
        }
    }
    None
}
fn xor_power_rec(exp: u128) -> u128 {
    if exp == 1 {
        return 11;
    }
    let tmp = xor_power_rec(exp / 2);
    let res = xor_product(tmp, tmp);
    if exp & 1 == 0 {
        res
    } else {
        xor_product(res, 11)
    }
}
const N: u128 = 1_000_000_001;
fn xor_power(mut exp: u128) -> u128 {
    let mut x = 11;
    let mut y = 1;
    while exp >= 1 {
        if exp & 1 == 1 {
            y = xor_product(x, y);
        }
        x = xor_product(x, x);
        exp >>= 1;
    }
    y
}
pub fn main() {
    dbg!(xor_power(3));
    dbg!(xor_power_rec(3));
    let exp: u128 = 8u128.pow(12) * 12u128.pow(8);
    dbg!(exp);
    dbg!(xor_power(exp));
    dbg!(xor_power_rec(exp));
    //assert_eq!(499, sieve_of_eratosthenes_xor_primes(1000).unwrap());
    dbg!(sieve_of_eratosthenes_xor_primes(150_000_000));
}
