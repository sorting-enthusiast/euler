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
fn is_xor_prime(x: usize) -> Result<(), (usize, usize)> {
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
                let multiple = xor_product(num, i);
                if multiple > new_lim {
                    continue;
                }
                sieve.set(multiple >> 1);
            }
        }
    }
    None
}

pub fn main() {
    dbg!(is_xor_prime(5));
    //assert_eq!(499, sieve_of_eratosthenes_xor_primes(1000).unwrap());
    dbg!(sieve_of_eratosthenes_xor_primes(150_000_000));
}
