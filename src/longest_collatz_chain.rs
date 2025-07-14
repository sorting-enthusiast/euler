//Which starting number, under one million, produces the longest chain?
use std::collections::HashMap;
#[derive(Default)]
pub struct DynamicProgrammingCollatz {
    lookup_table: HashMap<usize, usize>,
}
impl DynamicProgrammingCollatz {
    pub fn chain_length(&mut self, mut x: usize) -> usize {
        let original_x = x;
        let mut seq_length = 1;
        if x & 1 == 0 {
            let shift_amt = x.trailing_zeros() as usize;
            x >>= shift_amt;
            seq_length += shift_amt;
        }
        while x > 1 {
            if let Some(&len) = self.lookup_table.get(&x) {
                return seq_length + len - 1;
            }
            x += (x + 1) >> 1;
            let shift_amt = x.trailing_zeros() as usize;
            x >>= shift_amt;
            seq_length += shift_amt + 2;
        }
        self.lookup_table.insert(original_x, seq_length);
        seq_length
    }
}
//use for compile time evaluation
#[must_use]
pub const fn collatz_chain_length(mut x: usize) -> usize {
    let mut seq_length = 1;
    if x & 1 == 0 {
        let shift_amt = x.trailing_zeros() as usize;
        x >>= shift_amt; //equivalent to repeatedly dividing by 2 until reaching an odd number, adds shift_amt to seq_length

        seq_length += shift_amt;
    }
    while x > 1 {
        x += (x + 1) >> 1; //equivalent to (3x+1)/2, adds 2 to seq_length

        let shift_amt = x.trailing_zeros() as usize;
        x >>= shift_amt; //equivalent to repeatedly dividing by 2 until reaching an odd number, adds shift_amt to seq_length

        seq_length += shift_amt + 2;
    }
    seq_length
}
pub trait Collatz {
    fn collatz_chain_length(self) -> usize;
}
macro_rules! gen_collatz_impl {
    ($($type:ty),+) => { $(
        impl Collatz for $type {
            fn collatz_chain_length(self) -> usize {
                let mut x = self;
                let mut seq_length = 1;
                if x & 1 == 0 {
                    let shift_amt = x.trailing_zeros() as usize;
                    x >>= shift_amt; //equivalent to repeatedly dividing by 2 until reaching an odd number, adds shift_amt to seq_length

                    seq_length += shift_amt;
                }
                while x > 1 {
                    x += (x + 1) >> 1; //equivalent to (3x+1)/2, adds 2 to seq_length

                    let shift_amt = x.trailing_zeros() as usize;
                    x >>= shift_amt; //equivalent to repeatedly dividing by 2 until reaching an odd number, adds shift_amt to seq_length

                    seq_length += shift_amt + 2;
                }
                seq_length
            }
        }

    )+ };
}
gen_collatz_impl!(u8, u16, u32, u64, u128, usize);

const N: usize = 1_000_000;
const MAX_UNDER_MILLION: usize = 837_799;
const MAX_CHAIN_LENGTH: usize = collatz_chain_length(MAX_UNDER_MILLION);
pub fn main() {
    let mut max_start = 1;
    let mut max_chain = 1;
    //chain(2k) = 1 + chain(k) > chain(k) for all k, hence we can ignore all values < N/2
    for n in (N / 2)..N {
        let chain = n.collatz_chain_length();
        if chain > max_chain {
            max_chain = chain;
            max_start = n;
        }
        /* if n % 1_000_000_000 == 0 {
            println!("{n}, {max_start}, {max_chain}");
        } */
    }
    println!("{:?}", (max_start, max_chain));
    println!("{:?}", (MAX_UNDER_MILLION, MAX_CHAIN_LENGTH));
}
