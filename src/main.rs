pub mod longest_collatz_chain;
use crate::longest_collatz_chain::{collatz_chain_length, Collatz};
const N: usize = 1_000_000_000;

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
    }
    println!("{:?}", (max_start, max_chain));
    println!("{:?}", (837799, collatz_chain_length(837799)));
}
