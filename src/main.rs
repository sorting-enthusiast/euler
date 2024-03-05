pub mod longest_collatz_chain;
pub mod pandigital_products;
pub mod squarefree_square_plus_one;
pub mod xorprimes;
const fn modular_exponentiation(mut x: u128, mut exp: u128, modulo: u128) -> u128 {
    let mut r = 1;
    while exp > 1 {
        if exp & 1 == 1 {
            r = (r * x) % modulo;
        }
        x = (x * x) % modulo;
        exp >>= 1;
    }
    (r * x) % modulo
}
pub fn main() {
    //squarefree_square_plus_one::main();
    //longest_collatz_chain::main();
    xorprimes::main();
}
