use crate::utils::multiplicative_function_summation::totient_sieve;
// huh?
pub fn main() {
    const N: usize = 4e7 as _;
    let start = std::time::Instant::now();

    let mut totients = totient_sieve(N + 1);
    totients[1] = 1;
    totients[2] = 2;
    let mut sum = 0;

    for i in 3..=N {
        let phi = totients[i] as usize;
        totients[i] = 1 + totients[phi];

        if phi == i - 1 && totients[i] == 25 {
            sum += i;
        }
    }

    let end = start.elapsed();
    println!("res = {sum}, took {end:?}");
}
