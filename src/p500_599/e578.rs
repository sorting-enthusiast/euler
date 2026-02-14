use crate::utils::multiplicative_function_summation::count_squarefree;

pub fn main() {
    const N: usize = 1e6 as _;
    let sqf = count_squarefree(N);
    dbg!(sqf[N]);
    todo!()
}
