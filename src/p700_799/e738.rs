use crate::utils::{FIArray::FIArrayI64, multiplicative_function_summation::dirichlet_inv_i64};

const N: usize = 10;
const MOD: i64 = 1e9 as i64 + 7;

pub fn main() {
    let mut unit = FIArrayI64::unit(N);
    for e in &mut unit.arr {
        *e = 2 - *e;
    }
    dbg!(dirichlet_inv_i64(&unit, N as _));
}
