use std::time::Instant;

use crate::utils::{
    FIArray::FIArray, fast_divisor_sums::divisor_summatory, powerful_numbers::PowerfulExtSkipZero,
};
const N: usize = 1e14 as usize;

// powerful number trick
// f(p) = d(p) => g = d
// h(p^e) = f(p^e) - 2f(p^(e-1)) + f(p^(e-2)) = 2^(e-2)
pub fn main() {
    let start = Instant::now();

    let h = |_p, e| {
        unsafe { core::hint::assert_unchecked(e > 1) };
        1i64 << (e - 2)
    };
    let mut sum = 0;
    let mut cache = FIArray::new(N);
    for (n, hn) in PowerfulExtSkipZero::<_, { i64::MAX }>::new(N as _, h) {
        let i = cache.get_index(N / n as usize);
        if cache.arr[i] == 0 {
            cache.arr[i] = divisor_summatory(N / n as usize);
        }
        sum += hn * cache.arr[i] as i64;
    }

    let end = start.elapsed();
    println!("res = {sum}, took {end:?}");
}
