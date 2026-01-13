use std::time::Instant;

use crate::utils::{
    FIArray::FIArrayI64, fast_divisor_sums::sum_floors_fast, powerful_numbers::PowerfulExtSkipZero,
};
const N: i64 = 1e14 as i64;

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
    let mut cache = FIArrayI64::new(N);
    for (n, hn) in PowerfulExtSkipZero::<_, { i64::MAX }>::new(N, h) {
        let i = cache.get_index(N / n);
        if cache.arr[i] == 0 {
            cache.arr[i] = sum_floors_fast((N / n) as _) as i64;
        }
        sum += hn * cache.arr[i];
    }

    let end = start.elapsed();
    println!("res = {sum}, took {end:?}");
}
