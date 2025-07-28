use crate::utils::{
    FIArray::FIArrayI64, multiplicative_function_summation::divisor_summatory_i64,
    powerful_numbers::PowerfulExt,
};
const N: i64 = 1e15 as i64;

// powerful number trick
// f(p) = d(p) => g = d
// h(p^e) = f(p^e) - 2f(p^(e-1)) + f(p^(e-2)) = 0 for odd e, p^(e/2 - 1) * (p - 1) for even e
pub fn main() {
    let start = std::time::Instant::now();

    let h = |p: i64, e| {
        unsafe { core::hint::assert_unchecked(e > 1) };
        if e & 1 == 1 {
            0
        } else {
            p.pow((e as u32 >> 1) - 1) * (p - 1)
        }
    };
    let mut sum = 0;
    let mut cache = FIArrayI64::new(N);
    for (n, hn) in PowerfulExt::<_, { i64::MAX }>::new(N, h).filter(|&(_, hn)| hn != 0) {
        let index = cache.get_index(N / n);
        if cache.arr[index] == 0 {
            cache.arr[index] = divisor_summatory_i64(N / n);
        }
        sum += hn * cache.arr[index]; //divisor_summatory_i64(N / n);
    }
    let end = start.elapsed();
    println!("res = {sum}, took {end:?}");
}
