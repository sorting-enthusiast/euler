use crate::utils::{
    FIArray::FIArrayI64, fast_divisor_sums::divisor_summatory,
    powerful_numbers::PowerfulExtSkipZero,
};
const N: i64 = 1e15 as i64;

// powerful number trick
// f(p) = d(p) => g = d
// h(p^e) = f(p^e) - 2f(p^(e-1)) + f(p^(e-2)) = 0 for odd e, p^(e/2 - 1) * (p - 1) for even e, i.e. h = totient_sqrt
// can be optimized to O(n^1/3) space, using similar tricks to counting squarefree numbers
// DGF(f) = \frac{\zeta(s)^2\zeta(2s-1)}{\zeta(2s)}
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
    for (n, hn) in PowerfulExtSkipZero::<_, { i64::MAX }>::new(N, h) {
        let index = cache.get_index(N / n);
        if cache.arr[index] == 0 {
            cache.arr[index] = divisor_summatory((N / n) as _) as i64;
        }
        sum += hn * cache.arr[index];
    }
    let end = start.elapsed();
    println!("res = {sum}, took {end:?}");
}
