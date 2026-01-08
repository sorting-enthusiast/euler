use crate::utils::FIArray::FIArray;

#[must_use]
pub fn count_signature(sig: &[u8], lim: usize, pis: &FIArray, primes: &[u64]) -> usize {
    fn rec(sig: &[u8], lim: usize, pis: &FIArray, primes: &[u64], j: usize) -> usize {
        unsafe { core::hint::assert_unchecked(!sig.is_empty()) };
        let e = sig[0];
        if sig.len() == 1 {
            return pis[iroot(lim, e)] - j;
        }
        let rem: u32 = sig.iter().copied().map(u32::from).sum();
        let mut cnt = 0;
        for (k, &p) in primes[j..].iter().enumerate() {
            let p = p as usize;
            if p.pow(rem) > lim {
                break;
            }
            cnt += rec(&sig[1..], lim / p.pow(e.into()), pis, primes, j + k + 1);
        }
        cnt
    }
    rec(sig, lim, pis, primes, 0)
}
#[must_use]
pub const fn iroot(x: usize, k: u8) -> usize {
    let k = k as usize;
    let mut rt = 1usize << (1 + x.ilog2().div_ceil(k as _));
    let mut x_div_rtk1 = x / rt.pow(k as u32 - 1);
    while rt > x_div_rtk1 {
        rt = (rt * (k - 1) + x_div_rtk1) / k;
        x_div_rtk1 = x / rt.pow(k as u32 - 1);
    }
    rt
}
