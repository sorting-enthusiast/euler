use core::ptr::copy_nonoverlapping;
use std::cmp::min;
const MOD30_TO_BIT8: [u64; 30] = [
    !0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7,
];
const MOD30_TO_MASK: [u8; 30] = [
    0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 4, 4, 8, 8, 8, 8, 16, 16, 32, 32, 32, 32, 64, 64, 64, 64, 64,
    64, 128,
];
pub const BIT64TOVAL240: [u64; 64] = [
    1, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 49, 53, 59, 61, 67, 71, 73, 77, 79, 83, 89,
    91, 97, 101, 103, 107, 109, 113, 119, 121, 127, 131, 133, 137, 139, 143, 149, 151, 157, 161,
    163, 167, 169, 173, 179, 181, 187, 191, 193, 197, 199, 203, 209, 211, 217, 221, 223, 227, 229,
    233, 239,
];
const DIFF: [u64; 8] = [6, 4, 2, 4, 2, 4, 6, 2];

#[inline(always)]
const fn unmark(x: usize, bitset: *mut u8) {
    unsafe {
        *bitset.add(x / 30) &= !MOD30_TO_MASK[x % 30];
    }
}

#[inline(always)]
const fn marked(x: usize, bitset: *const u8) -> bool {
    unsafe { *bitset.add(x / 30) & MOD30_TO_MASK[x % 30] != 0 }
}

const fn extend(bitmap: *mut u8, length: u64, n: u64) -> u64 {
    let length_div_30 = length / 30;
    let mut offset = length_div_30;
    let mut i = 1;
    while i < (n / length) {
        unsafe { copy_nonoverlapping(bitmap, bitmap.add(offset as usize), length_div_30 as usize) };
        offset += length_div_30;
        i += 1;
    }
    let rem = n % length;
    unsafe { copy_nonoverlapping(bitmap, bitmap.add(offset as usize), (rem / 30) as usize) };
    offset += rem / 30;
    let rem_30 = rem % 30;
    if rem_30 > 0 {
        unsafe {
            *bitmap.add(offset as usize) =
                *bitmap.add((rem / 30) as usize) & (0xff >> (7 - MOD30_TO_BIT8[rem_30 as usize]));
        }
    }
    n
}

fn delete(bitmap: *mut u8, p: u64, length: u64) {
    let mut pr240_on_30 = [0u32; 64];
    let mut pr240bit8mask = [0u8; 64];
    for r in 0..64 {
        let t = p * BIT64TOVAL240[r];
        pr240_on_30[r] = (t / 30) as u32;
        pr240bit8mask[r] = !MOD30_TO_MASK[t as usize % 30];
    }
    let bitmap64: *const u64 = bitmap.cast(); // as *const u64;
    let mut maxf = length / p;
    if maxf & 1 == 0 {
        maxf -= 1;
    }
    let kmin = p / 240;
    let kmax = maxf / 240;
    let bit64 = (((maxf % 240) / 30) << 3).wrapping_add(MOD30_TO_BIT8[(maxf % 240 % 30) as usize]);
    let mut base_on_30 = (kmin * p) << 3;
    let length_to_1_on_3 = (length as f64).cbrt() as u64;

    if p > length_to_1_on_3 {
        for k in kmin..kmax {
            let mut bitset = unsafe { *bitmap64.add(k as usize) };
            while bitset != 0 {
                let r = bitset.trailing_zeros() as usize;
                let c_on_30 = base_on_30 + u64::from(pr240_on_30[r]);
                let cbit8mask = pr240bit8mask[r];
                bitset &= bitset - 1;
                unsafe {
                    *bitmap.add(c_on_30 as usize) &= cbit8mask;
                }
            }

            base_on_30 += p << 3;
        }
        let mut bitset = unsafe { *bitmap64.add(kmax as usize) & ((!0u64) >> (63 - bit64)) };
        while bitset != 0 {
            let r = bitset.trailing_zeros() as usize;
            let c_on_30 = base_on_30 + u64::from(pr240_on_30[r]);
            let cbit8mask = pr240bit8mask[r];
            bitset &= bitset - 1;
            unsafe {
                *bitmap.add(c_on_30 as usize) &= cbit8mask;
            }
        }

        unmark(p as usize, bitmap);
        return;
    }
    let mut maxf_on_p = maxf / p;
    if maxf_on_p & 1 == 0 {
        maxf_on_p -= 1;
    }
    let kmid = maxf_on_p / 240;
    let bit64_mid =
        (maxf_on_p % 240 / 30 * 8).wrapping_add(MOD30_TO_BIT8[(maxf_on_p % 240 % 30) as usize]);
    let mut c_stack = Vec::with_capacity((((maxf - 1) / 30 + 1) * 8) as usize);
    for k in kmin..kmid {
        let mut bitset = unsafe { *bitmap64.add(k as usize) };
        while bitset != 0 {
            let r = bitset.trailing_zeros() as usize;
            let c_on_30 = base_on_30 + u64::from(pr240_on_30[r]);
            let cbit8mask = pr240bit8mask[r];
            bitset &= bitset - 1;
            c_stack.push((c_on_30 << 8) | u64::from(cbit8mask));
        }

        base_on_30 += p << 3;
    }
    let mut bitset = unsafe { *bitmap64.add(kmid as usize) & ((!0u64) >> (63 - bit64_mid)) };
    while bitset != 0 {
        let r = bitset.trailing_zeros() as usize;
        let c_on_30 = base_on_30 + u64::from(pr240_on_30[r]);
        let cbit8mask = pr240bit8mask[r];
        bitset &= bitset - 1;
        c_stack.push((c_on_30 << 8) | u64::from(cbit8mask));
    }
    bitset = if bit64_mid == 63 {
        0
    } else {
        unsafe { *bitmap64.add(kmid as usize) & ((!0u64) << (1 + bit64_mid)) }
    };
    if kmax > kmid {
        while bitset != 0 {
            let r = bitset.trailing_zeros() as usize;
            let c_on_30 = base_on_30 + u64::from(pr240_on_30[r]);
            let cbit8mask = pr240bit8mask[r];
            bitset &= bitset - 1;
            unsafe {
                *bitmap.add(c_on_30 as usize) &= cbit8mask;
            }
        }
        base_on_30 += p << 3;
        for k in kmid + 1..kmax {
            let mut bitset = unsafe { *bitmap64.add(k as usize) };
            while bitset != 0 {
                let r = bitset.trailing_zeros() as usize;
                let c_on_30 = base_on_30 + u64::from(pr240_on_30[r]);
                let cbit8mask = pr240bit8mask[r];
                bitset &= bitset - 1;
                unsafe {
                    *bitmap.add(c_on_30 as usize) &= cbit8mask;
                }
            }
            base_on_30 += p << 3;
        }
        bitset = unsafe { *bitmap64.add(kmax as usize) & ((!0u64) >> (63 - bit64)) };
    } else {
        bitset &= unsafe { *bitmap64.add(kmid as usize) & ((!0u64) >> (63 - bit64)) };
    }
    while bitset != 0 {
        let r = bitset.trailing_zeros() as usize;
        let c_on_30 = base_on_30 + u64::from(pr240_on_30[r]);
        let cbit8mask = pr240bit8mask[r];
        bitset &= bitset - 1;
        unsafe {
            *bitmap.add(c_on_30 as usize) &= cbit8mask;
        }
    }
    while let Some(t) = c_stack.pop() {
        unsafe {
            *bitmap.add((t >> 8) as usize) &= (t & 0xff) as u8;
        }
    }

    unmark(p as usize, bitmap);
}

#[must_use]
pub fn sift(n: u64) -> Vec<u64> {
    let bitmap_size = ((n as usize / 30) >> 3) + 1;
    let mut primes = Vec::with_capacity((n as f64 / (n as f64).log(3.0)) as usize);
    primes.extend([2, 3, 5]);
    let mut bitmap64 = vec![0u64; bitmap_size].into_boxed_slice();
    bitmap64[0] = 0xFF;

    let bitmap = bitmap64.as_ptr().cast::<u8>();

    let mut length = 30;
    let mut p = 7;
    let mut p_squared = 49;
    let mut p_index = 1;
    if n < 30 {
        bitmap64[0] &= 0xFF >> (7 - MOD30_TO_BIT8[n as usize]);
        length = n;
    }
    while p_squared <= n {
        if length < n {
            //dbg!(p);
            length = extend(bitmap.cast_mut(), length, min(p * length, n));
            if length == n {
                unmark(1, bitmap.cast_mut());
            }
        }
        delete(bitmap.cast_mut(), p, length);
        primes.push(p);
        loop {
            p += DIFF[p_index as usize & 7];
            p_index += 1;
            if marked(p as usize, bitmap) {
                break;
            }
        }
        p_squared = p * p;
    }
    if length < n {
        _ = extend(bitmap.cast_mut(), length, n);
    }
    unmark(1, bitmap.cast_mut());
    //let bitmap64 = bitmap.as_ptr() as *const u64;
    let mut base = 0;
    for k in 0..bitmap_size {
        let mut bitset = bitmap64[k];
        while bitset != 0 {
            let r = bitset.trailing_zeros() as usize;
            primes.push(base + BIT64TOVAL240[r]);
            bitset &= bitset - 1;
        }
        base += 240;
    }
    primes
}
