use itertools::Itertools;
// change according to usecase
//3221225473 = (3<<30)+1
// (7<<26)+1
pub struct NTT<const MOD: i64, const ROOT: i64, const ROOT_EXP_PLUS_1: usize>; // ROOT_EXP_PLUS_1 is included because rust is dumb
impl<const MOD: i64, const ROOT: i64, const ROOT_EXP_PLUS_1: usize>
    NTT<MOD, ROOT, ROOT_EXP_PLUS_1>
{
    const fn powmod(mut x: i64, mut exp: i64) -> i64 {
        if exp == 0 {
            return 1;
        }
        let mut r = 1;
        x %= MOD;
        while exp > 1 {
            if exp & 1 == 1 {
                r = (r * x) % MOD;
            }
            x = (x * x) % MOD;
            exp >>= 1;
        }
        (r * x) % MOD
    }
    const fn modinv(x: i64) -> i64 {
        Self::powmod(x, MOD - 2)
    }
    const ROOT_1: i64 = Self::modinv(ROOT);
    //const ROOT_PW: usize = 1 << 26;
    const ROOT_POWS: [[i64; ROOT_EXP_PLUS_1]; 2] = {
        let mut ret = [[0; ROOT_EXP_PLUS_1]; 2];
        ret[0][0] = 1;
        ret[1][0] = 1;

        ret[0][1] = ROOT;
        ret[1][1] = Self::ROOT_1;

        let mut i = 2;
        while i < ROOT_EXP_PLUS_1 {
            ret[0][i] = (ret[0][i - 1] * ret[0][i - 1]) % MOD;
            ret[1][i] = (ret[1][i - 1] * ret[1][i - 1]) % MOD;

            i += 1;
        }
        ret
    };
    pub fn ntt(a: &mut [i32], invert: bool) {
        let n = a.len();
        let mut j = 0;
        for i in 1..n {
            let mut bit = n >> 1;
            while j & bit != 0 {
                j ^= bit;
                bit >>= 1;
            }
            j ^= bit;

            if i < j {
                a.swap(i, j);
            }
        }
        let mut len = 2;
        let mut len2 = 1;
        let mut pow = ROOT_EXP_PLUS_1 - 1;
        while len <= n {
            let wlen = Self::ROOT_POWS[usize::from(invert)][pow];
            for i in (0..n).step_by(len) {
                let mut w = 1;
                for j in 0..len2 {
                    let u = a[i | j];
                    let v = ((a[i | j | len2] as i64 * w) % MOD) as i32;
                    a[i | j] = if u + v < MOD as i32 {
                        u + v
                    } else {
                        u + v - MOD as i32
                    };
                    a[i | j | len2] = if u >= v { u - v } else { u + MOD as i32 - v };
                    w = (w * wlen) % MOD;
                }
            }
            len <<= 1;
            len2 <<= 1;
            pow -= 1;
        }

        if invert {
            let n_1 = Self::modinv(n as i64);

            for x in a.iter_mut() {
                *x = ((*x as i64 * n_1) % MOD) as i32;
            }
        }
    }
    #[must_use]
    pub fn multiply(a: &[i32], b: &[i32]) -> Vec<i32> {
        let mut fa = a.iter().copied().collect_vec();
        let mut fb = b.iter().copied().collect_vec();
        let Some(trunc_a) = a.iter().rposition(|&e| e != 0) else {
            return vec![];
        };
        let Some(trunc_b) = b.iter().rposition(|&e| e != 0) else {
            return vec![];
        };

        let n = (trunc_a + 1 + trunc_b + 1 - 1).next_power_of_two();
        fa.resize(n, 0);
        fb.resize(n, 0);

        Self::ntt(&mut fa, false);
        Self::ntt(&mut fb, false);

        for i in 0..n {
            fa[i] = ((i64::from(fa[i]) * i64::from(fb[i])) % MOD) as i32;
        }
        Self::ntt(&mut fa, true);
        fa.truncate(1 + unsafe { fa.iter().rposition(|&e| e != 0).unwrap_unchecked() });
        fa
    }
}
/*const MOD: i64 = (7 << 26) + 1; //998_244_353; // u64::MAX - u32::MAX as u64 + 1
const ROOT: i64 = powmod(3, 7); //powmod(3, 119);
const ROOT_1: i64 = modinv(ROOT);
const ROOT_PW: usize = 1 << 26;
const ROOT_EXP: usize = 26;
const ROOT_POWS: [[i64; ROOT_EXP + 1]; 2] = {
    let mut ret = [[0; ROOT_EXP + 1]; 2];
    ret[0][0] = 1;
    ret[1][0] = 1;

    ret[0][1] = ROOT;
    ret[1][1] = ROOT_1;

    let mut i = 2;
    while i <= ROOT_EXP {
        ret[0][i] = (ret[0][i - 1] * ret[0][i - 1]) % MOD;
        ret[1][i] = (ret[1][i - 1] * ret[1][i - 1]) % MOD;

        i += 1;
    }
    ret
};

const fn powmod(mut x: i64, mut exp: i64) -> i64 {
    if exp == 0 {
        return 1;
    }
    let mut r = 1;
    x %= MOD;
    while exp > 1 {
        if exp & 1 == 1 {
            r = (r * x) % MOD;
        }
        x = (x * x) % MOD;
        exp >>= 1;
    }
    (r * x) % MOD
}
const fn modinv(x: i64) -> i64 {
    powmod(x, MOD - 2)
}
pub fn ntt(a: &mut [i32], invert: bool) {
    let n = a.len();
    let mut j = 0;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;

        if i < j {
            a.swap(i, j);
        }
    }
    let mut len = 2;
    let mut len2 = 1;
    let mut pow = ROOT_EXP;
    while len <= n {
        let wlen = ROOT_POWS[usize::from(invert)][pow];
        for i in (0..n).step_by(len) {
            let mut w = 1;
            for j in 0..len2 {
                let u = a[i | j];
                let v = ((a[i | j | len2] as i64 * w) % MOD) as i32;
                a[i | j] = if u + v < MOD as i32 {
                    u + v
                } else {
                    u + v - MOD as i32
                };
                a[i | j | len2] = if u >= v { u - v } else { u + MOD as i32 - v };
                w = (w * wlen) % MOD;
            }
        }
        len <<= 1;
        len2 <<= 1;
        pow -= 1;
    }

    if invert {
        let n_1 = modinv(n as i64);

        for x in a.iter_mut() {
            *x = ((*x as i64 * n_1) % MOD) as i32;
        }
    }
}
#[must_use]
pub fn polymul(a: &[i32], b: &[i32]) -> Vec<i32> {
    let mut fa = a.iter().copied().collect_vec();
    let mut fb = b.iter().copied().collect_vec();
    let Some(trunc_a) = a.iter().rposition(|&e| e != 0) else {
        return vec![];
    };
    let Some(trunc_b) = b.iter().rposition(|&e| e != 0) else {
        return vec![];
    };

    let n = (trunc_a + 1 + trunc_b + 1 - 1).next_power_of_two();
    fa.resize(n, 0);
    fb.resize(n, 0);

    ntt(&mut fa, false);
    ntt(&mut fb, false);

    for i in 0..n {
        fa[i] = ((i64::from(fa[i]) * i64::from(fb[i])) % MOD) as i32;
    }
    ntt(&mut fa, true);
    //fa.truncate(unsafe { fa.iter().rposition(|&e| e != 0).unwrap_unchecked() });
    fa
}
*/
/*
pub fn polypow(a: &[i32], mut exp: u32) -> Vec<i32> {
    if exp == 0 {
        return vec![1];
    }
    let mut r_store = vec![1];
    let mut fa = a.iter().copied().collect_vec();
    let mut tmp_store = vec![];
    let (mut x, mut r, mut tmp) = (&mut fa, &mut r_store, &mut tmp_store);

    while exp > 1 {
        if exp & 1 == 1 {
            polymul(r, x, tmp);
            std::mem::swap(&mut tmp, &mut r);
            //r = (r * x) ;
        }
        polymul(x, x, tmp);
        std::mem::swap(&mut tmp, &mut x);
        //x = (x * x);
        exp >>= 1;
    }
    // r * x
    polymul(r, x, tmp)
}
 */
