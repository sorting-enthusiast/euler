use itertools::Itertools;

use crate::{
    divisor_sums,
    p300_399::e362::mult_sparse,
    utils::{FIArray::FIArray, math::iroot, primes::wheel_sieve},
};
use std::ops::{Add, AddAssign, Sub, SubAssign};
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Pair(pub usize, pub usize);
impl Add for Pair {
    type Output = Pair;
    fn add(self, rhs: Self) -> Self::Output {
        Self(self.0 + rhs.0, self.1 + rhs.1)
    }
}
impl AddAssign for Pair {
    fn add_assign(&mut self, rhs: Self) {
        self.0 += rhs.0;
        self.1 += rhs.1;
    }
}
impl Sub for Pair {
    type Output = Pair;
    fn sub(self, rhs: Self) -> Self::Output {
        Self(self.0 - rhs.0, self.1 - rhs.1)
    }
}
impl SubAssign for Pair {
    fn sub_assign(&mut self, rhs: Self) {
        self.0 -= rhs.0;
        self.1 -= rhs.1;
    }
}
impl Add<usize> for Pair {
    type Output = Pair;
    fn add(self, rhs: usize) -> Self::Output {
        Self(self.0 + rhs, self.1 + rhs)
    }
}
impl AddAssign<usize> for Pair {
    fn add_assign(&mut self, rhs: usize) {
        self.0 += rhs;
        self.1 += rhs;
    }
}
impl Sub<usize> for Pair {
    type Output = Pair;
    fn sub(self, rhs: usize) -> Self::Output {
        Self(self.0 - rhs, self.1 - rhs)
    }
}
impl SubAssign<usize> for Pair {
    fn sub_assign(&mut self, rhs: usize) {
        self.0 -= rhs;
        self.1 -= rhs;
    }
}

crate::incremental_flattening::DynamicPrefixSum_impl_for!(Pair);
// 1e16: 477.575608s
// 1e15: 118.9152936s
// 1e14: 28.5571859s
pub fn log_zeta_fast_partial_flatten_v2(n: usize) -> FIArray {
    fn mult_correction(d: &FIArray, primes: &[u64]) -> FIArray {
        struct Correction(FIArray, usize);
        impl Correction {
            fn fill(&mut self, primes: &[u64], lim: usize, x: usize, y: usize) {
                self.0[x] += y;
                self.1 += 1;
                for (i, &p) in primes.iter().enumerate() {
                    let p = p as usize;
                    if p > lim / p {
                        break;
                    }
                    let mut pp = p * p;
                    let mut new_lim = lim / pp;
                    for e in 2.. {
                        let hp = 1 << (e - 2);
                        if hp != 0 {
                            self.fill(&primes[i + 1..], new_lim, x * pp, y * hp);
                        }
                        if p > new_lim {
                            break;
                        }
                        pp *= p;
                        new_lim /= p;
                    }
                }
            }
        }
        let mut correction = Correction(FIArray::new(d.x), 0);
        correction.fill(primes, d.x, 1, 1);
        for i in 1..correction.0.arr.len() {
            correction.0.arr[i] += correction.0.arr[i - 1];
        }
        //dbg!(correction.1);
        mult_sparse(d, &correction.0)
    }
    const fn inv_odd(mut k: usize) -> usize {
        let mut exp = (1u64 << 63) - 1;

        let mut r: usize = 1;
        while exp > 1 {
            r = r.overflowing_mul(k).0;
            k = k.overflowing_mul(k).0;
            exp >>= 1;
        }
        r.overflowing_mul(k).0
    }

    const INVS: [usize; 4] = [0, 2, 1, inv_odd(3) << 1];
    let start = std::time::Instant::now();
    let x = iroot::<4>(n) + 1;
    let primes = wheel_sieve(x as u64 - 1);
    let mut zeta_2 = divisor_sums(n);
    dbg!(start.elapsed());
    zeta_2 = mult_correction(&zeta_2, &primes);
    dbg!(start.elapsed());
    let rt = zeta_2.isqrt;
    let len = zeta_2.arr.len();
    let pairs = FIArray::keys(n)
        .zip(zeta_2.arr)
        .map(|(a, b)| Pair(a, b))
        .collect_vec()
        .into_boxed_slice();
    let mut zeta_fenwick = DynamicPrefixSumPair(pairs, len);

    zeta_fenwick.shrink_flattened_prefix(1);
    zeta_fenwick.dec(0);
    let get_index = |v| -> usize {
        if v == 0 {
            unsafe { core::hint::unreachable_unchecked() };
        } else if v <= rt {
            v - 1
        } else {
            len - (n / v)
        }
    };

    for p in primes {
        let p = p as usize;
        zeta_fenwick.extend_flattened_prefix(1 + get_index(p * p));

        let lim = n / p;
        let mut j = 1;
        let mut cur = zeta_fenwick.sum(get_index(lim));
        while (j + 1) <= lim / (j + 1) {
            let next = zeta_fenwick.sum(get_index(lim / (j + 1)));
            if next != cur {
                zeta_fenwick.sub(len - j, Pair(cur.0 - next.0, 2 * (cur.1 - next.1)));
                cur = next;
            }
            j += 1;
        }
        for i in (p..=lim / j).rev() {
            let next = zeta_fenwick.sum(i - 2);
            if next != cur {
                zeta_fenwick.sub(get_index(p * i), Pair(cur.0 - next.0, 2 * (cur.1 - next.1)));
                cur = next;
            }
        }
    }
    let mut zeta_1_2 = zeta_fenwick.flatten();
    dbg!(start.elapsed());

    let mut ret = FIArray::new(n);
    for i in (1..len).rev() {
        zeta_1_2[i] -= zeta_1_2[i - 1];
    }
    for i in 0..x - 1 {
        ret.arr[i] = zeta_1_2[i].0;
        zeta_1_2[i] = Pair(0, 0);
    }
    for i in x..len {
        zeta_1_2[i] += zeta_1_2[i - 1];
    }
    for i in x..len {
        zeta_1_2[i].1 -= 2 * zeta_1_2[i].0;
    }
    dbg!(start.elapsed());

    // zeta now equals zeta_t - 1, and zeta_2 (zeta_t - 1)^2
    // compute log(zeta_t) using log(x + 1) = x^3 / 3 - x^2 / 2 + x
    // in order to not have to deal with rational numbers, we compute 2 * log(zeta_t)
    // and adjust later

    for i in x..=len {
        ret.arr[i - 1] = zeta_1_2[i - 1].0 * INVS[1];
    }

    let pa = {
        let mut vec = vec![];
        for i in x..=rt {
            if zeta_1_2[i - 1].0 != zeta_1_2[i - 2].0 {
                vec.push(i);
            }
        }
        vec.push(rt + 1);
        vec
    };
    let va = &pa[..pa.len() - 1];

    /* let ind = zeta_2.get_index(x.pow(2));
    assert!(pow_zeta.arr[..ind].iter().all(|&e| e == 0));
    dbg!(ind, x.pow(2), rt - 1, len, len - ind); */

    for i in rt + 1..=len {
        ret.arr[i - 1] -= zeta_1_2[i - 1].1 * INVS[2];
    }

    //pow_zeta = mult_sparse(&zeta, &pow_zeta);

    //zeta.arr[rt..].fill(0);
    for i in rt..len {
        zeta_1_2[i].0 = 0;
    }
    for &i in va {
        for j in 1..=rt / i {
            //zeta[n / j] += pow_zeta[n / (i * j)];
            zeta_1_2[len - j].0 += zeta_1_2[len - i * j].1;
        }
    }
    //zeta.arr[..rt].fill(0);

    let ind = get_index(x.pow(3));
    //dbg!(ind, len, len - ind);
    //assert!(pow_zeta.arr[..ind].iter().all(|&e| e == 0)); */
    for i in ind + 1..=len {
        ret.arr[i - 1] += zeta_1_2[i - 1].0 * INVS[3];
    }
    //dbg!(start.elapsed());
    // correction phase: get rid of contributions of prime powers
    for i in (x..len).rev() {
        ret.arr[i] -= ret.arr[i - 1];
    }

    for x in x..=rt {
        let v = ret.arr[x - 1] / INVS[1];
        if v == 0 {
            continue;
        }
        let mut e = 1;
        let mut px = x;
        while px <= n / x {
            e += 1;
            px *= x;

            ret[px] -= INVS[e];
        }
    }
    for i in 1..x - 1 {
        ret.arr[i] += ret.arr[i - 1];
    }
    for i in x - 1..len {
        ret.arr[i] /= INVS[1];
        ret.arr[i] += ret.arr[i - 1];
    }
    ret
}
