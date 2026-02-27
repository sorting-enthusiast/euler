use crate::utils::{
    FIArray::{DirichletFenwickU64Mod, FIArrayI64, FIArrayU64},
    math::iroot,
    multiplicative_function_summation::{sum_n_u64, totient_sum},
};
use std::{collections::HashMap, time::Instant};

const MOD: i64 = 1e9 as i64;
const N: usize = 510_510;
const PHI_N: usize = 92160;

const M: usize = 1e11 as _;

const PRIMES: [usize; 7] = [2, 3, 5, 7, 11, 13, 17];
const fn lookup_create() -> [(usize, usize); 128] {
    let mut ret = [(1, 1); 128];
    let mut i = 1_usize;
    while i < 128 {
        let mut b = i;
        while b != 0 {
            let p = b.trailing_zeros() as usize;
            ret[i].0 *= PRIMES[p];
            ret[i].1 *= PRIMES[p] - 1;
            b &= b - 1;
        }
        i += 1;
    }
    ret
}
const LOOKUP: [(usize, usize); 128] = lookup_create();
const fn mobius(pmask: i8) -> i64 {
    1 - ((pmask.count_ones() as i64 & 1) << 1)
}

fn original_s(
    n: i8,
    m: usize,
    totient_sums: &FIArrayI64,
    cache: &mut HashMap<(i8, usize), i64>,
) -> i64 {
    let phi_n = LOOKUP[n as usize].1;
    if m == 1 {
        return phi_n as _;
    }

    if let Some(&v) = cache.get(&(n, m)) {
        return v;
    }
    assert_ne!(n, 0); // only times when n would be 0 are unrolled

    let mut sum = 0;
    let mut d_submask = n;
    while d_submask != 0 {
        let mut inner = 0;
        let complement = n & !d_submask;
        let mut d_prime_submask = complement;
        while d_prime_submask != 0 {
            let ndp = d_prime_submask | d_submask;
            let prod = LOOKUP[ndp as usize].0;
            if m >= prod {
                inner += mobius(d_prime_submask) * original_s(ndp, m / prod, totient_sums, cache);
                inner %= MOD;
                if inner < 0 {
                    inner += MOD;
                }
            }
            d_prime_submask = (d_prime_submask - 1) & complement;
        }

        let (d, phi_d) = LOOKUP[d_submask as usize];

        if m >= d {
            inner += original_s(d_submask, m / d, totient_sums, cache);
            inner %= MOD;
        }

        let coeff = (phi_n * d) / phi_d;
        sum += coeff as i64 * inner;
        sum %= MOD;
        d_submask = (d_submask - 1) & n;
    }
    {
        let mut inner = 0;
        let complement = n;
        let mut d_prime_submask = complement;
        while d_prime_submask != 0 {
            let ndp = d_prime_submask;
            let prod = LOOKUP[ndp as usize].0;
            if m >= prod {
                inner += mobius(d_prime_submask) * original_s(ndp, m / prod, totient_sums, cache);
                inner %= MOD;
                if inner < 0 {
                    inner += MOD;
                }
            }
            d_prime_submask = (d_prime_submask - 1) & complement;
        }

        inner += totient_sums[m];
        if inner >= MOD {
            inner -= MOD;
        }

        sum += phi_n as i64 * inner;
        sum %= MOD;
    }
    cache.insert((n, m), sum);
    sum
}

fn s(n: i8, m: usize, totient_sums: &FIArrayI64) -> i64 {
    let phi_n = LOOKUP[n as usize].1;
    if m == 0 {
        0
    } else if m == 1 {
        phi_n as _
    } else {
        assert_ne!(n, 0); // only times when n would be 0 are unrolled
        let mut sum = 0;
        let mut d_submask = n;
        loop {
            let (d, phi_d) = LOOKUP[d_submask as usize];
            sum += ((phi_n / phi_d) as i64 * s(d_submask, m / d, totient_sums)) % MOD;
            sum %= MOD;
            d_submask = (d_submask - 1) & n;
            if d_submask == 0 {
                break;
            }
        }
        sum += phi_n as i64 * totient_sums[m];
        sum % MOD
    }
}

fn s_prime_rec(n: i8, m: usize, totient_sums: &FIArrayI64) -> i64 {
    let phi_n = LOOKUP[n as usize].1;
    if m == 1 {
        phi_n as _
    } else if n == 0 {
        totient_sums[m]
    } else {
        let p = PRIMES[n.trailing_zeros() as usize];
        let mut sum = (p - 1) as i64 * s_prime_rec(n & (n - 1), m, totient_sums);
        if m >= p {
            sum += s_prime_rec(n, m / p, totient_sums);
        }
        sum % MOD
    }
}

pub fn main() {
    let start = Instant::now();
    let sums = totient_sum::<MOD>(M);
    let end = start.elapsed();
    println!("all totient summations precomputed in {end:?}");
    // FIArrayI64::keys(M) = all possible values of M div n, hence all possible totient queries are already in the cache
    // note, that only queries where n is 17-smooth might occur, but inserting only those values is a pain in the ass

    let start = Instant::now();
    let sum = s_prime_rec(0x7f, M, &sums);
    let end = start.elapsed();
    println!("Recursion over primes: {sum}, took {end:?}");

    let start = Instant::now();
    let sum = s(0x7f, M, &sums);
    let end = start.elapsed();
    println!("Recursion over divisors: {sum}, took {end:?}");

    let start = Instant::now();
    let mut cache = HashMap::new();
    let sum = original_s(0x7f, M, &sums, &mut cache);
    let end = start.elapsed();
    println!("Recursion over divisors (original): {sum}, took {end:?}");
    //dbg!(cache.len());
}
pub fn solve_lucy() {
    const MOD: u64 = 1e9 as _;
    fn mult_correction(
        d: &FIArrayU64,
        primes: &[usize],
        f: impl Fn(usize, usize, u8) -> u64,
    ) -> FIArrayU64 {
        struct Correction(FIArrayU64);
        impl Correction {
            fn fill(
                &mut self,
                primes: &[usize],
                lim: usize,
                x: usize,
                y: u64,
                f: &impl Fn(usize, usize, u8) -> u64,
            ) {
                let entry = &mut self.0[x];
                *entry += y as u64;
                if *entry >= MOD {
                    *entry -= MOD;
                }
                for (i, &p) in primes.iter().enumerate() {
                    if p > lim / p {
                        break;
                    }
                    let fp = u64::from(f(p, p, 1));
                    let mut prev = fp;
                    let mut pp = p * p;
                    let mut new_lim = lim / pp;
                    for e in 2.. {
                        let cur = u64::from(f(pp, p, e));
                        let mut hp = cur + MOD - (fp * prev) % MOD;
                        if hp >= MOD {
                            hp -= MOD;
                        }
                        if hp != 0 {
                            self.fill(&primes[i + 1..], new_lim, x * pp, (y * hp) % MOD, f);
                        }
                        prev = cur;
                        if p > new_lim {
                            break;
                        }
                        pp *= p;
                        new_lim /= p;
                    }
                }
            }
        }
        let mut correction = Correction(FIArrayU64::new(d.x));
        correction.fill(primes, d.x, 1, 1, &f);
        for i in 1..correction.0.arr.len() {
            correction.0.arr[i] += correction.0.arr[i - 1];
            if correction.0.arr[i] >= MOD {
                correction.0.arr[i] -= MOD;
            }
        }
        mult_sparse_mod(d, &correction.0)
    }
    const fn add(v: &mut u64, a: u64, b: u64) {
        *v += (((a as u64) * (b as u64)) % MOD) as u64;
        if *v >= MOD {
            *v -= MOD;
        }
    }
    const fn sub(v: &mut u64, a: u64, b: u64) {
        *v += MOD - (((a as u64) * (b as u64)) % MOD) as u64;
        if *v >= MOD {
            *v -= MOD;
        }
    }
    // taken from https://loj.ac/s/1214183
    fn mult_sparse_mod(a: &FIArrayU64, b: &FIArrayU64) -> FIArrayU64 {
        unsafe { core::hint::assert_unchecked(a.x == b.x) };
        let mut res = a.clone();

        res.arr.fill(0);
        let R2 = a.isqrt;
        let n = a.x;
        let s1 = |ds: &FIArrayU64| {
            let mut vec = vec![];
            if ds.arr[0] != 0 {
                vec.push((1, ds.arr[0]));
            }
            for i in 1..R2 {
                if ds.arr[i] != ds.arr[i - 1] {
                    let mut v = ds.arr[i];
                    sub(&mut v, 1, ds.arr[i - 1]);
                    vec.push((i + 1, v));
                }
            }
            vec.push((R2 + 1, 0));
            vec
        };
        let pa = s1(a);
        let pb = s1(b);
        let va = &pa[..pa.len() - 1];
        let vb = &pb[..pb.len() - 1];
        let len = res.arr.len();

        let mut r = va.len();
        for &(x, y) in vb {
            let Nx = n / x;
            while r != 0 && Nx / pa[r - 1].0 < r - 1 {
                r -= 1;
            }
            let mut i = 0;
            let X = R2 / x;
            let Nx = n / x;
            while i != r && pa[i].0 <= X {
                add(&mut res.arr[x * pa[i].0 - 1], y, pa[i].1);
                i += 1;
            }
            while i != r {
                add(&mut res.arr[len - Nx / pa[i].0], y, pa[i].1);
                i += 1;
            }
            if r != 0 && pa[r].0 <= Nx {
                sub(&mut res[x * pa[r].0], y, a.arr[pa[r - 1].0 - 1]);
            }
        }
        for &(x, y) in va {
            sub(&mut res.arr[len - R2 / x], y, b.arr[R2 - 1]);
        }
        for i in 1..len {
            res.arr[i] += res.arr[i - 1];
            if res.arr[i] >= MOD {
                res.arr[i] -= MOD;
            }
        }
        let mut r = va.len();
        for &(x, y) in vb {
            let Nx = n / x;
            while r != 0 && Nx / pa[r - 1].0 < r - 1 {
                r -= 1;
            }
            let mut i = Nx / pa[r].0;
            let X = R2 / x;
            while i > X {
                add(&mut res.arr[len - i], y, a.arr[Nx / i - 1]);
                i -= 1;
            }
            while i > 0 {
                add(&mut res.arr[len - i], y, a.arr[len - x * i]);
                i -= 1;
            }
        }
        for &(x, y) in va {
            for j in 1..=R2 / x {
                add(&mut res.arr[len - j], y, b.arr[len - x * j]);
            }
        }
        res
    }
    fn div_mod(a: &FIArrayU64, b: &FIArrayU64) -> FIArrayU64 {
        unsafe { core::hint::assert_unchecked(a.x == b.x) };
        let R2 = a.isqrt;
        let n = a.x;
        let mut res = FIArrayU64::new(n);

        let s1 = |ds: &FIArrayU64| {
            let mut vec = vec![];
            if ds.arr[0] != 0 {
                vec.push((1, ds.arr[0]));
            }
            for i in 1..R2 {
                if ds.arr[i] != ds.arr[i - 1] {
                    let mut v = ds.arr[i];
                    sub(&mut v, 1, ds.arr[i - 1]);
                    vec.push((i + 1, v));
                }
            }
            vec.push((R2 + 1, 0));
            vec
        };
        let len = res.arr.len();

        let mut pa = vec![];
        let pb = s1(b);
        let vb = &pb[..pb.len() - 1];
        res.arr[0] = a.arr[0];
        for i in 1..R2 {
            res.arr[i] = a.arr[i]; //- a.arr[i - 1];
            sub(&mut res.arr[i], 1, a.arr[i - 1]);
        }
        let mut sum = 0;
        for i in 1..=R2 {
            let val = res.arr[i - 1];
            //sum += val;
            add(&mut sum, 1, val);
            res.arr[i - 1] = sum;
            //dbg!(sum);
            if val == 0 {
                continue;
            }
            pa.push((i, val));
            for &(y, fy) in &vb[1..] {
                if y * i > R2 {
                    break;
                }
                sub(&mut res.arr[i * y - 1], val, fy);
                //res.arr[i * y - 1] -= val * fy;
            }
        }
        //dbg!(&res.arr[..R2]);
        pa.push((R2 + 1, 0));
        let va = &pa[..pa.len() - 1];

        let mut r = vb.len();
        let mut l0 = r;
        let mut l = 0;
        for &(x, fx) in va {
            while pb[l].0 <= x {
                l += 1;
            }
            let Nx = n / x;
            let X = R2 / x;

            while l0 > l && X < pb[l0 - 1].0 {
                l0 -= 1;
            }
            while r > l && Nx < (r - 1) * pb[r - 1].0 {
                r -= 1;
            }

            if l.max(l0) < r {
                for &(y, fy) in &pb[l.max(l0)..r] {
                    //res.arr[len - Nx / y] += fx * fy;
                    add(&mut res.arr[len - Nx / y], fx, fy);
                }
            }
            r = r.max(l);

            if r > 0 && pb[r].0 <= Nx {
                sub(&mut res.arr[len - Nx / pb[r].0], fx, b.arr[pb[r - 1].0 - 1]);
                //res.arr[len - Nx / pb[r].0] -= fx * b.arr[pb[r - 1].0 - 1];
            }
        }
        r = va.len();
        l0 = r;
        l = 0;
        let mut bound_z = n / (R2 + 1);
        for &(y, fy) in vb {
            while pa[l].0 < y {
                l += 1;
            }
            let Ny = n / y;
            let bound_y = Ny / (bound_z + 1);
            while l0 > l && R2 < y * pa[l0 - 1].0 {
                l0 -= 1;
            }
            while r > l && pa[r - 1].0 > bound_y && Ny < (r - 1) * pa[r - 1].0 {
                r -= 1;
            }
            if l.max(l0) < r {
                for &(x, fx) in &va[l.max(l0)..r] {
                    add(&mut res.arr[len - Ny / x], fy, fx);
                    //res.arr[len - Ny / x] += fy * fx;
                }
            }

            r = r.max(l);

            if r > 0 && pa[r].0 <= Ny {
                let val = res.arr[pa[r - 1].0 - 1];
                sub(&mut res.arr[len - Ny / pa[r].0], fy, val);
                //res.arr[len - Ny / pa[r].0] -= fy * res.arr[pa[r - 1].0 - 1];
            }
            bound_z = Ny / pa[r].0;
        }

        add(&mut res.arr[R2], 1, a.arr[R2 - 1]);
        //res.arr[R2] += a.arr[R2 - 1];
        for i in R2 + 1..len {
            //res.arr[i] += res.arr[i - 1];
            let val = res.arr[i - 1];
            add(&mut res.arr[i], 1, val);
        }
        l = 0;
        r = vb.len();
        for &(x, fx) in va {
            while pb[l].0 <= x {
                l += 1;
            }
            let Nx = n / x;
            while r > l && Nx < (r - 1) * pb[r - 1].0 {
                r -= 1;
            }
            //assert!(r >= l);
            //mul_sparse_large(x, Nx, fx, pairs2[std::max(l, r)].first, block2, block_out);
            if r < l {
                r = l;
            }

            let mut i = Nx / pb[r].0;
            let X = R2 / x;
            while i > X {
                //res.arr[len - i] += fx * b.arr[Nx / i - 1];
                add(&mut res.arr[len - i], fx, b.arr[Nx / i - 1]);
                i -= 1;
            }
            while i > 0 {
                //res.arr[len - i] += fx * b.arr[len - x * i];
                add(&mut res.arr[len - i], fx, b.arr[len - x * i]);
                i -= 1;
            }
        }
        let mut c_y = 0;
        l = 0;
        r = va.len();
        bound_z = n / (R2 + 1);
        for i in (1..=bound_z).rev() {
            while bound_z >= i {
                c_y += 1;
                let y = pb[c_y].0;
                while pa[l].0 < y {
                    l += 1;
                }
                let Ny = n / y;
                let bound_y = Ny / (bound_z + 1);
                while r > l && pa[r - 1].0 > bound_y && (r - 1) * (pa[r - 1].0) > Ny {
                    r -= 1;
                }
                r = r.max(l);
                bound_z = Ny / pa[r].0;
            }
            let Nz = n / i;
            let mut ans = a.arr[len - i]; // - res.arr[len - i];
            sub(&mut ans, 1, res.arr[len - i]);
            for &(y, fy) in &pb[1..c_y] {
                //ans -= fy * res[Nz / y];
                sub(&mut ans, fy, res[Nz / y]);
            }
            res.arr[len - i] = ans;
        }
        res
    }

    const LIM: usize = iroot::<8>(M) + 1;

    let start = std::time::Instant::now();
    let mut zeta_mod = FIArrayU64::new(M);
    for (i, v) in FIArrayU64::keys(M).enumerate() {
        zeta_mod.arr[i] = (v % MOD as usize) as u64;
    }
    let mut id_mod = FIArrayU64::new(M);
    for (i, v) in FIArrayU64::keys(M).enumerate() {
        id_mod.arr[i] = sum_n_u64::<MOD>(v) as u64;
    }

    let mut id_mod = DirichletFenwickU64Mod::<MOD>::from(id_mod);
    let mut zeta_mod = DirichletFenwickU64Mod::<MOD>::from(zeta_mod);
    let mut primes = vec![];
    for p in 2..LIM {
        if zeta_mod.get_bucket_prefix(p - 1) - 1 == 0 {
            continue;
        }
        primes.push(p);
        id_mod.sparse_mul_at_most_one(p, p as _);
        zeta_mod.sparse_mul_at_most_one(p, 1);
    }
    let id_lim = FIArrayU64::from(id_mod);
    let zeta_lim = FIArrayU64::from(zeta_mod);

    let mut approx = DirichletFenwickU64Mod::<MOD>::from(div_mod(&id_lim, &zeta_lim));
    for &p in primes.iter().rev() {
        approx.sparse_mul_unlimited(p, (p - usize::from(p > 17)) as _);
    }
    let approx = FIArrayU64::from(approx);

    let accurate = mult_correction(&approx, &primes[7..], |pe, p, _| (pe - pe / p) as u64 % MOD);
    let res = (PHI_N as u64 * accurate[M] as u64) % MOD;
    println!(
        "Sparse division based: res = {res}, took {:?}",
        start.elapsed()
    );
}
