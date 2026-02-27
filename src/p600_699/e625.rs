use crate::utils::{
    FIArray::{DirichletFenwickI128, FIArrayI128},
    math::iroot,
    multiplicative_function_summation::{sum_n_i64, sum_n_i128, totient_sum},
};

// sum of convolution of N and totient
// dirichlet hyperbola method
pub fn main() {
    const N: usize = 1e11 as _;
    const MOD: i64 = 998_244_353;
    let start = std::time::Instant::now();

    let sqrtn = N.isqrt();
    let totient_sums = totient_sum::<MOD>(N);
    let mut sum = (totient_sums[N] + sum_n_i64::<MOD>(N)) % MOD;
    for i in 2..=sqrtn {
        sum += (i as i64 * totient_sums[N / i]) % MOD;
        sum += ((totient_sums.arr[i - 1] - totient_sums.arr[i - 2]) % MOD
            * sum_n_i64::<MOD>(N / i))
            % MOD;
        sum %= MOD;
    }
    sum -= (totient_sums[sqrtn] * sum_n_i64::<MOD>(sqrtn)) % MOD;
    if sum < 0 {
        sum += MOD;
    }
    let end = start.elapsed();
    println!("res = {sum}, took {end:?}");
}
// 131.5644394s
pub fn solve_ext() {
    const N: usize = 1 << 50; // 3812934143065599435357584102851032

    let start = std::time::Instant::now();
    let mut id = DirichletFenwickI128::from(FIArrayI128::id::<0>(N));
    let lim = iroot::<6>(2 * N) + 1;
    let mut primes = vec![];
    println!("started removal of primes < {lim}: {:?}", start.elapsed());
    for p in 2..lim {
        if id.get_bucket_prefix(p - 1) == 1 {
            continue;
        }
        primes.push(p);
        id.sparse_mul_at_most_one(p, p as _);
    }
    let mut id_lim = FIArrayI128::from(id);

    println!(
        "Finished removal of primes < {lim}, started convolution: {:?}",
        start.elapsed()
    );
    let id2_lim = mult_i128(&id_lim, &id_lim);
    println!(
        "finished convolution, started removing primes < {lim}: {:?}",
        start.elapsed()
    );
    for (i, v) in FIArrayI128::keys(N).enumerate() {
        id_lim.arr[i] = v as i128;
    }
    let mut zeta = DirichletFenwickI128::from(id_lim);
    for &p in &primes {
        zeta.sparse_mul_at_most_one(p, 1);
    }
    let zeta_lim = FIArrayI128::from(zeta);
    println!(
        "Finished removal of primes < {lim}, started division: {:?}",
        start.elapsed()
    );
    let mut approx = DirichletFenwickI128::from(div_i128(&id2_lim, &zeta_lim));
    drop(id2_lim);
    drop(zeta_lim);
    println!(
        "Finished division, started adding back primes < {lim}: {:?}",
        start.elapsed()
    );
    for &p in primes.iter().rev() {
        approx.sparse_mul_unlimited(p, 2 * p as i128 - 1);
    }
    let approx = FIArrayI128::from(approx);
    println!(
        "Finished adding back primes < {lim}, started correction: {:?}",
        start.elapsed()
    );
    let res = mult_correction(&approx, &primes, |pp, p, e| {
        pp as i128 * (i128::from(e) + 1) - i128::from(e) * (pp / p) as i128
    })[N];
    let end = start.elapsed();
    println!("res = {res}, took {end:?}");
}
pub fn solve_ext_alt() {
    const N: usize = 1 << 50; // 3812934143065599435357584102851032

    let start = std::time::Instant::now();
    let mut id = DirichletFenwickI128::from(FIArrayI128::id::<0>(N));
    let mut zeta = DirichletFenwickI128::zeta(N);
    let lim = iroot::<6>(2 * N) + 1;
    let mut primes = vec![];
    println!("started removal of primes < {lim}: {:?}", start.elapsed());
    for p in 2..lim {
        if zeta.get_bucket_prefix(p - 1) == 1 {
            continue;
        }
        primes.push(p);
        zeta.sparse_mul_at_most_one(p, 1);
        id.sparse_mul_at_most_one(p, p as _);
    }
    let id_lim = FIArrayI128::from(id);
    let zeta_lim = FIArrayI128::from(zeta);

    println!(
        "Finished removal of primes < {lim}, started division: {:?}",
        start.elapsed()
    );
    let mut approx = DirichletFenwickI128::from(div_i128(&id_lim, &zeta_lim));
    drop(id_lim);
    drop(zeta_lim);
    println!(
        "Finished division, started adding back primes < {lim}: {:?}",
        start.elapsed()
    );
    for &p in primes.iter().rev() {
        approx.sparse_mul_unlimited(p, p as i128 - 1);
    }
    let approx = FIArrayI128::from(approx);
    println!(
        "Finished adding back primes < {lim}, started correction: {:?}",
        start.elapsed()
    );
    let totient = mult_correction(&approx, &primes, |pp, p, _e| pp as i128 - (pp / p) as i128);
    let rt_n = totient.isqrt;
    let len = totient.arr.len();
    let mut res =
        totient.arr[len - 1] + sum_n_i128::<0>(N) - totient.arr[rt_n - 1] * sum_n_i128::<0>(rt_n);
    for i in 2..=rt_n {
        res += i as i128 * totient.arr[len - i]
            + (totient.arr[i - 1] - totient.arr[i - 2]) * sum_n_i128::<0>(N / i);
    }
    let end = start.elapsed();
    println!("res = {res}, took {end:?}");
}

#[must_use]
pub fn div_i128(a: &FIArrayI128, b: &FIArrayI128) -> FIArrayI128 {
    unsafe { core::hint::assert_unchecked(a.x == b.x) };
    let R2 = a.isqrt;
    let n = a.x;
    let mut res = FIArrayI128::new(n);

    let s1 = |ds: &FIArrayI128| {
        let mut vec = vec![];
        if ds.arr[0] != 0 {
            vec.push((1, ds.arr[0]));
        }
        for i in 1..R2 {
            if ds.arr[i] != ds.arr[i - 1] {
                vec.push((i + 1, ds.arr[i] - ds.arr[i - 1]));
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
        res.arr[i] = a.arr[i] - a.arr[i - 1];
    }
    let mut sum = 0;
    for i in 1..=R2 {
        let val = res.arr[i - 1];
        sum += val;
        res.arr[i - 1] = sum;
        if val == 0 {
            continue;
        }
        pa.push((i, val));
        for (y, fy) in &vb[1..] {
            if y * i > R2 {
                break;
            }
            res.arr[i * y - 1] -= val * fy;
        }
    }
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
            for (y, fy) in &pb[l.max(l0)..r] {
                res.arr[len - Nx / y] += fx * fy;
            }
        }
        r = r.max(l);

        if r > 0 && pb[r].0 <= Nx {
            res.arr[len - Nx / pb[r].0] -= fx * b.arr[pb[r - 1].0 - 1];
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
            for (x, fx) in &va[l.max(l0)..r] {
                res.arr[len - Ny / x] += fy * fx;
            }
        }

        r = r.max(l);

        if r > 0 && pa[r].0 <= Ny {
            res.arr[len - Ny / pa[r].0] -= fy * res.arr[pa[r - 1].0 - 1];
        }
        bound_z = Ny / pa[r].0;
    }

    res.arr[R2] += a.arr[R2 - 1];
    for i in R2 + 1..len {
        res.arr[i] += res.arr[i - 1];
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
            res.arr[len - i] += fx * b.arr[Nx / i - 1];
            i -= 1;
        }
        while i > 0 {
            res.arr[len - i] += fx * b.arr[len - x * i];
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
        let mut ans = a.arr[len - i] - res.arr[len - i];
        for (y, fy) in &pb[1..c_y] {
            ans -= fy * res[Nz / y];
        }
        res.arr[len - i] = ans;
    }
    res
}
#[must_use]
pub fn mult_i128(a: &FIArrayI128, b: &FIArrayI128) -> FIArrayI128 {
    unsafe { core::hint::assert_unchecked(a.x == b.x) };
    let R2 = a.isqrt as usize;
    let n = a.x as usize;
    let mut res = FIArrayI128::new(n as _);

    let s1 = |ds: &FIArrayI128| {
        let mut vec = vec![];
        if ds.arr[0] != 0 {
            vec.push((1, ds.arr[0]));
        }
        for i in 1..R2 {
            if ds.arr[i] != ds.arr[i - 1] {
                vec.push((i + 1, ds.arr[i] - ds.arr[i - 1]));
            }
        }
        vec.push((R2 + 1, 0));
        vec
    };
    let len = res.arr.len();
    if a == b {
        let pa = s1(a);
        let va = &pa[..pa.len() - 1];
        let mut r = va.len();
        let mut l = 0;
        for &(x, fx) in va {
            res[x * x] += fx * fx;

            l += 1;
            let Nx = n / x;
            while r > l && Nx / pa[r - 1].0 < r - 1 {
                r -= 1;
            }
            if r < l {
                r = l;
            }
            let X = R2 / x;
            let Nx = n / x;

            let mut i = l;

            while i != r && pa[i].0 <= X {
                res.arr[x * pa[i].0 - 1] += 2 * fx * pa[i].1;
                i += 1;
            }
            while i != r {
                res.arr[len - Nx / pa[i].0] += 2 * fx * pa[i].1;
                i += 1;
            }
            if r != 0 && pa[r].0 <= Nx {
                res[x * pa[r].0] -= fx * a.arr[pa[r - 1].0 - 1] * 2;
            }
        }
        for i in 1..len {
            res.arr[i] += res.arr[i - 1];
        }
        r = va.len();
        l = 0;
        for &(x, fx) in va {
            l += 1;
            let Nx = n / x;
            while r > l && Nx / pa[r - 1].0 < r - 1 {
                r -= 1;
            }
            if r < l {
                r = l;
            }

            let mut i = Nx / pa[r].0;
            let X = R2 / x;
            while i > X {
                res.arr[len - i] += 2 * fx * a.arr[Nx / i - 1];
                i -= 1;
            }
            while i > 0 {
                res.arr[len - i] += 2 * fx * a.arr[len - x * i];
                i -= 1;
            }
        }
    } else {
        let pa = s1(a);
        let va = &pa[..pa.len() - 1];
        let pb = s1(b);
        let vb = &pb[..pb.len() - 1];
        res.arr[0] += a.arr[0] * b.arr[0];
        for i in 2..=R2 {
            res[(i * i) as _] += (a.arr[i - 1] - a.arr[i - 2]) * (b.arr[i - 1] - b.arr[i - 2]);
        }
        let mut r = vb.len();
        let mut l = 0;
        for &(x, y) in va {
            while pb[l].0 <= x {
                l += 1;
            }
            let Nx = n / x;
            while r > l && Nx / pb[r - 1].0 < r - 1 {
                r -= 1;
            }
            if r < l {
                r = l;
            }
            let X = R2 / x;
            let Nx = n / x;

            let mut i = l;

            while i != r && pb[i].0 <= X {
                res.arr[x * pb[i].0 - 1] += y * pb[i].1;
                i += 1;
            }
            while i != r {
                res.arr[len - Nx / pb[i].0] += y * pb[i].1;
                i += 1;
            }
            if r != 0 && pb[r].0 <= Nx {
                res[(x * pb[r].0) as _] -= y * b.arr[pb[r - 1].0 - 1];
            }
        }
        r = va.len();
        l = 0;
        for &(x, y) in vb {
            while pa[l].0 <= x {
                l += 1;
            }
            let Nx = n / x;
            while r > l && Nx / pa[r - 1].0 < r - 1 {
                r -= 1;
            }
            if r < l {
                r = l;
            }
            let X = R2 / x;
            let Nx = n / x;

            let mut i = l;

            while i != r && pa[i].0 <= X {
                res.arr[x * pa[i].0 - 1] += y * pa[i].1;
                i += 1;
            }
            while i != r {
                res.arr[len - Nx / pa[i].0] += y * pa[i].1;
                i += 1;
            }
            if r != 0 && pa[r].0 <= Nx {
                res[(x * pa[r].0) as _] -= y * a.arr[pa[r - 1].0 - 1];
            }
        }
        for i in 1..len {
            res.arr[i] += res.arr[i - 1];
        }
        r = vb.len();
        l = 0;
        for &(x, y) in va {
            while pb[l].0 <= x {
                l += 1;
            }
            let Nx = n / x;
            while r > l && Nx / pb[r - 1].0 < r - 1 {
                r -= 1;
            }
            if r < l {
                r = l;
            }

            let mut i = Nx / pb[r].0;
            let X = R2 / x;
            while i > X {
                res.arr[len - i] += y * b.arr[Nx / i - 1];
                i -= 1;
            }
            while i > 0 {
                res.arr[len - i] += y * b.arr[len - x * i];
                i -= 1;
            }
        }
        r = va.len();
        l = 0;
        for &(x, y) in vb {
            while pa[l].0 <= x {
                l += 1;
            }
            let Nx = n / x;
            while r > l && Nx / pa[r - 1].0 < r - 1 {
                r -= 1;
            }
            if r < l {
                r = l;
            }

            let mut i = Nx / pa[r].0;
            let X = R2 / x;
            while i > X {
                res.arr[len - i] += y * a.arr[Nx / i - 1];
                i -= 1;
            }
            while i > 0 {
                res.arr[len - i] += y * a.arr[len - x * i];
                i -= 1;
            }
        }
    }
    res
}

#[must_use]
pub fn mult_correction_single(
    d: &FIArrayI128,
    primes: &[usize],
    f: impl Fn(usize, usize, u8) -> i128,
) -> i128 {
    fn dfs(
        d: &FIArrayI128,
        primes: &[usize],
        lim: usize,
        x: usize,
        hx: i128,
        f: &impl Fn(usize, usize, u8) -> i128,
    ) -> i128 {
        let mut res = hx * d[d.x / x];
        for (i, &p) in primes.iter().enumerate() {
            if p > lim / p {
                break;
            }
            let fp = f(p, p, 1);
            let mut prev = fp;
            let mut pp = p * p;
            let mut new_lim = lim / pp;
            for e in 2.. {
                let cur = f(pp, p, e);
                let hp = cur - fp * prev;
                if hp != 0 {
                    res += dfs(d, &primes[i + 1..], new_lim, x * pp, hx * hp, f);
                }
                prev = cur;
                if p > new_lim {
                    break;
                }
                pp *= p;
                new_lim /= p;
            }
        }
        res
    }

    dfs(d, primes, d.x, 1, 1, &f)
}

#[must_use]
pub fn mult_correction(
    d: &FIArrayI128,
    primes: &[usize],
    f: impl Fn(usize, usize, u8) -> i128,
) -> FIArrayI128 {
    struct Correction(FIArrayI128, usize);
    impl Correction {
        fn fill(
            &mut self,
            primes: &[usize],
            lim: usize,
            x: usize,
            y: i128,
            f: &impl Fn(usize, usize, u8) -> i128,
        ) {
            self.0[x] += y;
            self.1 += 1;
            for (i, &p) in primes.iter().enumerate() {
                if p > lim / p {
                    break;
                }
                let fp = f(p, p, 1);
                let mut prev = fp;
                let mut pp = p * p;
                let mut new_lim = lim / pp;
                for e in 2.. {
                    let cur = f(pp, p, e);
                    let hp = cur - fp * prev;
                    if hp != 0 {
                        self.fill(&primes[i + 1..], new_lim, x * pp, y * hp, f);
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
    let mut correction = Correction(FIArrayI128::new(d.x), 0);
    correction.fill(primes, d.x, 1, 1, &f);
    for i in 1..correction.0.arr.len() {
        correction.0.arr[i] += correction.0.arr[i - 1];
    }
    dbg!(correction.1);
    mult_sparse_i128(d, &correction.0)
}

pub fn mult_sparse_with_buffer_i128(a: &FIArrayI128, b: &FIArrayI128, res: &mut FIArrayI128) {
    unsafe { core::hint::assert_unchecked(a.x == b.x && a.x == res.x) };
    res.arr.fill(0);
    let R2 = a.isqrt;
    let n = a.x;
    let s1 = |ds: &FIArrayI128| {
        let mut vec = vec![];
        if ds.arr[0] != 0 {
            vec.push((1, ds.arr[0]));
        }
        for i in 1..R2 {
            if ds.arr[i] != ds.arr[i - 1] {
                vec.push((i + 1, ds.arr[i] - ds.arr[i - 1]));
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
            res.arr[x * pa[i].0 - 1] += y * pa[i].1;
            i += 1;
        }
        while i != r {
            res.arr[len - Nx / pa[i].0] += y * pa[i].1;
            i += 1;
        }
        if r != 0 && pa[r].0 <= Nx {
            res[(x * pa[r].0) as _] -= y * a.arr[pa[r - 1].0 - 1];
        }
    }
    for &(x, y) in va {
        res.arr[len - R2 / x] -= y * b.arr[R2 - 1];
    }
    for i in 1..len {
        res.arr[i] += res.arr[i - 1];
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
            res.arr[len - i] += y * a.arr[Nx / i - 1];
            i -= 1;
        }
        while i > 0 {
            res.arr[len - i] += y * a.arr[len - x * i];
            i -= 1;
        }
    }
    for &(x, y) in va {
        for j in 1..=R2 / x {
            res.arr[len - j] += y * b.arr[len - x * j];
        }
    }
}
#[must_use]
pub fn mult_sparse_i128(a: &FIArrayI128, b: &FIArrayI128) -> FIArrayI128 {
    let mut res = a.clone();
    mult_sparse_with_buffer_i128(a, b, &mut res);
    res
}
