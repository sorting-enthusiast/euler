use crate::{
    inverse_pseudo_euler_transform_fraction_i64, mult_correction_single,
    pseudo_euler_transform_fraction_i64, utils::FIArray::FIArrayI64,
};

const N: usize = 1e14 as _;
// sum of coefficients of \frac{\zeta(s)\beta(s)}{\zeta(2s)\beta(2s)}
// TODO: optimize to better time complexity than O\left(n^{2/3}\log^{-4/3}n\right)
// can easily be done in \tilde{O}(\sqrt n) time, just need to implement a linear sieve for the dirichlet inverse of u * \chi_4
pub fn main() {
    original_approach();
}
fn original_approach() {
    let start = std::time::Instant::now();

    let mut chi4 = FIArrayI64::new(N);
    for (i, v) in FIArrayI64::keys(N).enumerate() {
        chi4.arr[i] = [0, 1, 1, 0][v & 3];
    }

    println!("Started first inverse transform: {:?}", start.elapsed());
    let chi4_p = inverse_pseudo_euler_transform_fraction_i64(chi4);

    println!("Started second inverse transform: {:?}", start.elapsed());
    let pi = inverse_pseudo_euler_transform_fraction_i64(FIArrayI64::unit(N));

    let mut primes = vec![];
    for p in 2..=pi.isqrt {
        if pi.arr[p - 1] != pi.arr[p - 2] {
            primes.push(p);
        }
    }

    let mut sum_over_primes = pi;
    for i in 0..sum_over_primes.arr.len() {
        sum_over_primes.arr[i] += chi4_p.arr[i];
    }

    println!("Started transform: {:?}", start.elapsed());
    let approx = pseudo_euler_transform_fraction_i64(sum_over_primes);

    println!("Started correction: {:?}", start.elapsed());
    let res = mult_correction_single(&approx, &primes, |_, p, e| {
        if e > 2 {
            return 0;
        }
        match p & 3 {
            1 => 1 + i64::from(e == 1),
            2 => i64::from(e == 1),
            3 => i64::from(e == 2),
            _ => unsafe { core::hint::unreachable_unchecked() },
        }
    });

    println!("res = {res}, took {:?}", start.elapsed());
}
