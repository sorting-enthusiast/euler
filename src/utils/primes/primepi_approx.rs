use crate::utils::multiplicative_function_summation::mobius_sieve;

/// Calculate the logarithmic integral using
/// Ramanujan's formula:
/// <https://en.wikipedia.org/wiki/Logarithmic_integral_function#Series_representation>
///
#[must_use]
pub fn li(x: f64) -> f64 {
    const GAMMA: f64 = 0.577215664901532860606512090082402431_f64;
    if x <= 1. {
        return 0.;
    }
    let mut sum = 0.;
    let mut inner_sum = 0.;
    let mut factorial = 1.;
    let mut p = -1.;
    let mut q;
    let mut power2 = 1.;
    let logx = x.ln();
    let mut k = 0;

    // The condition n < ITERS is required in case the computation
    // does not converge. This happened on Linux i386 where
    // the precision of the libc math functions is very limited.
    for n in 1..1000 {
        p *= -logx;
        factorial *= f64::from(n);
        q = factorial * power2;
        power2 *= 2.;

        while k <= (n - 1) / 2 {
            inner_sum += f64::from(2 * k + 1).recip();
            k += 1;
        }

        let old_sum = sum;
        sum += (p / q) * inner_sum;

        // Not converging anymore
        if (sum - old_sum).abs() <= f64::EPSILON {
            break;
        }
    }

    GAMMA + logx.ln() + x.sqrt() * sum
}

/// Calculate the Eulerian logarithmic integral which is a very
/// accurate approximation of the number of primes <= x.
/// Li(x) > pi(x) for 24 <= x <= ~ 10^316
///
#[must_use]
pub fn Li(x: usize) -> usize {
    const li2: f64 = 1.045163780117492784844588889194613136_f64;

    if x <= 2 {
        0
    } else {
        (li(x as f64) - li2) as usize
    }
}

#[must_use]
pub fn R(x: usize) -> usize {
    /* let Li = |x: f64| {
        const li2: f64 = 1.045163780117492784844588889194613136_f64;
        if x <= 2. { 0. } else { li(x as f64) - li2 }
    }; */
    let x = x as f64;
    let mobius = mobius_sieve(50);
    let mut sum = 0.;
    for i in 1..50 {
        let inv_i = (i as f64).recip();
        let li = li(x.powf(inv_i));
        if li == 0. {
            break;
        }
        sum += f64::from(mobius[i]) * inv_i * li;
    }
    sum as usize
}
