/// Algorithm 3 in <https://arxiv.org/pdf/1206.3369>
const fn S_Q(n: i64, x1: i64, x2: i64) -> i64 {
    if x1 > x2 {
        return 0;
    }
    let mut x = x2;
    let mut S = 0;
    let mut beta = n / (x + 1);
    let mut epsilon = n % (x + 1);
    let mut delta = (n / x) - beta;
    let mut gamma = beta - x * delta;
    while x >= x1 {
        epsilon += gamma;
        if epsilon >= x {
            delta += 1;
            gamma -= x;
            epsilon -= x;
            if epsilon >= x {
                delta += 1;
                gamma -= x;
                epsilon -= x;
                if epsilon >= x {
                    break;
                }
            }
        } else if epsilon < 0 {
            delta -= 1;
            gamma += x;
            epsilon += x;
        }
        gamma += 2 * delta;
        beta += delta;
        S += beta;
        x -= 1;
    }

    epsilon = n % (x + 1);
    delta = (n / x) - beta;
    gamma = beta - x * delta;
    while x >= x1 {
        epsilon += gamma;
        let delta_2 = epsilon / x;
        delta += delta_2;
        epsilon -= x * delta_2;
        gamma += 2 * delta - x * delta_2;
        beta += delta;
        S += beta;
        x -= 1;
    }

    while x >= x1 {
        S += n / x;
        x -= 1;
    }
    S
}

const fn triangle(i: i64) -> i64 {
    if i & 1 == 0 {
        (i >> 1) * (i + 1)
    } else {
        ((i + 1) >> 1) * i
    }
}
const fn ceil_sqrt(n: u128) -> u128 {
    let r = n.isqrt();
    if r * r == n { r } else { r + 1 }
}
const fn ceil_cbrt(n: i64) -> i64 {
    let r = icbrt(n);
    if r * r * r == n { r } else { r + 1 }
}
const C1: i64 = 10;
const C2: i64 = 10;
// Sadly not particularly useful without arbitrary precision arithmetic, but aside from overflow it does indeed work
fn s(n: i64, x_max: i64) -> i64 {
    let mut s = 0;
    let y_min = n / x_max;
    let x_min = x_max.min(C1 * ceil_cbrt(n << 1));
    let mut a_2 = 1;
    let mut x_2 = x_max;
    let mut y_2 = y_min;
    let mut c_2 = a_2 * x_2 + y_2;

    let mut a_1;
    let mut x_4;
    let mut y_4;
    let mut c_4;
    let mut x_5;
    let mut y_5;
    let mut c_5;
    loop {
        a_1 = a_2 + 1;

        x_4 = (n / a_1).isqrt();
        y_4 = n / x_4;
        c_4 = a_1 * x_4 + y_4;

        x_5 = x_4 + 1;
        y_5 = n / x_5;
        c_5 = a_1 * x_5 + y_5;

        if x_4 < x_min {
            break;
        }

        s += triangle(c_4 - c_2 - x_min) - triangle(c_4 - c_2 - x_5) + triangle(c_5 - c_2 - x_5);
        s += S_R(
            n,
            a_1 * x_2 + y_2 - c_5,
            a_2 * x_5 + y_5 - c_2,
            a_1,
            1,
            c_5,
            a_2,
            1,
            c_2,
        );
        a_2 = a_1;
        x_2 = x_4;
        y_2 = y_4;
        c_2 = c_4;
    }

    s += S_Q(n, 1, x_min - 1); // S1 
    s += (x_max - x_min + 1) * y_min + triangle(x_max - x_min); // S2
    s += (x_min..x_2)
        .map(|x| (n / x) - (a_2 * (x_2 - x) + y_2))
        .sum::<i64>(); // S3

    s
}

fn S_R(
    n: i64,
    mut w: i64,
    mut h: i64,
    a_1: i64,
    b_1: i64,
    mut c_1: i64,
    a_2: i64,
    b_2: i64,
    mut c_2: i64,
) -> i64 {
    let mut s = 0;
    let a_3 = a_1 + a_2;
    let b_3 = b_1 + b_2;
    if h > 0 && {
        let H = |u: i64, v: i64| {
            (b_2 * (u + c_1) - b_1 * (v + c_2)) * (a_1 * (v + c_2) - a_2 * (u + c_1))
        };
        H(w, 1)
    } <= n
    {
        s += w;
        c_2 += 1;
        h -= 1;
    }
    if w > 0 && {
        let H = |u: i64, v: i64| {
            (b_2 * (u + c_1) - b_1 * (v + c_2)) * (a_1 * (v + c_2) - a_2 * (u + c_1))
        };
        H(1, h)
    } <= n
    {
        s += h;
        c_1 += 1;
        w -= 1;
    }
    let U_tan = || {
        (((a_1 * b_2 + a_2 * b_1 + 2 * a_1 * b_1).pow(2) as u128 * n as u128) / (a_3 * b_3) as u128)
            .isqrt() as i64
            - c_1
    };
    let V_floor = |u: i64| {
        (((a_1 * b_2 + a_2 * b_1) * (u + c_1)
            - ceil_sqrt(((u + c_1) as u128).pow(2) - 4 * a_1 as u128 * b_1 as u128 * n as u128)
                as i64)
            / (2 * a_1 * b_1))
            - c_2
    };
    let U_floor = |v: i64| {
        (((a_1 * b_2 + a_2 * b_1) * (v + c_2)
            - ceil_sqrt(((v + c_2) as u128).pow(2) - 4 * a_2 as u128 * b_2 as u128 * n as u128)
                as i64)
            / (2 * a_2 * b_2))
            - c_1
    };

    let S_w = || (1..w).map(V_floor).sum::<i64>();
    let S_h = || (1..h).map(U_floor).sum::<i64>();

    if w <= C2 {
        return s + S_w();
    }
    if h <= C2 {
        return s + S_h();
    }
    let u_4 = U_tan();
    let v_4 = V_floor(u_4);
    let u_5 = u_4 + 1;
    let v_5 = V_floor(u_5);
    let v_6 = u_4 + v_4;
    let u_7 = u_5 + v_5;
    s + triangle(v_6 - 1) - triangle(v_6 - u_5)
        + triangle(u_7 - u_5)
        + S_R(n, u_4, h - v_6, a_1, b_1, c_1, a_3, b_3, c_1 + c_2 + v_6)
        + S_R(n, w - u_7, v_5, a_3, b_3, c_1 + c_2 + u_7, a_2, b_2, c_2)
}

fn d(n: i64) -> i64 {
    let x_max = n.isqrt();
    2 * s(n, x_max) - x_max * x_max
}
#[must_use]
pub const fn icbrt(x: i64) -> i64 {
    unsafe { core::hint::assert_unchecked(x >= 1) };
    let mut rt = 1 << (1 + x.ilog2().div_ceil(3));
    let mut x_div_rt2 = (x / rt) / rt;
    while rt > x_div_rt2 {
        rt = ((rt << 1) + x_div_rt2) / 3;
        x_div_rt2 = (x / rt) / rt;
    }
    rt
}
/// Based on identity for `T_3(n)` in <https://arxiv.org/pdf/1206.3369>
/// O(n^2/3) time, O(1) space
#[must_use]
pub fn d3(n: i64) -> i64 {
    let cbrtn = icbrt(n);
    let mut ret = 0;
    let mut z = 1;
    while z <= cbrtn {
        let nz = n / z;
        let sqrtnz = nz.isqrt();
        ret += 2 * S_Q(nz, z + 1, sqrtnz) - sqrtnz * sqrtnz + nz / z;
        //ret += 2 * (s(nz, sqrtnz) - s(nz, z)) - sqrtnz * sqrtnz + nz / z;
        z += 1;
    }
    ret *= 3;
    ret += cbrtn * cbrtn * cbrtn;
    ret
}
