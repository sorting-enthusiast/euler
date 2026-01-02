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
        let d = epsilon / x;
        delta += d;
        epsilon -= x * d;
        gamma += 2 * delta - x * d;
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
#[must_use]
pub const fn icbrt(x: i64) -> i64 {
    let mut rt = 0;
    let mut rt_squared = 0;
    let mut rt_cubed = 0;
    while rt_cubed <= x {
        rt += 1;
        rt_squared += 2 * rt - 1;
        rt_cubed += 3 * rt_squared - 3 * rt + 1;
    }
    rt - 1
}
/// Based on identity for `T_3(n)` in <https://arxiv.org/pdf/1206.3369>
#[must_use]
pub const fn d3(n: i64) -> i64 {
    let cbrtn = icbrt(n);
    let mut ret = 0;
    let mut z = 1;
    while z <= cbrtn {
        let nz = n / z;
        let sqrtnz = nz.isqrt();
        ret += 2 * S_Q(nz, z + 1, sqrtnz) - sqrtnz * sqrtnz + nz / z;
        z += 1;
    }
    ret *= 3;
    ret += cbrtn * cbrtn * cbrtn;
    ret
}
