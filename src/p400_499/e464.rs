use itertools::Itertools;

use crate::utils::{fenwick::FenwickTreeI64, multiplicative_function_summation::mobius_sieve};
const N: usize = 2e7 as _;

#[derive(Clone, Copy, Default)]
struct Pt {
    x: i32,
    yr: i32,
    //idx: usize,
}

// pts are initially sorted by index, and returned sorted by x
fn cdq(pts: &mut [Pt], tmp: &mut [Pt], bit: &mut FenwickTreeI64) -> i64 {
    let len = pts.len();
    if len < 2 {
        return 0;
    }
    let m = len >> 1;
    let mut res = cdq(&mut pts[..m], &mut tmp[..m], bit) + cdq(&mut pts[m..], &mut tmp[m..], bit);

    let mut i = 0;
    let mut j = m;
    let mut k = 0;

    loop {
        if pts[i].x <= pts[j].x {
            bit.add(pts[i].yr as _, 1);
            tmp[k] = pts[i];
            i += 1;
            if i >= m {
                tmp[k + 1..].copy_from_slice(&pts[j..]);
                for pt in &pts[j..] {
                    res += bit.sum(pt.yr as _);
                }
                break;
            }
        } else {
            // No more left points will have x <= pts[j].x,
            // or in other words, contribute to the points dominated by pts[j]
            // Therefore, we can safely forget about pts[j] from here on out
            res += bit.sum(pts[j].yr as _);
            tmp[k] = pts[j];
            j += 1;
            if j == len {
                tmp[k + 1..].copy_from_slice(&pts[i..m]);
                /* for pt in &pts[i..m] {
                    bit.add(pt.yr, 1);
                } */
                break;
            }
        }
        k += 1;
    }

    // Roll back BIT updates: one addition for each point from the left side (or until the right points ran out)
    for pt in &pts[..i] {
        bit.add(pt.yr as _, -1);
    }
    //bit.0.fill(0); // will accomplish the same, but is slower asymptotically

    // Copy back so parent sees pts sorted by x
    // Alternatively, could take inspiration from ping-pong merging
    // to get rid of unproductive memcpy calls
    pts.copy_from_slice(tmp);
    res
}

// 0.3 seconds faster lmao
fn pingpong_cdq(pts: &mut [Pt], tmp: &mut [Pt], bit: &mut FenwickTreeI64) -> i64 {
    let len = pts.len();
    if len < 4 {
        return cdq(pts, tmp, bit);
    }
    let m = len >> 1;
    let q1 = m >> 1;
    let q2 = (len - m) >> 1;
    let mut res = pingpong_cdq(&mut pts[..q1], &mut tmp[..q1], bit)
        + pingpong_cdq(&mut pts[q1..m], &mut tmp[q1..m], bit)
        + pingpong_cdq(&mut pts[m..][..q2], &mut tmp[m..][..q2], bit)
        + pingpong_cdq(&mut pts[m..][q2..], &mut tmp[m..][q2..], bit);
    {
        let mut i = 0;
        let mut j = q1;
        let mut k = 0;
        loop {
            if pts[..m][i].x <= pts[..m][j].x {
                bit.add(pts[..m][i].yr as _, 1);
                tmp[..m][k] = pts[..m][i];
                i += 1;
                if i >= q1 {
                    tmp[..m][k + 1..].copy_from_slice(&pts[..m][j..]);
                    for pt in &pts[..m][j..] {
                        res += bit.sum(pt.yr as _);
                    }
                    break;
                }
            } else {
                // No more left points will have x <= pts[j].x,
                // or in other words, contribute to the points dominated by pts[j]
                // Therefore, we can safely forget about pts[j] from here on out
                res += bit.sum(pts[..m][j].yr as _);
                tmp[..m][k] = pts[..m][j];
                j += 1;
                if j == m {
                    tmp[..m][k + 1..].copy_from_slice(&pts[i..q1]);
                    break;
                }
            }
            k += 1;
        }
        for pt in &pts[..m][..i] {
            bit.add(pt.yr as _, -1);
        }
    }
    {
        let mut i = 0;
        let mut j = q2;
        let mut k = 0;

        loop {
            if pts[m..][i].x <= pts[m..][j].x {
                bit.add(pts[m..][i].yr as _, 1);
                tmp[m..][k] = pts[m..][i];
                i += 1;
                if i >= q2 {
                    tmp[m..][k + 1..].copy_from_slice(&pts[m..][j..]);
                    for pt in &pts[m..][j..] {
                        res += bit.sum(pt.yr as _);
                    }
                    break;
                }
            } else {
                // No more left points will have x <= pts[j].x,
                // or in other words, contribute to the points dominated by pts[j]
                // Therefore, we can safely forget about pts[j] from here on out
                res += bit.sum(pts[m..][j].yr as _);
                tmp[m..][k] = pts[m..][j];
                j += 1;
                if j == len - m {
                    tmp[m..][k + 1..].copy_from_slice(&pts[m..][i..q2]);
                    break;
                }
            }
            k += 1;
        }

        for pt in &pts[m..][..i] {
            bit.add(pt.yr as _, -1);
        }
    }
    {
        let mut i = 0;
        let mut j = m;
        let mut k = 0;

        loop {
            if tmp[i].x <= tmp[j].x {
                bit.add(tmp[i].yr as _, 1);
                pts[k] = tmp[i];
                i += 1;
                if i >= m {
                    pts[k + 1..].copy_from_slice(&tmp[j..]);
                    for pt in &tmp[j..] {
                        res += bit.sum(pt.yr as _);
                    }
                    break;
                }
            } else {
                // No more left points will have x <= pts[j].x,
                // or in other words, contribute to the points dominated by pts[j]
                // Therefore, we can safely forget about pts[j] from here on out
                res += bit.sum(tmp[j].yr as _);
                pts[k] = tmp[j];
                j += 1;
                if j == len {
                    pts[k + 1..].copy_from_slice(&tmp[i..m]);
                    break;
                }
            }
            k += 1;
        }

        for pt in &tmp[..i] {
            bit.add(pt.yr as _, -1);
        }
    }
    res
}

fn solve3() {
    {
        let start = std::time::Instant::now();

        let mob = mobius_sieve(N + 1);
        let mut xys: Vec<(i32, i32)> = vec![(0, 0); N + 1];
        xys[1] = (100, -99);
        for i in 2..=N {
            let (mut k1, mut k2) = xys[i - 1];
            if mob[i] == 1 {
                k1 += 100;
                k2 -= 99;
            } else if mob[i] == -1 {
                k1 -= 99;
                k2 += 100;
            }
            xys[i] = (k1, k2);
        }

        // Coordinate-compress y
        let mut ys = xys.iter().map(|&(_, y)| y).collect_vec();
        ys.sort_unstable();
        ys.dedup();
        let y_rank = |y| ys.partition_point(|&v| v < y) as i32;
        let y_max = ys.len();

        let mut pts = (0..=N)
            .map(|i| Pt {
                x: xys[i].0,
                yr: y_rank(xys[i].1),
                //idx: i,
            })
            .collect_vec();
        let mut tmp = vec![Pt::default(); N + 1]; // can halve the size of tmp easily, using similar tricks to my mergesorts
        //let mut ans = vec![0; N + 1];
        let mut bit = FenwickTreeI64::new(y_max, 0);
        //dbg!(start.elapsed());
        // 4) CDQ on index
        let count = cdq(&mut pts, &mut tmp, &mut bit);

        //let count = ans[2..=N].iter().sum::<i64>();

        dbg!(count, start.elapsed());
    }

    {
        let start = std::time::Instant::now();

        let mob = mobius_sieve(N + 1);
        let mut xys = vec![(0, 0); N + 1];

        xys[1] = (100, -99);
        for i in 2..=N {
            let (mut k1, mut k2) = xys[i - 1];
            if mob[i] == 1 {
                k1 += 100;
                k2 -= 99;
            } else if mob[i] == -1 {
                k1 -= 99;
                k2 += 100;
            }
            xys[i] = (k1, k2);
        }

        // Coordinate-compress y
        let mut ys = xys.iter().map(|&(_, y)| y).collect_vec();
        ys.sort_unstable();
        ys.dedup();
        let y_rank = |y| ys.partition_point(|&v| v < y) as i32;
        let y_max = ys.len();

        let mut pts = (0..=N)
            .map(|i| Pt {
                x: xys[i].0,
                yr: y_rank(xys[i].1),
                //idx: i,
            })
            .collect_vec();
        let mut tmp = vec![Pt::default(); N + 1];
        //let mut ans = vec![0; N + 1];
        let mut bit = FenwickTreeI64::new(y_max, 0);
        // dbg!(start.elapsed());
        // 4) CDQ on index
        let count = pingpong_cdq(&mut pts, &mut tmp, &mut bit);

        //let count = ans[2..=N].iter().sum::<i64>();

        dbg!(count, start.elapsed());
    }
}

pub fn main() {
    solve2();
    solve3();
    //solve();
}
fn solve() {
    let start = std::time::Instant::now();

    let mob = mobius_sieve(N + 1);
    let mut xys = vec![(0, 0); N + 1];

    xys[1] = (100, -99);
    let mut count = 0;
    let mut prev = 0;
    for i in 2..=N {
        /* if i & 65535 == 0 {
            dbg!((i, count, start.elapsed()));
        } */
        let (mut k1, mut k2) = xys[i - 1];
        if mob[i] == 1 {
            k1 += 100;
            k2 -= 99;
        } else if mob[i] == -1 {
            k1 -= 99;
            k2 += 100;
        } else {
            prev += 1;
            count += prev;
            xys[i] = xys[i - 1];

            continue;
        }
        prev = xys[..i]
            .iter()
            .filter(|&&(x, y)| x <= k1 && y <= k2)
            .count();
        count += prev;
        xys[i] = (k1, k2);
    }
    dbg!(count, start.elapsed());
}
fn solve2() {
    {
        let start = std::time::Instant::now();

        let mob = mobius_sieve(N + 1);
        let mut xs = vec![0; N + 1];
        let mut ys = vec![0; N + 1];

        (xs[1], ys[1]) = (100, -99);
        for i in 2..=N {
            let (mut k1, mut k2) = (xs[i - 1], ys[i - 1]);
            if mob[i] == 1 {
                k1 += 100;
                k2 -= 99;
            } else if mob[i] == -1 {
                k1 -= 99;
                k2 += 100;
            }
            (xs[i], ys[i]) = (k1, k2);
        }
        let mut count = (N * (N + 1)) >> 1;
        let mut tmp = vec![0; N + 1];
        count -= pingpong_inversions(&mut xs, &mut tmp);
        count -= pingpong_inversions(&mut ys, &mut tmp);
        /*
            for i in 1..=N {
                let (k1, k2) = (xs[i], ys[i]);

                count -= xs[..i].iter().filter(|&&x| x > k1).count() // inversions in x
                    +ys[..i].iter().filter(|&&y| y > k2).count() // inversions in y
                    - both, but no points will fit this criterion
        ;
            } */
        dbg!(count, start.elapsed());
    }
    {
        let start = std::time::Instant::now();

        let mob = mobius_sieve(N + 1);
        let mut xs = vec![0; N + 1];
        let mut ys = vec![0; N + 1];

        (xs[1], ys[1]) = (100, -99);
        for i in 2..=N {
            let (mut k1, mut k2) = (xs[i - 1], ys[i - 1]);
            if mob[i] == 1 {
                k1 += 100;
                k2 -= 99;
            } else if mob[i] == -1 {
                k1 -= 99;
                k2 += 100;
            }
            (xs[i], ys[i]) = (k1, k2);
        }
        let mut count = (N * (N + 1)) >> 1;
        let mut tmp = vec![0; N + 1];
        count -= inversions(&mut xs, &mut tmp);
        count -= inversions(&mut ys, &mut tmp);
        /*
            for i in 1..=N {
                let (k1, k2) = (xs[i], ys[i]);

                count -= xs[..i].iter().filter(|&&x| x > k1).count() // inversions in x
                    +ys[..i].iter().filter(|&&y| y > k2).count() // inversions in y
                    - both, but no points will fit this criterion
        ;
            } */
        dbg!(count, start.elapsed());
    }
}
fn inversions(pts: &mut [i32], tmp: &mut [i32]) -> usize {
    let len = pts.len();
    if len < 2 {
        return 0;
    }
    let m = len >> 1;
    let mut res =
        inversions(&mut pts[..m], &mut tmp[..m]) + inversions(&mut pts[m..], &mut tmp[m..]);

    let mut i = 0;
    let mut j = m;
    let mut k = 0;
    loop {
        if pts[i] <= pts[j] {
            tmp[k] = pts[i];
            i += 1;
            if i >= m {
                tmp[k + 1..].copy_from_slice(&pts[j..]);
                break;
            }
        } else {
            res += m - i; // current is smaller than all remaining left points
            tmp[k] = pts[j];
            j += 1;
            if j == len {
                tmp[k + 1..].copy_from_slice(&pts[i..m]);
                break;
            }
        }
        k += 1;
    }

    // Copy back so parent sees pts sorted by x
    // Alternatively, could take inspiration from ping-pong merging
    // to get rid of unproductive memcpy calls
    pts.copy_from_slice(tmp);
    res
}
fn pingpong_inversions(pts: &mut [i32], tmp: &mut [i32]) -> usize {
    let len = pts.len();
    if len < 4 {
        return inversions(pts, tmp);
    }
    let m = len >> 1;
    let q1 = m >> 1;
    let q2 = (len - m) >> 1;
    let mut res = pingpong_inversions(&mut pts[..q1], &mut tmp[..q1])
        + pingpong_inversions(&mut pts[q1..m], &mut tmp[q1..m])
        + pingpong_inversions(&mut pts[m..][..q2], &mut tmp[m..][..q2])
        + pingpong_inversions(&mut pts[m..][q2..], &mut tmp[m..][q2..]);
    {
        let mut i = 0;
        let mut j = q1;
        let mut k = 0;
        loop {
            if pts[..m][i] <= pts[..m][j] {
                tmp[..m][k] = pts[..m][i];
                i += 1;
                if i >= q1 {
                    tmp[..m][k + 1..].copy_from_slice(&pts[..m][j..]);
                    break;
                }
            } else {
                res += q1 - i;
                tmp[..m][k] = pts[..m][j];
                j += 1;
                if j == m {
                    tmp[..m][k + 1..].copy_from_slice(&pts[i..q1]);
                    break;
                }
            }
            k += 1;
        }
    }
    {
        let mut i = 0;
        let mut j = q2;
        let mut k = 0;

        loop {
            if pts[m..][i] <= pts[m..][j] {
                tmp[m..][k] = pts[m..][i];
                i += 1;
                if i >= q2 {
                    tmp[m..][k + 1..].copy_from_slice(&pts[m..][j..]);
                    break;
                }
            } else {
                res += q2 - i;
                tmp[m..][k] = pts[m..][j];
                j += 1;
                if j == len - m {
                    tmp[m..][k + 1..].copy_from_slice(&pts[m..][i..q2]);
                    break;
                }
            }
            k += 1;
        }
    }
    {
        let mut i = 0;
        let mut j = m;
        let mut k = 0;

        loop {
            if tmp[i] <= tmp[j] {
                pts[k] = tmp[i];
                i += 1;
                if i >= m {
                    pts[k + 1..].copy_from_slice(&tmp[j..]);
                    break;
                }
            } else {
                res += m - i;
                pts[k] = tmp[j];
                j += 1;
                if j == len {
                    pts[k + 1..].copy_from_slice(&tmp[i..m]);
                    break;
                }
            }
            k += 1;
        }
    }
    res
}
