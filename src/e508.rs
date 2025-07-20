use std::time::Instant;
const fn fill_lookup() -> [(i64, i64); 256] {
    const BASE: [(i64, i64); 8] = [
        (1, 0),
        (-1, 1),
        (0, -2),
        (2, 2),
        (-4, 0),
        (4, -4),
        (0, 8),
        (-8, -8),
    ];
    let mut ret = [(0, 0); 256];
    let mut i: usize = 0;
    while i != 256 {
        let (mut a, mut b) = (0, 0);
        let mut j = i;
        while j != 0 {
            let (a1, b1) = BASE[j.trailing_zeros() as usize];
            (a, b) = (a + a1, b + b1);
            j &= j - 1;
        }
        ret[i] = (a, b);
        i += 1;
    }
    ret
}
const fn from(mut rep: u128) -> (i64, i64) {
    const LOOKUP: [(i64, i64); 256] = fill_lookup();

    let (mut a, mut b) = (0, 0);
    let mut shift = 0;
    while rep != 0 {
        let (a2, b2) = LOOKUP[(rep & 0xFF) as usize];
        a += a2 << shift;
        b += b2 << shift;
        rep >>= 8;
        shift += 4;
    }
    (a, b)
}
const fn b_to_nth(n: u32) -> (i64, i64) {
    const BASE_NORMALIZED: [(i64, i64); 8] = [
        (1, 0),
        (-1, 1),
        (0, -1),
        (1, 1),
        (-1, 0),
        (1, -1),
        (0, 1),
        (-1, -1),
    ];
    let shift = n >> 1;
    let (a, b) = BASE_NORMALIZED[n as usize & 7];
    (a << shift, b << shift)
}
const fn to((mut a, mut b): (i64, i64)) -> u128 {
    let mut rep = 0;
    while a != 0 || b != 0 {
        let next =
            ((a | b).trailing_zeros() << 1) | (a.trailing_zeros() == b.trailing_zeros()) as u32;
        rep |= 1 << next;
        let (c, d) = b_to_nth(next);
        a -= c;
        b -= d;
    }
    rep
}

const fn popcnt((mut a, mut b): (i64, i64)) -> u128 {
    let mut cnt = 0;
    while a != 0 || b != 0 {
        let next =
            ((a | b).trailing_zeros() << 1) | (a.trailing_zeros() == b.trailing_zeros()) as u32;
        let (c, d) = b_to_nth(next);
        a -= c;
        b -= d;
        cnt += 1;
    }
    cnt
}

const fn is_divisible_by_b_to_nth((a, b): (i64, i64), n: usize) -> bool {
    const LUT: [(i64, i64); 8] = [
        (1, 0),
        (-1, -1),
        (0, 1),
        (1, -1),
        (-1, 0),
        (1, 1),
        (0, -1),
        (-1, 1),
    ];
    let (a, b) = mul_gaussian(LUT[n & 7], (a, b));
    (a | b).trailing_zeros() as usize >= (n >> 1) + (n & 1)
    /* let new = [(a | b), (a - b) | (a + b)];
    new[n & 1].trailing_zeros() as usize >= (n >> 1) + (n & 1) */
}
pub const fn add_gaussian((a, b): (i64, i64), (c, d): (i64, i64)) -> (i64, i64) {
    (a + c, b + d)
}
pub const fn mul_gaussian((a, b): (i64, i64), (c, d): (i64, i64)) -> (i64, i64) {
    (a * c - b * d, a * d + b * c)
}
/*
// tetris lmfao
fn B_recursive(l: i64, depth: u8) -> u128 {
    if l < 16 {
        return [
            0, 20, 75, 175, 310, 515, 750, 1065, 1380, 1795, 2265, 2815, 3355, 4045, 4755, 5560,
        ][l as usize];
    }
    println!("{depth}: l = {l}, in");
    let ret = 4 * B_recursive(l >> 1, depth + 1)
        + if l & 1 == 1 {
            (4 * l * l) as u128
                + (0..l)
                    .flat_map(|i| [(l, l - 2 * i), (-l, -l + 1 + 2 * i)])
                    .map(popcnt)
                    .sum::<u128>()
                + 2 * (0..l)
                    .map(|i| (1 + (l / 2), -l / 2 + i))
                    .map(popcnt)
                    .sum::<u128>()
                + 3 * l as u128
                + popcnt((l, -l))
        } else {
            (4 * (l + 1) * (l + 1)) as u128
                - (0..=l)
                    .flat_map(|i| [(l + 1, l - 2 * i), (-(l + 1), -l + 1 + 2 * i)])
                    .map(popcnt)
                    .sum::<u128>()
                - 2 * (0..=l)
                    .map(|i| (-l / 2, -l / 2 + i))
                    .map(popcnt)
                    .sum::<u128>() // can be optimized: half of the points here are added back next iteration
                - 3 * (l + 1) as u128
                + popcnt((-l - 1, l + 1))
        };
    println!("{depth}: l = {l}, out");
    ret
}
 */
// linear in l, more specifically ~6l, seems to be easy to optimize to log_2(l)
fn B(mut L: i64) -> u128 {
    const PRECOMPUTED: [i128; 16] = [
        0, 20, 75, 175, 310, 515, 750, 1065, 1380, 1795, 2265, 2815, 3355, 4045, 4755, 5560,
    ];

    let expected_popcnt_calls = 6 * L;
    let mut popcnt_calls = 0;
    let mut ret = 0;
    let mut shift = 0;
    while L >= 16 {
        println!("l = {L}");
        ret += (if L & 1 == 1 {
            (4 * L * L) as u128
                + (0..L)
                    .flat_map(|i| [(L, L - 2 * i), (-L, -L + 1 + 2 * i)])
                    .map(popcnt)
                    .sum::<u128>()
                + 2 * (0..L)
                    .map(|i| (1 + (L / 2), -L / 2 + i))
                    .map(popcnt)
                    .sum::<u128>()
                + 3 * L as u128
                + popcnt((L, -L))
        } else {
            (4 * (L + 1) * (L + 1)) as u128
                    - (0..=L)
                        .flat_map(|i| [(L + 1, L - 2 * i), (-(L + 1), -L + 1 + 2 * i)])
                        .map(popcnt)
                        .sum::<u128>()
                    - 2 * (0..=L)
                        .map(|i| (-L / 2, -L / 2 + i))
                        .map(popcnt)
                        .sum::<u128>() // can be optimized: half of the points here are added back next iteration
                    - 3 * (L + 1) as u128
                + popcnt((-L - 1, L + 1))
        } as i128)
            << shift;
        popcnt_calls += 1 + 3 * (L + 1 - (L & 1));
        shift += 2;
        L >>= 1;
    }
    println!("{popcnt_calls} ~ {expected_popcnt_calls}");
    ret += PRECOMPUTED[L as usize] << shift;
    ret as u128
}

pub fn main() {
    const L: i64 = 1e15 as i64;
    const MOD: u128 = 1e9 as u128 + 7;

    let start = Instant::now();
    let res = B_rec_v3(L).0;
    let end = start.elapsed();
    println!("res = {:?}, took {end:?}", res % MOD);

    let start = Instant::now();
    let res = B_rec_v4(L);
    let end = start.elapsed();
    println!("res = {:?}, took {end:?}", res % MOD);
}
/*
fn B_rec_v2(L: i64, depth: u8) -> (u128, [[u128; 2]; 2]) {
    if L <= 1 {
        return [(0, [[0; 2]; 2]), (20, [[7, 8], [2, 3]])][L as usize];
    }
    let (square, [[l, r], [t, b]]) = B_rec_v2(L >> 1, depth + 1);
    println!("{depth}: {L}");

    let top_right = popcnt((L, L));
    let bottom_left = popcnt((-L, -L));
    let top_left = popcnt((-L, L));
    let bottom_right = popcnt((L, -L));
    if L & 1 == 0 {
        let smaller_square = square - l - r - t - b;
        /* let mut current_square =
        (smaller_square + ((L - 1) as u128 * (L - 1) as u128) as u64 % MOD) << 2; */
        let mut current_square = (L as u128 - 1)
            .checked_mul(L as u128 - 1)
            .and_then(|x| smaller_square.checked_add(x))
            .and_then(|x| x.checked_shl(2))
            .unwrap();

        let bottom2rows =
            4 * (r + 1 + (L as u128)) - popcnt((L + 1, -L)) - popcnt((-L - 1, -L + 1));
        current_square += bottom2rows;

        let bottom_row_without_edges =
            2 * r + (L as u128) + 1 - bottom_left - bottom_right - popcnt((L + 1, -L));

        let top_row = 2 * l + (L as u128) + 1 - popcnt((L + 1, L));
        current_square += top_row;

        let top_row_without_edges = top_row - top_right - top_left;

        let left_triples = 3 * (b + L as u128 - 1);
        let right_triples = 3 * (t + L as u128 - 1);
        current_square += left_triples + right_triples;

        let left = 2 * (b + L as u128 - 1) + top_left + popcnt((-L, -L + 1)) + bottom_left;
        let right = 2 * (t + L as u128 - 1) + top_right + popcnt((L, -L + 1)) + bottom_right;

        (
            current_square,
            [
                [left, right],
                [top_row_without_edges, bottom_row_without_edges],
            ],
        )
    } else {
        //something here overflows
        let mut current_square = //4 * (square + (L * L) as u64 % MOD);
        (L as u128)
            .checked_mul(L as u128)
            .and_then(|x| square.checked_add(x))
            .and_then(|x| x.checked_shl(2))
            .unwrap();

        let top_row_without_edges = 2 * l + (3 * L as u128) - top_left;

        let bottom_row = 2
            * (0..L)
                .into_par_iter()
                .map(|i| (1 + (L / 2), -L / 2 + i))
                .map(popcnt)
                .reduce(|| 0, |acc, p| acc.checked_add(p).unwrap())
            + (3 * L as u128)
            + bottom_right;
        let bottom_row_without_edges = bottom_row - bottom_right - bottom_left;

        current_square += bottom_row;

        let left_holes = (0..L)
            .into_par_iter()
            .map(|i| (-L, -L + 1 + 2 * i))
            .map(popcnt)
            .reduce(|| 0, |acc, p| acc.checked_add(p).unwrap());
        let left = left_holes + b + (L as u128) - 2 + popcnt((-L, -L + 2)) + top_left + bottom_left;

        let right_holes = (0..L)
            .into_par_iter()
            .map(|i| (L, L - 2 * i))
            .map(popcnt)
            .reduce(|| 0, |acc, p| acc.checked_add(p).unwrap());
        let right = right_holes + t + (L as u128) - 2
            + bottom_right
            + popcnt((L, -L + 1))
            + popcnt((L, L - 1));
        current_square += left_holes + right_holes;

        (
            current_square,
            [
                [left, right],
                [top_row_without_edges, bottom_row_without_edges],
            ],
        )
    }
}
 */
//calculating for odd sizes is trivial, but maintaining extra data becomes expensive. nvm, is now the fastest way
const fn B_rec_v3(L: i64) -> (u128, [[u128; 2]; 2]) {
    if L <= 1 {
        return [(0, [[0; 2]; 2]), (20, [[7, 8], [2, 3]])][L as usize];
    }

    if L & 1 == 0 {
        let (square, [[l, r], [t, b]]) = B_rec_v3(L >> 1);

        let top_right = popcnt((L, L));
        let bottom_left = popcnt((-L, -L));
        let top_left = popcnt((-L, L));
        let bottom_right = popcnt((L, -L));
        let smaller_square = square - l - r - t - b;

        let mut current_square = ((L as u128 - 1) * (L as u128 - 1) + smaller_square) << 2;

        let bottom2rows =
            4 * (r + 1 + (L as u128)) - popcnt((L + 1, -L)) - popcnt((-L - 1, -L + 1));
        current_square += bottom2rows;

        let bottom_row_without_edges =
            2 * r + (L as u128) + 1 - bottom_left - bottom_right - popcnt((L + 1, -L));

        let top_row = 2 * l + (L as u128) + 1 - popcnt((L + 1, L));
        current_square += top_row;

        let top_row_without_edges = top_row - top_right - top_left;

        let left_triples = 3 * (b + L as u128 - 1);
        let right_triples = 3 * (t + L as u128 - 1);
        current_square += left_triples + right_triples;

        let left = 2 * (b + L as u128 - 1) + top_left + popcnt((-L, -L + 1)) + bottom_left;
        let right = 2 * (t + L as u128 - 1) + top_right + popcnt((L, -L + 1)) + bottom_right;

        (
            current_square,
            [
                [left, right],
                [top_row_without_edges, bottom_row_without_edges],
            ],
        )
    } else {
        let ceil_div_2 = 1 + (L >> 1);
        let (square, [[l, r], [t, b]]) = B_rec_v3(ceil_div_2);
        //println!("{depth}: {L}");

        let top_right = popcnt((L, L));
        //let bottom_left = popcnt((-L, -L));
        let top_left = popcnt((-L, L));
        //let bottom_right = popcnt((L, -L));

        let mut current_square = ((L as u128 + 2) * (L as u128 + 2) + square) << 2;
        current_square -= 4 * (l + L as u128 + 2);
        current_square -= 2 * r
            + L as u128
            + 2
            + 2 * popcnt((ceil_div_2, -ceil_div_2))
            + 3
            + popcnt((ceil_div_2, ceil_div_2))
            + 2;
        current_square -= 3 * b + 3 * L as u128;
        current_square -= 3 * t + 3 * L as u128;

        let left = left_column(L);
        /* b
        + L as u128
        + bottom_left
        + (0..L)
            .into_par_iter()
            .map(|i| (-L, L - 2 * i))
            .map(popcnt)
            .reduce(|| 0, |acc, p| acc.checked_add(p).unwrap()); */
        let right = right_column(L);
        /* t
        + L as u128
        + bottom_right
        + (0..L)
            .into_par_iter()
            .map(|i| (L, L - 1 - 2 * i))
            .map(popcnt)
            .reduce(|| 0, |acc, p| acc.checked_add(p).unwrap()); */
        let bottom_row_without_edges = 2 * r + 3 * (L as u128 + 2)
            - 2 * popcnt((ceil_div_2, -ceil_div_2))
            - 3
            - 2 * popcnt((ceil_div_2, ceil_div_2))
            - 3
            - popcnt((-L, -L));
        let top_row_without_edges = top_row(L)/* 2
            * (0..L)
                .into_par_iter()
                .map(|i| (-L / 2, i - L / 2))
                .map(popcnt)
                .reduce(|| 0, |acc, p| acc.checked_add(p).unwrap())
            + 3 * L as u128 */
            - top_left - top_right;

        (
            current_square,
            [
                [left, right],
                [top_row_without_edges, bottom_row_without_edges],
            ],
        )
    }
}

const fn top_row(L: i64) -> u128 {
    if L == 0 {
        return 0;
    }
    2 * left_column(L >> 1)
        + if L & 1 == 0 {
            L as u128 - popcnt((-L / 2, L / 2))
        } else {
            3 * L as u128 + 1 + popcnt((-L / 2, 1 + L / 2))
        }
}

const fn bottom_row(L: i64) -> u128 {
    if L == 0 {
        0
    } else if L & 1 == 0 {
        2 * right_column(L >> 1) + L as u128 - popcnt((L / 2, L / 2))
    } else {
        2 * right_column(1 + (L >> 1)) + 3 * L as u128 + 1
            - 2 * popcnt((1 + L / 2, -1 - L / 2))
            - popcnt((1 + L / 2, 1 + L / 2))
    }
}
const fn left_column(L: i64) -> u128 {
    if L == 0 {
        return 0;
    }
    2 * L as u128
        + if L & 1 == 0 {
            2 * bottom_row(L >> 1) - popcnt((-L / 2, -L / 2))
        } else {
            bottom_row(L >> 1) + bottom_row(1 + (L >> 1)) + 1 + popcnt((1 + L / 2, -L / 2))
                - popcnt((1 + L / 2, -L / 2 - 1))
                - popcnt((-L / 2 - 1, -L / 2 - 1))
        }
}

const fn right_column(L: i64) -> u128 {
    if L == 0 {
        return 0;
    }
    2 * L as u128
        + if L & 1 == 0 {
            2 * top_row(L >> 1) - popcnt((-L / 2, L / 2))
        } else {
            top_row(L >> 1) + top_row(1 + (L >> 1)) + 1 - popcnt((-L / 2 - 1, L / 2 + 1))
        }
}
// can actually run this at compile time lmfao
const fn B_rec_v4(L: i64) -> u128 {
    if L == 0 {
        return 0;
    }
    let ret = 4 * B_rec_v4(L >> 1) + 4 * (L as u128) * (L as u128) + 3 * L as u128;
    if L & 1 == 0 {
        ret - top_row(L >> 1) - 2 * left_column(L >> 1) - bottom_row(L >> 1)
            + popcnt((-L / 2, -L / 2))
    } else {
        ret + 2 * L as u128
            + 1
            + bottom_row(1 + L / 2)
            + 2 * right_column(1 + L / 2)
            + top_row(1 + L / 2)
            - popcnt((-1 - L / 2, -1 - L / 2))
            - 3 * popcnt((1 + L / 2, -1 - L / 2))
            - popcnt((-1 - L / 2, 1 + L / 2))
            - 2 * popcnt((1 + L / 2, 1 + L / 2))
    }
}
