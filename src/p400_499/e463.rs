use std::collections::HashMap;

// plugging a bunch of values into oeis gives that f is just the bit-reversal of n: A030101
const fn f(n: u64) -> u64 {
    n.reverse_bits() >> n.leading_zeros()
}
fn test(n: u64) -> u64 {
    assert!((n + 1).is_power_of_two());
    (4u64.pow(n.trailing_ones()) - 1) / 3
}

// \sum_{i\le 2^k - 1} f(i) = \frac{4^k - 1}{3}
pub fn main() {
    assert_eq!(22, (1..=8).map(f).sum::<u64>());
    assert_eq!(3604, (1..=100).map(f).sum::<u64>());
    let k = 10;
    assert_eq!(
        (1..=(1 << k) - 1).step_by(2).map(f).sum::<u64>(),
        1u64 << ((k - 1) * 2)
    );
    dbg!(dbg!(test(63)) + (64..=100).map(f).sum::<u64>());

    println!("{:064b}", 3u64.pow(37));
    let mut cache = HashMap::new();
    dbg!(S(3u64.pow(37), &mut cache));
}
const MOD: u64 = 1e9 as _;
fn S(n: u64, cache: &mut HashMap<u64, u64>) -> u64 {
    if n < 3 {
        return n;
    }
    //if (n + 1) & n == 0 {
    //println!("{n}");
    //return (4u64.pow(n.trailing_ones()) - 1) / 3;
    //}
    if let Some(&v) = cache.get(&n) {
        return v;
    }
    let mut idk = [
        (n >> 1, 1),
        (((n - 1) >> 1) | 1, 2),
        ((((n - 1) >> 1) - 1) | 1, 3),
    ];
    idk.sort_by_key(|&(k, _)| k);
    let mut sum = (6 * S(idk[0].0, cache)) % MOD;
    for &(v, m) in &idk[1..] {
        for i in idk[0].0 + 1..=v {
            sum += m * (f(i) % MOD);
            sum %= MOD;
        }
    }
    sum += MOD - (8 * S((n - 3) >> 2, cache)) % MOD;
    for i in ((n - 3) >> 2) + 1..=(n - 1) >> 2 {
        sum += MOD - (3 * f(i)) % MOD;
    }
    sum += MOD - 1; // 4k + 1 branch definition does not align with actual value, leading to overcount
    sum %= MOD;

    cache.insert(n, sum);
    sum
}
