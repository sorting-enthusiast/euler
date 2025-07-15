use itertools::Itertools;

use crate::utils::bit_array::BitArray;

const fn modexp(mut x: usize, mut exp: usize, modulo: usize) -> usize {
    if exp == 0 {
        return 1;
    }
    let mut r = 1;
    while exp > 1 {
        if exp & 1 == 1 {
            r = (r * x) % modulo;
        }
        x = (x * x) % modulo;
        exp >>= 1;
    }
    (r * x) % modulo
}
fn find_divisor(p: usize) -> (usize, usize) {
    assert_eq!(
        p & 3,
        1,
        "cannot find divisor, prime is not equivalent to 1 mod 4"
    );
    let r = (2..p)
        .map(|i| modexp(i, (p - 1) >> 2, p))
        .find(|r| p - 1 == (r * r) % p)
        .expect("divisor was not found, something is very wrong");
    let mut a = p;
    let mut b = r;

    // Euclidean algorithm
    while b * b > p {
        let tmp = a % b;
        a = b;
        b = tmp;
    }

    (b, (p - b * b).isqrt())
}
fn sieve(n: usize) -> Vec<usize> {
    let mut res = vec![0; n];
    let mut pow = vec![0; n];
    //let mut composite = BitArray::zeroed(n);
    let mut primes = vec![];
    let mut prime_power_divisors = vec![];
    res[1] = 1;
    for i in 2..n {
        if res[i] == 0
        /* !composite.get(i) */
        {
            //dbg!(!composite.get(i) && (res[i] == 0));
            primes.push(i);
            pow[i] = i;
            res[i] = i
                + 1
                + match i & 3 {
                    1 => {
                        let (x, y) = find_divisor(i);
                        assert_eq!(x * x + y * y, i);
                        if i * i < n {
                            prime_power_divisors.push(Some((
                                (x as i64, y as i64),
                                vec![
                                    (1, 0),
                                    (x as i64, y as i64),
                                    (x as i64, -(y as i64)),
                                    (y as i64, x as i64),
                                    (y as i64, -(x as i64)),
                                    (i as i64, 0),
                                ],
                            )));
                        }
                        (x + y) << usize::from(x != y)
                    }
                    2 => {
                        if i * i < n {
                            prime_power_divisors.push(None); /* 
                            prime_power_divisors
                            .push(Some(((1, 1), vec![(1, 0), (1, 1), (1, -1), (2, 0)]))); */
                        }
                        2
                    }
                    3 => {
                        if i * i < n {
                            prime_power_divisors.push(None);
                        }
                        0
                    }
                    _ => unsafe { core::hint::unreachable_unchecked() },
                }
        }
        for (j, &p) in primes.iter().enumerate() {
            if i * p >= n {
                break;
            }
            //composite.set(i * p);
            if i % p == 0 {
                pow[i * p] = pow[i] * p;
                let v = i / pow[i];
                if v != 1 {
                    assert_eq!(v * pow[i], i);
                    assert_ne!(v % p, 0);
                    res[i * p] = res[v] * res[p * pow[i]];
                } else
                //prime power
                if p & 3 != 3 {
                    //has complex divisors :(
                    if p == 2 {
                        res[i << 1] = res[i] + (i << 2);
                    } else {
                        let Some(((x, y), ref mut v)) = prime_power_divisors[j] else {
                            unsafe { core::hint::unreachable_unchecked() }
                        };
                        println!("{p}: {x} {y}, {}", i * p);
                        for t in v.iter() {
                            print!("{t:?} ");
                        }
                        println!();
                        let len = v.len();
                        for (a, b) in [
                            /* (1, 0),  */ (x, y),
                            (x, -y),
                            (y, x),
                            (y, -x),
                            (p as i64, 0),
                        ] {
                            for k in 0..len {
                                let (c, d) = v[k];
                                let acbd = a * c + b * d;
                                let bcad = b * c - a * d;
                                if acbd <= 0 || v.contains(&(acbd, bcad)) {
                                    continue;
                                }
                                v.push((acbd, bcad));
                            }
                        }
                        v[len..].sort_by_key(|(a, b)| a.pow(2) + b.pow(2));

                        for t in &v[len..] {
                            print!("{t:?} ");
                        }
                        println!("\n");
                        assert_eq!(v.iter().map(|t| t.1).sum::<i64>(), 0);
                        assert_eq!(res[i], v[..len].iter().map(|t| t.0 as usize).sum::<usize>());
                        res[i * p] = res[i] + v[len..].iter().map(|t| t.0 as usize).sum::<usize>();
                    }
                } else {
                    res[i * p] = res[i] + i * p;
                }
                break;
            }
            assert_eq!(res[i * p], 0);
            res[i * p] = res[i] * res[p];
            pow[i * p] = p;
        }
    }
    res
}

pub fn main() {
    const N: usize = 2 as usize;

    let start = std::time::Instant::now();
    let s = sieve(N + 1);
    //dbg!(s[64]);
    let count = s.iter().sum::<usize>();
    let end = start.elapsed();
    println!("res = {count}, took {end:?}");
    let measured = (2..=50)
        .map(sieve)
        .map(|s| s.iter().sum::<usize>())
        .collect_vec();
    let expected = [
        1, 6, 10, 23, 35, 55, 63, 92, 105, 161, 173, 225, 249, 289, 337, 398, 426, 491, 511, 655,
        687, 747, 771, 887, 968, 1080, 1120, 1224, 1268, 1492, 1524, 1649, 1697, 1833, 1929, 2098,
        2150, 2250, 2346, 2666, 2726, 2886, 2930, 3086, 3242, 3362, 3410, 3654, 3711,
    ];
    for t in measured {
        print!("{t} ");
    }
    println!();
    for t in expected {
        print!("{t} ");
    }
    println!();

    dbg!(sieve(10 + 1));
    dbg!(17_924_657_155_usize.abs_diff(count));
}
