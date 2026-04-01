// can be done much better using euler transform in flint, using flint from rust on windows is a bitch though
// polynomial division using flint would also be great
pub fn main() {
    const N: usize = 1e7 as _;
    const MOD: i64 = 1e9 as i64 + 7;
    let start = std::time::Instant::now();
    let mut s = vec![0; N + 1];
    let mut n = 0;
    let mut tn = 0;
    while tn <= N {
        s[tn] = 1;
        n += 1;
        tn += n;
    }
    for i in 1..=N {
        s[i] += s[i - 1];
    }
    for i in 0..=N {
        let mut val = 0;
        for j in (2..=i.isqrt()).step_by(2) {
            if j & 3 == 0 {
                val -= s[i - j * j];
            } else {
                val += s[i - j * j];
            }
        }
        s[i] = (s[i] + (val << 1)) % MOD;
        if s[i] < 0 {
            s[i] += MOD;
        }
    }
    let res = s[N] - s[0];
    println!("res = {res}, took {:?}", start.elapsed());
}
