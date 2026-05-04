const MOD: u64 = 1e9 as u64 + 7;
pub fn main() {
    let start = std::time::Instant::now();
    // same idea as cses counting bits, but much uglier
    let mut pow2 = [1; 129]; // could instead maintain 2^k, 2^(k+1) and 2^(2k)
    for i in 0..128 {
        pow2[i + 1] = pow2[i] << 1;
        if pow2[i + 1] >= MOD {
            pow2[i + 1] -= MOD;
        }
    }
    let update = |k: usize, s2k1: [u64; 9], s: [u64; 9]| -> [u64; 9] {
        let mut res = s;
        for i in 0..9 {
            res[i] += s2k1[i];
            if res[i] >= MOD {
                res[i] -= MOD;
            }
        }
        res[1] = (res[1] + pow2[k] * s[0]) % MOD;
        res[2] = (res[2] + pow2[2 * k] * s[0] + pow2[k + 1] * s[1]) % MOD;
        res[3] += s[0];
        if res[3] >= MOD {
            res[3] -= MOD;
        }
        res[4] = (res[4] + s[1] + pow2[k] * (s[0] + s[3])) % MOD;
        res[5] = (res[5] + s[2] + pow2[k + 1] * (s[1] + s[4]) + pow2[2 * k] * (s[0] + s[3])) % MOD;
        res[6] = (res[6] + s[0] + 2 * s[3]) % MOD;
        res[7] = (res[7] + s[1] + 2 * s[4] + pow2[k] * (s[0] + 2 * s[3] + s[6])) % MOD;
        res[8] = (res[8]
            + s[2]
            + 2 * s[5]
            + pow2[2 * k] * (s[0] + 2 * s[3] + s[6])
            + (pow2[k + 1] * (s[1] + 2 * s[4] + s[7])) % MOD)
            % MOD;
        res
    };
    let mut s = [0; 9];
    s[0] = 1;
    let mut s2k1 = s;
    let mut n = 1e16 as usize;
    let mut k = 0;
    while n != 0 {
        if n & 1 == 1 {
            s = update(k, s2k1, s);
        }
        s2k1 = update(k, s2k1, s2k1);
        n >>= 1;
        k += 1;
    }
    println!("res = {}, took {:?}", s[8], start.elapsed());
}
