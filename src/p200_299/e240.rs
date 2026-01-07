const MOD: usize = 15 + 1;
const N: usize = 3;
const SIDES: usize = 6;
fn polymulmod(a: &[u128; MOD], b: &[u128; MOD], c: &mut [u128; MOD]) {
    c.fill(0);
    for i in 0..MOD {
        for j in 0..MOD - i {
            c[i + j] += a[i] * b[j];
        }
    }
}

pub fn main() {
    let start = std::time::Instant::now();
    let mut exp = N;
    let mut x_store = [0; MOD];
    let mut r_store = [0; MOD];
    let mut tmp_store = [0; MOD];
    let (mut x, mut r, mut tmp) = (&mut x_store, &mut r_store, &mut tmp_store);
    for i in 1..=SIDES {
        x[i] = 1;
    }
    r[0] = 1;
    while exp > 1 {
        if exp & 1 == 1 {
            polymulmod(r, x, tmp);
            std::mem::swap(&mut tmp, &mut r);
            //r = (r * x) % MODULO;
        }
        polymulmod(x, x, tmp);
        std::mem::swap(&mut tmp, &mut x);
        //x = (x * x) % MODULO;
        exp >>= 1;
    }
    polymulmod(r, x, tmp);
    dbg!(tmp[MOD - 1]);
}
