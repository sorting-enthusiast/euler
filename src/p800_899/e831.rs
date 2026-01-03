const M: usize = 142_857;
fn polymulmod(a: &[u128; 6], b: &[u128; 6], c: &mut [u128; 6]) {
    c.fill(0);
    for i in 0..6 {
        for j in 0..6 - i {
            c[i + j] += a[i] * b[j];
        }
    }
}

pub fn main() {
    let start = std::time::Instant::now();
    let mut exp = M;
    let mut x_store = [1, 3, 5, 5, 3, 1];
    let mut r_store = [1, 5, 10, 10, 5, 1];
    let mut tmp_store = [0; 6];
    let (mut x, mut r, mut tmp) = (&mut x_store, &mut r_store, &mut tmp_store);
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
    let mut res = tmp[5];
    let mut max_pow = 1;
    while max_pow <= res {
        max_pow *= 7;
    }
    max_pow /= 7;
    print!("res = ");
    for _ in 0..10 {
        print!("{}", res / max_pow);
        res %= max_pow;
        max_pow /= 7;
    }
    println!(", took {:?}", start.elapsed());
}
