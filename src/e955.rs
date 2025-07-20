const fn find_next(n: u128, nn: u128) -> (u128, u128, u128) {
    let mut k = n + 1;
    let mut kk = nn + k;

    let mut m = 1;
    let mut mm = 1;

    while kk != mm + nn {
        if kk > mm + nn {
            m += 1;
            mm += m;
        } else {
            k += 1;
            kk += k;
        }
    }
    (kk, k, m)
}
pub fn main() {
    let start = std::time::Instant::now();
    let mut n = 2;
    let mut t_n = 3;
    let mut index = 0;
    for i in 0..69 {
        let (new_t_n, new_n, incr) = find_next(n, t_n);
        println!(
            "{i}: {new_t_n}, {new_n}, {incr}, took {:?}",
            start.elapsed()
        );

        index += incr;
        n = new_n;
        t_n = new_t_n;
    }
    let end = start.elapsed();
    println!("res = {index}, took {end:?}");
}
