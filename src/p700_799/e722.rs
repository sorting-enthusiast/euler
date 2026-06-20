// better ideas are based on the theory of modular forms, or on polylogarithms
pub fn main() {
    let start = std::time::Instant::now();
    //let x = (1. - 0.5f64.powi(25)).recip();
    let x_ln = 0.5f64.powi(25)
        + 0.5f64.powi(51)
        + 0.5f64.powi(75) / 3.
        + 0.5f64.powi(102)
        + 0.5f64.powi(125) / 5.;
    let mut sum = 0.;
    let mut c = 0.;
    for n in 1..1usize << 32 {
        let n = n as f64;
        let e = n.powi(15) / (n * x_ln).exp_m1();
        let y = e - c;
        let t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    println!("res = {sum:.12e}, took {:?}", start.elapsed());
}
