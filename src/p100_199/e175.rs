// f is the function from e169
// recurrence can be extracted from the generating function of f:
// let R_n = \frac{f(n)}{f(n-1)}. Then R_{2n} = 1 + R_n, R_{2n+1} = \frac{R_n}{1 + R_n}
// therefore if R_n = \frac{p}{q} > 1, n = 2k and R_k = \frac{p - q}{q}
// otherwise n = 2k + 1 and R_k = \frac{p}{q - p}
// using this the SBR of n can be gleaned from the continued fraction representation of \frac{p}{q}
// for 123456789/987654321 the continued fraction is [0; 8, 13717421] so the solution is 1,13717420,8
pub fn main() {
    println!("1,13717420,8");
}
