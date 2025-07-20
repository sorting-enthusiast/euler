use crate::utils::powerful_numbers::PowerfulExt;

// first 100%!!!
// (p^e)' = e*p^(e-1): easy enough to prove, also shown on wiki page
// let p be a prime not dividing n:
// gcd(p^e*n, (p^e*n)') = gcd(p^e*n, e*p^(e-1)*n + p^e * n') = p^(e-1) * gcd(pn, en + pn')
// gcd(pn, en + pn'):
// p|e: gcd(pn, en + pn') = gcd(pn, p*k*n + pn') = p * gcd(n, k*n+n') = p * gcd(n,n') => gcd(p^e*n, (p^e*n)') = p^e * gcd(n,n')
// p!|e: gcd(pn, en + pn') = gcd(n, n')

// let v_p(n) be the p-adic valuation of n.
// define f(n) = gcd(n,n'): product over p|n of {p^v_p(n) if p|v_p(n), else p^(v_p(n)-1)}
// notice that f(n) is multiplicative, and that f(p) = gcd(p,1) = 1 = u(p)
// find h such that f = h * u <=> h = f * mobius
// h(p^e) = f(p^e) - f(p^(e-1))
// use powerful number trick to evaluate the sum of f(n) as the sum over powerful numbers k of h(k)*(sum up to N/k of u(i)) = h(k) * (N/k)
// subtract 1 because the qustion asked to sum f starting at n=2, not 1.
pub fn main() {
    const N: i64 = 5e15 as _;
    let start = std::time::Instant::now();
    let mut sum = 0;
    let h = |p: i64, e| {
        let pe2 = p.pow(e as u32 - 2);
        (if e % p == 0 { pe2 * p * p } else { pe2 * p })
            - if (e - 1) % p == 0 { pe2 * p } else { pe2 }
    };
    for (n, hn) in PowerfulExt::<_, { i64::MAX }>::new(N, h).filter(|&(_, hn)| hn != 0) {
        sum += hn * (N / n);
    }
    sum -= 1; // sum requested starts at 2, not 1
    let end = start.elapsed();
    println!("res = {sum}, took {end:?}");
}
