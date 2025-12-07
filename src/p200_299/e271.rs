use crate::utils::primes::primality::{mulmod, powmod};

// f_3(n) counts # of cube roots of unity mod n, including 1
// f_3(n) is multiplicative
// if p = 3k+1: f_3(p^e) = 3;
// if p = 3k+2: f_3(p^e) = 1;
// f_3(3) = 1, f_3(3^e), e>=2 = 3
const N: u64 = 13_082_761_331_670_030;
const factors_1mod3: [u64; 6] = [7, 13, 19, 31, 37, 43];
const factors_rest: [u64; 8] = [2, 3, 5, 11, 17, 23, 29, 41];

const roots: [[u64; 3]; 6] = [
    [1, 2, 4],
    [1, 3, 9],
    [1, 7, 11],
    [1, 5, 25],
    [1, 10, 26],
    [1, 6, 36],
];
pub fn main() {
    let mut sum = 0;
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                for l in 0..3 {
                    for m in 0..3 {
                        for n in 0..3 {
                            sum += crt(&[
                                roots[0][i],
                                roots[1][j],
                                roots[2][k],
                                roots[3][l],
                                roots[4][m],
                                roots[5][n],
                            ]);
                        }
                    }
                }
            }
        }
    }
    dbg!(sum - 1);
}
fn crt(a_s: &[u64; 6]) -> u64 {
    let mut solution = 0;
    for (&a, p) in a_s.iter().zip(factors_1mod3) {
        let M_i = N / p;
        let N_i = modinv(M_i, p);
        solution += mulmod(a, mulmod(M_i, N_i, N), N);
        if solution >= N {
            solution -= N;
        }
    }
    for p in factors_rest {
        let M_i = N / p;
        let N_i = modinv(M_i, p);
        solution += mulmod(M_i, N_i, N);
        if solution >= N {
            solution -= N;
        }
    }
    solution
}
const fn modinv(x: u64, p: u64) -> u64 {
    powmod(x, p - 2, p)
}
