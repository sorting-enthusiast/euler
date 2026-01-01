/// Bivariate OGF coefficient extraction
/// Any single grouping of B and W (i,j) where at least one of i,j != 0,
/// can appear any number of times. The contribution of (i,j) is therefore 1 + x^i y^j + (x^i y^j)^2 + ... = (1 - x^i y^j)^-1
/// Then, the OGF of the number of groupings of i B's and j W's F(x, y) is the product over i,j of (1 - x^i y^j)^-1
/// The answer is [x^60 y^40]F(x, y)
pub fn main() {
    const B: usize = 60;
    const W: usize = 40;

    let start = std::time::Instant::now();
    let mut gf = vec![[0; W + 1]; B + 1].into_boxed_slice();
    gf[0][0] = 1u64;

    for i in 0..=B {
        for j in 0..=W {
            if i | j == 0 {
                continue;
            }
            for b in i..=B {
                for w in j..=W {
                    gf[b][w] += gf[b - i][w - j];
                }
            }
        }
    }
    println!("res = {}, took {:?}", gf[B][W], start.elapsed());
}
