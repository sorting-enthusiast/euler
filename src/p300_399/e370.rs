const N: usize = 1e6 as _;

// https://en.wikipedia.org/wiki/Coprime_integers#Generating_all_coprime_pairs
fn calkin_wilf(x: usize, y: usize) -> usize {
    assert!(x < y);
    let perimeter = x * x + x * y + y * y;
    if perimeter > N {
        return 0;
    }
    (if perimeter - y * y >= y * y {
        (N / perimeter).isqrt()
    } else {
        0
    }) + calkin_wilf(x, x + y)
        + calkin_wilf(y, x + y)
}

pub fn main() {
    dbg!((N / 3).isqrt() + calkin_wilf(1, 2));
}
