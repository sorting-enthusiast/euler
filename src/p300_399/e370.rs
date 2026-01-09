use rayon::iter::{IntoParallelIterator, ParallelIterator};

const N: usize = 2.5e13 as _;

// https://en.wikipedia.org/wiki/Coprime_integers#Generating_all_coprime_pairs
fn calkin_wilf(x: usize, y: usize) -> usize {
    assert!(x < y);
    let perimeter = x * x + x * y + y * y;
    if perimeter > N {
        return 0;
    }
    (if perimeter - y * y >= y * y {
        N / perimeter
    } else {
        0
    }) + calkin_wilf(x, x + y)
        + calkin_wilf(y, x + y)
}
//1e12:1518637586385
// TODO: optimize the crap out of this
pub fn main() {
    let mut res = N / 3;
    let mut stack = vec![(1, 2)];
    let mut i = 0usize;
    while let Some((x, y)) = stack.pop() {
        let perimeter = x * x + x * y + y * y;
        if perimeter <= N {
            if i.trailing_zeros() >= 31 {
                println!("{i}: {x}, {y}");
            }
            if perimeter >= 2 * y * y {
                res += N / perimeter;
            }
            stack.push((x, x + y));
            stack.push((y, x + y));
            if i == 24 {
                dbg!(stack.len());
                break;
            }
            i += 1;
        }
    }
    res += stack
        .into_par_iter()
        .map(|(x, y)| {
            let mut res = 0;
            let mut stack = vec![(x, y)];
            while let Some((x, y)) = stack.pop() {
                let perimeter = x * x + x * y + y * y;
                if perimeter <= N {
                    if perimeter >= 2 * y * y {
                        res += N / perimeter;
                    }
                    stack.push((x, x + y));
                    stack.push((y, x + y));
                }
            }
            res
        })
        .sum::<usize>();
    dbg!(res);
    //dbg!(12996874312u64);
    dbg!(i);
}
