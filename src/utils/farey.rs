// https://en.wikipedia.org/wiki/Farey_sequence#Next_term
pub fn farey(n: u32) -> impl Iterator<Item = (u32, u32)> {
    let (mut a, mut b, mut c, mut d) = (0, 1, 1, n);
    core::iter::once((a, b)).chain(core::iter::from_fn(move || {
        (c <= n).then(|| {
            let k = (n + b) / d;
            (a, b, c, d) = (c, d, k * c - a, k * d - b);
            (a, b)
        })
    }))
}
