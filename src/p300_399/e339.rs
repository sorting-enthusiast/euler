use std::collections::HashMap;

pub fn main() {
    const N: i32 = 5;
    let mut cache = HashMap::new();
    dbg!(opt(N, N, &mut cache));
}
fn opt(w: i32, b: i32, cache: &mut HashMap<(i32, i32), f32>) -> f32 {
    if w < 0 || b < 0 {
        return 0.;
    }
    if w == 0 {
        cache.insert((w, b), b as f32);
        return b as f32;
    }
    dbg!((w, b));
    if let Some(&v) = cache.get(&(w, b)) {
        return v;
    }
    let remove_some = (0..b).map(|i| opt(i, b, cache)).max_by(f32::total_cmp);
    let n = w + b;
    let n = n as f32;
    let do_nothing =
        (w as f32 / n) * opt(w + 1, b - 1, cache) + (b as f32 / n) * opt(w - 1, b + 1, cache);
    let ret = remove_some.map_or(do_nothing, |v| do_nothing.max(v));
    cache.insert((w, b), ret);
    ret
}
