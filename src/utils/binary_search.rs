// Because rust decided not to include lower_bound/upper_bound functions
use core::hint::select_unpredictable;
pub fn upper_bound<T, F, P, U>(a: &[T], value: &U, mut f: F, mut p: P) -> usize
where
    P: FnMut(&T) -> U,
    F: FnMut(&U, &U) -> bool,
{
    let mut size = a.len();
    if size == 0 {
        return 0;
    }
    let mut base = &raw const a[0];

    while size > 1 {
        let half = size >> 1;
        let mid = unsafe { base.add(half) };

        let cmp = !f(value, &p(unsafe { &*(mid.sub(1)) }));

        // Binary search interacts poorly with branch prediction, so force
        // the compiler to use conditional moves if supported by the target
        // architecture.
        base = select_unpredictable(cmp, mid, base);

        size -= half;
    }

    let cmp = !f(value, &p(unsafe { &*base }));
    usize::from(cmp) + unsafe { base.offset_from(&raw const a[0]) as usize }
}

pub fn lower_bound<T, F, P, U>(a: &[T], value: &U, mut f: F, mut p: P) -> usize
where
    P: FnMut(&T) -> U,
    F: FnMut(&U, &U) -> bool,
{
    let mut size = a.len();
    if size == 0 {
        return 0;
    }
    let mut base = &raw const a[0];

    while size > 1 {
        let half = size >> 1;
        let mid = unsafe { base.add(half) };

        let cmp = f(&p(unsafe { &*(mid.sub(1)) }), value);

        // Binary search interacts poorly with branch prediction, so force
        // the compiler to use conditional moves if supported by the target
        // architecture.
        base = select_unpredictable(cmp, mid, base);

        size -= half;
    }

    unsafe { base.offset_from(&raw const a[0]) as usize }
}

pub fn binsearch<T, P, F>(n: usize, mut map: P, mut pred: F) -> usize
where
    P: FnMut(usize) -> T,
    F: FnMut(T) -> bool,
{
    if n == 0 {
        return 0;
    }
    let two_k = 1usize << (usize::BITS - 1 - n.leading_zeros());
    let mut b = if pred(map(n >> 1)) {
        n - two_k
    } else {
        usize::MAX
    };
    let mut bit = two_k >> 1;
    while bit != 0 {
        if pred(map(b + bit)) {
            b += bit;
        }
        bit >>= 1;
    }
    b + 1
}
