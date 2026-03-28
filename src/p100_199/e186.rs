pub fn main() {
    let mut buf = [0; 55];
    let mut k = 0;
    let mut s_k = std::iter::from_fn(move || {
        k += 1;
        let val = if k <= 55 {
            let k = k as i64;
            (100003 - 200003 * k + 300007 * k * k * k).rem_euclid(1_000_000)
        } else {
            let s_k_24 = buf[(k - 24 - 1) % 55];
            let s_k_55 = buf[(k - 55 - 1) % 55];
            (s_k_24 + s_k_55) % 1_000_000
        };

        buf[(k - 1) % 55] = val;
        Some(val)
    });

    let mut components = DSU::<1_000_000>::new();
    let mut calls = 0usize;
    loop {
        let caller = unsafe { s_k.next().unwrap_unchecked() };
        let called = unsafe { s_k.next().unwrap_unchecked() };
        if caller == called {
            continue;
        }
        calls += 1;
        if components.merge(caller as _, called as _) != 0
            && components.get_size(524_287) >= 990_000
        {
            break;
        }
    }
    dbg!(calls);
}
struct DSU<const LEN: usize>(Box<[i32; LEN]>);
impl<const LEN: usize> DSU<LEN> {
    fn new() -> Self {
        Self(Box::new([-1; LEN]))
    }
    const fn find(&mut self, mut a: i32) -> i32 {
        let mut rep = a;
        while self.0[rep as usize] >= 0 {
            rep = self.0[rep as usize];
        }
        while a != rep {
            let tmp = self.0[a as usize];
            self.0[a as usize] = rep;
            a = tmp;
        }
        rep
    }
    const fn get_size(&mut self, a: i32) -> i32 {
        -self.0[self.find(a) as usize]
    }
    const fn merge(&mut self, mut a: i32, mut b: i32) -> i32 {
        a = self.find(a);
        b = self.find(b);
        if a == b {
            return 0;
        }
        if self.0[a as usize] > self.0[b as usize] {
            core::mem::swap(&mut a, &mut b);
        }
        self.0[b as usize] += self.0[a as usize];
        self.0[a as usize] = b;
        -self.0[b as usize]
    }
}
