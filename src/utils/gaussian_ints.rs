#[derive(Clone, Copy, Debug, Default)]
pub struct GaussianI64 {
    pub re: i64,
    pub im: i64,
}
impl GaussianI64 {
    #[must_use]
    pub const fn conj(&self) -> Self {
        Self {
            re: self.re,
            im: -self.im,
        }
    }
    #[must_use]
    pub const fn norm(&self) -> u64 {
        let a = self.re.unsigned_abs();
        let b = self.im.unsigned_abs();
        a * a + b * b
    }
}

impl std::ops::Add for GaussianI64 {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        Self {
            re: self.re + rhs.re,
            im: self.im + rhs.im,
        }
    }
}
impl std::ops::Mul for GaussianI64 {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        let Self { re: a, im: b } = self;
        let Self { re: c, im: d } = rhs;
        Self {
            re: a * c - b * d,
            im: a * d + b * c,
        }
    }
}
impl std::ops::Div for GaussianI64 {
    type Output = Self;
    fn div(self, rhs: Self) -> Self::Output {
        let norm = rhs.norm();
        let Self { re: a, im: b } = self * rhs.conj();
        Self {
            re: a / norm as i64,
            im: b / norm as i64,
        }
    }
}

impl From<(i64, i64)> for GaussianI64 {
    fn from(value: (i64, i64)) -> Self {
        Self {
            re: value.0,
            im: value.1,
        }
    }
}
