const ELEM_BIT_WIDTH: u8 = (7 - (usize::BITS as u8).leading_zeros()) as u8;
const MASK: usize = usize::BITS as usize - 1;
pub struct BitArray {
    pub bits: Box<[usize]>,
}
impl BitArray {
    pub fn zeroed(n: usize) -> BitArray {
        BitArray {
            bits: vec![0; (n >> ELEM_BIT_WIDTH) + 1].into_boxed_slice(),
        }
    }
    #[inline]
    pub fn zero(&mut self) {
        self.bits.iter_mut().for_each(|x| *x = 0);
    }
    #[inline(always)]
    pub fn set(&mut self, i: usize) {
        self.bits[i >> ELEM_BIT_WIDTH] |= 1 << (i & MASK);
    }
    #[inline(always)]
    pub fn get(&self, i: usize) -> bool {
        self.bits[i >> ELEM_BIT_WIDTH] & (1 << (i & MASK)) != 0
    }
    #[inline(always)]
    pub fn clear(&mut self, i: usize) {
        self.bits[i >> ELEM_BIT_WIDTH] &= !(1 << (i & MASK));
    }
    #[inline(always)]
    pub fn flip(&mut self, i: usize) {
        self.bits[i >> ELEM_BIT_WIDTH] ^= 1 << (i & MASK);
    }
}
