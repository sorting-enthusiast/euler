// counts all numbers in range s.t. all their prime factors are from the prime set given
pub fn sum_over_smooth(
    acc: u64,
    f: &mut impl FnMut(u64) -> u64,
    lo: u64,
    hi: u64,
    primes: &[u64],
) -> u64 {
    let mut sum = if acc >= lo { f(acc) } else { 0 };
    for (i, &p) in primes.iter().enumerate() {
        if p * acc > hi {
            break;
        }
        let mut power = p;
        loop {
            sum += sum_over_smooth(acc * power, f, lo, hi, &primes[i + 1..]);
            power *= p;
            if power * acc > hi {
                break;
            }
        }
    }
    sum
}
