use super::utils::binary_search::binsearch;
use super::utils::sieve_of_pritchard::sift;
// 501
// ez: only look at composites with either 3 distinct prime factors,
// or 2 distinct prime factors, with one of them being of multplicity 3 and the other of 1,
// or p^7
// TODO: fix
pub fn main() {
    const N: u64 = 1e12 as u64;
    let block_size = (N as f64).sqrt() as u64 + 1;
    let primes = sift(block_size);
    let len = primes.len();
    let mut count = primes.iter().take_while(|&&p| p.pow(7) <= N).count();

    for (i, p) in primes.iter().enumerate() {
        let p_cubed = p.pow(3);
        let k = binsearch(len, |i| primes[i] * p_cubed, |e| e <= N);
        //assert_eq!(k, upper_bound(&primes, &N, u64::lt, |&p| p * p_cubed));
        if k == 0 {
            break;
        }
        count += k - usize::from(k >= i);
    }
    for (i, p1) in primes.iter().enumerate() {
        for j in i + 1..len {
            let p2 = primes[j];
            let acc = p1 * p2;
            let k = binsearch(
                primes[j + 1..].len(),
                |i| primes[j + 1..][i] * acc,
                |e| e <= N,
            );
            //assert_eq!(k, upper_bound(&primes[j + 1..], &N, u64::lt, |&p| p * acc));
            if k == 0 {
                break;
            }
            count += k;
        }
    }
    let mut sieve = vec![0u64; 1 + ((block_size as usize + 63) >> 7)].into_boxed_slice();
    let mut low = block_size;
    let mut high = low + block_size;
    let mut primes_in_block = vec![];
    let mut i = 0u64;
    while low < N {
        if i.trailing_zeros() >= 12 {
            println!("{}", i >> 12);
        }
        i += 1;
        if high > N {
            high = N;
        }
        low |= 1;
        for &p in &primes[1..] {
            let mut multiple = low.next_multiple_of(p);
            if multiple & 1 == 0 {
                multiple += p;
            }
            while multiple < high {
                sieve[(multiple - low) as usize >> 7] |= 1 << ((multiple >> 1) & 63);
                multiple += 2 * p;
            }
        }

        primes_in_block.extend(
            (low..high)
                .step_by(2)
                .filter(|n| sieve[(n - low) as usize >> 7] & 1 << ((n >> 1) & 63) == 0),
        );
        if !primes_in_block.is_empty() {
            for p in &primes {
                let p_cubed = p.pow(3);
                let k = binsearch(
                    primes_in_block.len(),
                    |i| primes_in_block[i] * p_cubed,
                    |e| e <= N,
                );
                /* assert_eq!(
                    k,
                    upper_bound(&primes_in_block, &N, u64::lt, |&p| p * p_cubed)
                ); */
                if k == 0 {
                    break;
                }
                count += k;
            }
            for (i, p1) in primes.iter().enumerate() {
                if p1 * p1 * primes_in_block[0] > N {
                    break;
                }
                for p2 in &primes[i + 1..] {
                    let acc = p1 * p2;
                    let k = binsearch(
                        primes_in_block.len(),
                        |i| primes_in_block[i] * acc,
                        |e| e <= N,
                    );
                    //assert_eq!(k, upper_bound(&primes_in_block, &N, u64::lt, |&p| p * acc));
                    if k == 0 {
                        break;
                    }
                    count += k;
                }
            }
            primes_in_block.clear();
        }
        low = high;
        high += block_size;
        sieve.fill(0);
    }

    dbg!(count);
}
