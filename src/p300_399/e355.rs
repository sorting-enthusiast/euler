use crate::utils::binary_search::binsearch;
use crate::utils::sieve_of_pritchard::sift as sieve;
use std::{
    collections::{BTreeMap, BTreeSet, HashMap, HashSet},
    fmt::Write,
    time::Instant,
};
fn lp(lim: u64, base: u64) -> u64 {
    base.pow((lim as f64).log(base as f64).floor() as u32)
}

fn find_factors(mut n: u64, primes: &[u64]) -> ([u32; 6], usize) {
    let mut factors = [0u32; 6];
    let mut i = 0;
    for &p in primes.iter().take_while(move |&&p| p * p <= n) {
        if n == 1 {
            return (factors, i);
        }
        if n % p == 0 {
            factors[i] = p as u32;
            n /= p;
            i += 1;
            while n % p == 0 {
                n /= p;
            }
        }
    }
    if n > 1 {
        factors[i] = n as u32;
        i += 1;
    }
    (factors, i)
}

#[derive(Debug, Clone, Copy)]
struct Candidate {
    weight: u32,
    cliques: [u32; 4],
    cliques_count: u8,
}
impl Candidate {
    fn new(weight: u32, cliques: [u32; 4], cliques_count: usize) -> Self {
        Self {
            weight,
            cliques,
            cliques_count: cliques_count as u8,
        }
    }
}

#[derive(Debug)]
enum State {
    In,
    DependsOn(Vec<u32>),
}

fn resolve_uncertainty(solution: &mut HashMap<u32, State>, node: u32) {
    let Some(state) = solution.remove(&node) else {
        return;
    };

    if let State::DependsOn(others) = &state {
        for &other in others {
            resolve_uncertainty(solution, other);
        }

        if others.iter().all(|d| !solution.contains_key(d)) {
            solution.insert(node, State::In);
        }
    } else {
        solution.insert(node, State::In);
    }
}

fn reduce_graph(
    candidates: &mut BTreeMap<u32, Candidate>,
    primes_to_candidates: &mut HashMap<u32, Vec<u32>>,
) -> HashMap<u32, State> {
    let mut solution = HashMap::new();

    let decrement_clique = |c,
                            p,
                            candidates: &mut BTreeMap<u32, Candidate>,
                            solution: &mut HashMap<u32, State>|
     -> bool {
        let node = candidates.get_mut(&c).unwrap();
        let old = node.cliques;
        let old_count = node.cliques_count;
        node.cliques = [0; 4];
        node.cliques_count = 0;
        for clique in old
            .into_iter()
            .take(old_count as _)
            .filter(|&clique| clique != p)
        {
            node.cliques[node.cliques_count as usize] = clique;
            node.cliques_count += 1;
        }
        if node.cliques_count == 0 {
            candidates.remove(&c);
            solution.insert(c, State::In);
            true
        } else {
            false
        }
    };

    primes_to_candidates.retain(|&p, cs| {
        if cs.len() == 1 {
            decrement_clique(cs[0], p, candidates, &mut solution);
            false
        } else {
            true
        }
    });

    dbg!(primes_to_candidates.len());

    dbg!(candidates.len());

    // easy to do since cliques are based exclusively on prime factors
    // TODO: implement isolated weight transfer and isolated vertex removal, until further reduction is impossible, then switch to ILP
    let mut valid_cliques = primes_to_candidates
        .keys()
        .copied()
        .collect::<BTreeSet<_>>();

    while !candidates.is_empty() {
        let len_before = candidates.len();

        // isolated weight transfer and isolated vertex removal, domination reduction proved to be of little use
        valid_cliques.retain(|&p| {
            let cs = primes_to_candidates.get_mut(&p).unwrap();
            if cs.len() == 1 {
                decrement_clique(cs[0], p, candidates, &mut solution);
                return false;
            }

            if let Some((&c, &node)) = cs
                .iter()
                .map(|c| (c, candidates.get(c).unwrap()))
                .filter(|&(_, c)| c.cliques_count == 1)
                .max_by_key(|&(_, c)| c.weight)
            {
                let cs_to_erase = cs
                    .extract_if(.., |&mut other| {
                        candidates.get(&other).unwrap().weight <= node.weight
                    })
                    .filter(|&other| other != c)
                    .collect::<Vec<_>>();
                let keep_clique = if cs.is_empty() {
                    solution.insert(c, State::In);
                    candidates.remove(&c);
                    // delete clique
                    false
                } else {
                    solution.insert(c, State::DependsOn(cs.clone()));
                    candidates.remove(&c);
                    for c in cs.iter() {
                        candidates.get_mut(c).unwrap().weight -= node.weight;
                    }
                    if cs.len() == 1 {
                        // decrement clique count, delete clique
                        decrement_clique(cs[0], p, candidates, &mut solution);
                        false
                    } else {
                        true
                    }
                };

                for other in cs_to_erase {
                    let other_node = candidates.remove(&other).unwrap();

                    for &other_clique in other_node
                        .cliques
                        .iter()
                        .take(other_node.cliques_count as usize)
                        .filter(|&&other_p| p != other_p)
                    {
                        primes_to_candidates
                            .entry(other_clique)
                            .and_modify(|cs| cs.retain(|&c| c != other));
                    }
                }
                keep_clique
            } else {
                true
            }
        });
        valid_cliques.retain(|p| primes_to_candidates.get(p).is_none_or(|v| !v.is_empty()));
        primes_to_candidates.retain(|p, _| valid_cliques.contains(p));

        let len_after = candidates.len();
        dbg!(len_after);
        if len_before == len_after {
            break;
        }
    }
    solution
}

// 355:
// TODO: implement improved domination filter
pub fn main() {
    let start = Instant::now();
    let n = 200_000;
    let primes = sieve(n);

    let j = binsearch(dbg!(primes.len()), |i| primes[i].pow(2), |p2| p2 <= n);

    let ps_to_powers = primes[..j]
        .iter()
        .map(|&p| (p as u32, lp(n, p) as u32))
        .collect::<HashMap<_, _>>();

    let mut sum = dbg!(
        1 + primes[..j]
            .iter()
            .map(|&p| u64::from(ps_to_powers.get(&(p as u32)).copied().unwrap()))
            .sum::<u64>()
            + primes[j..].iter().sum::<u64>()
    );

    let last_prime = binsearch(primes[j..].len(), |i| primes[j..][i] << 1, |p| p <= n);
    /* assert_eq!(
        last_prime,
        upper_bound(&primes[j..], &n, u64::lt, |&e| e << 1)
    ); */
    println!("started candidate sifting");
    let mut candidates = (2..=n)
        .filter(|x| primes[last_prime..].binary_search(x).is_err()) // filter primes
        .map(|i| (i as u32, find_factors(i, &primes[..j])))
        .filter(|&(_, (_, f))| f > 1) // filter prime powers
        .filter_map(|(i, (pfs, count))| {
            let penalty = pfs
                .into_iter()
                .take(count)
                .fold(0, |acc, pf| acc + ps_to_powers.get(&pf).unwrap_or(&pf));
            (i > penalty).then(|| (i, i - penalty, pfs, count))
        })
        .map(|(i, weight, pfs, count)| {
            (
                i,
                Candidate::new(weight, [pfs[0], pfs[1], pfs[2], pfs[3]], count),
            )
        })
        .collect::<BTreeMap<_, _>>();
    println!("stopped candidate sifting");
    let mut primes_to_candidates: HashMap<u32, Vec<u32>> = HashMap::new();
    for (&n, c) in &candidates {
        c.cliques
            .into_iter()
            .take(c.cliques_count as _)
            .for_each(|p| {
                primes_to_candidates.entry(p).or_default().push(n);
            });
    }

    let mut solution = reduce_graph(&mut candidates, &mut primes_to_candidates);

    dbg!(solution.values().filter(|e| matches!(e, State::In)).count());
    dbg!(
        solution
            .values()
            .filter(|e| matches!(e, State::DependsOn(_)))
            .count()
    );
    dbg!(candidates.values().all(|info| info.cliques_count == 2));

    let map_to_index = candidates
        .keys()
        .enumerate()
        .map(|(i, &c)| (c, i))
        .collect::<HashMap<_, _>>();
    /*
       let mut nodes = String::new();
       for w in candidates.values().map(|info| info.weight) {
           writeln!(&mut nodes, "{w}").unwrap();
       }
       std::fs::write(
           r"C:\Users\delta\my_coding_projects\euler_python\nodes.txt",
           nodes,
       )
       .unwrap();
       let mut edges = String::new();
       for list in primes_to_candidates.values() {
           let len = list.len();
           for i in 0..len - 1 {
               let &a = map_to_index.get(&list[i]).unwrap();
               for j in &list[i + 1..] {
                   let &b = map_to_index.get(j).unwrap();
                   writeln!(&mut edges, "{a} {b}").unwrap();
               }
           }
       }
       std::fs::write(
           r"C:\Users\delta\my_coding_projects\euler_python\edges.txt",
           edges,
       )
       .unwrap();
    */
    let new_in: HashSet<usize> = HashSet::from([
        654, 737, 806, 812, 816, 824, 826, 827, 829, 831, 835, 839, 842, 843, 844, 846, 847, 848,
        850, 852, 853, 854, 855, 856, 857, 861, 862, 863, 864, 865, 866, 867, 868, 869, 870, 871,
    ]);
    /*
    // 1_000_000
    let new_in: HashSet<usize> = HashSet::from([
        2190, 2300, 2381, 2426, 2450, 2451, 2467, 2471, 2484, 2499, 2501, 2506, 2516, 2526, 2528,
        2533, 2538, 2541, 2545, 2547, 2548, 2551, 2554, 2556, 2562, 2563, 2567, 2568, 2572, 2574,
        2575, 2576, 2577, 2578, 2580, 2581, 2582, 2583, 2584, 2585, 2586, 2587, 2588, 2589, 2592,
        2593, 2594, 2595, 2596, 2597, 2598, 2599, 2600, 2601, 2602, 2603, 2604, 2605, 2606, 2607,
        2608, 2609, 2610,
    ]); */
    solution.extend(
        map_to_index
            .into_iter()
            .filter_map(|(c, i)| (new_in.contains(&i)).then_some((c, State::In))),
    );
    dbg!(solution.values().filter(|e| matches!(e, State::In)).count());
    dbg!(
        solution
            .values()
            .filter(|e| matches!(e, State::DependsOn(_)))
            .count()
    );
    while let Some((&c, _)) = solution
        .iter()
        .find(|(_, s)| matches!(s, State::DependsOn(_)))
    {
        resolve_uncertainty(&mut solution, c);
    }
    dbg!(&solution);
    for c in solution.into_keys().map(u64::from) {
        let (factors, cnt) = find_factors(c, &primes);
        for &p in factors.iter().take(cnt) {
            sum -= u64::from(ps_to_powers.get(&p).copied().unwrap_or(p));
        }
        sum += c;
    }

    let end = start.elapsed();
    println!("sum = {sum}, took {end:?}");
}
