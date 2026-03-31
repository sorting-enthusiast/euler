pub fn main() {
    const N: usize = 1e17 as _;
    let start = std::time::Instant::now();
    let mut fibs = [0; 82];
    let mut dfs_fibs = [0; 82];

    fibs[0] = 1;
    fibs[1] = 2;
    dfs_fibs[1] = 1;

    for i in 2..82 {
        fibs[i] = fibs[i - 1] + fibs[i - 2];
        dfs_fibs[i] = fibs[i - 2] + dfs_fibs[i - 1] + dfs_fibs[i - 2];
    }
    fibs.reverse();
    dfs_fibs.reverse();

    let res = dfs(N, &fibs, &dfs_fibs);
    println!("res = {res}, took {:?}", start.elapsed());
}
fn dfs(n: usize, fibs: &[usize], dfs_fibs: &[usize]) -> usize {
    if fibs.is_empty() {
        return 0;
    }
    let p = fibs.partition_point(|&f| f > n);
    if p == fibs.len() {
        return 0;
    }
    (n - fibs[p]) + dfs_fibs[p] + dfs(n - fibs[p], &fibs[p + 1..], &dfs_fibs[p + 1..])
}
