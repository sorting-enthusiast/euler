struct DSU<const LEN: usize>([i32; LEN]);
impl<const LEN: usize> DSU<LEN> {
    fn new() -> Self {
        Self([-1; LEN])
    }
    fn find(&mut self, mut a: i32) -> i32 {
        const MASK: i32 = 1 << 31;
        let mut rep = a;
        while (self.0[rep as usize] & MASK) == 0 {
            rep = self.0[rep as usize];
        }
        while a != rep {
            let tmp = self.0[a as usize];
            self.0[a as usize] = rep;
            a = tmp;
        }
        rep
    }
    fn merge(&mut self, mut a: i32, mut b: i32) -> i32 {
        a = self.find(a);
        b = self.find(b);
        if a == b {
            return 0;
        }
        if self.0[a as usize] < self.0[b as usize] {
            core::mem::swap(&mut a, &mut b);
        }
        self.0[b as usize] += self.0[a as usize];
        self.0[a as usize] = b;
        -self.0[b as usize]
    }
}

// MST
pub fn main() {
    let s = std::fs::read_to_string(r"C:\Users\delta\my_coding_projects\euler\0107_network.txt")
        .unwrap();
    let mut adjacency_matrix = [[0; 40]; 40];
    for (i, l) in s.lines().enumerate() {
        for (j, n) in l.split(',').enumerate() {
            if n == "-" {
                continue;
            }
            adjacency_matrix[i][j] = n.parse().unwrap();
        }
    }
    let mut edges = vec![];
    let mut initial_cost = 0;
    for i in 0..40 {
        for j in i + 1..40 {
            if adjacency_matrix[i][j] != 0 {
                edges.push((i, j));
                initial_cost += adjacency_matrix[i][j];
            }
        }
    }
    edges.sort_unstable_by_key(|&(i, j)| adjacency_matrix[i][j]);
    dbg!(edges.len());
    let mut dsu = DSU::<40>::new();
    let mut cost = 0;
    let mut res = vec![];
    for (u, v) in edges {
        let w = adjacency_matrix[u][v];
        let (u, v) = (u as i32, v as i32);
        if dsu.find(u) != dsu.find(v) {
            cost += w;
            res.push((u, v));
            dsu.merge(u, v);
        }
        if res.len() == 40 - 1 {
            break;
        }
    }
    dbg!(initial_cost - cost);
}
