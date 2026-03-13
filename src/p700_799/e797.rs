// cyclogenic polynomials can be constructed from cyclotomic polynomials:
// the 10 6-cyclogenic polynomials are \Phi_6(x) multiplied by each by the product of all subsets of {\Phi_1(x),\Phi_2(x),\Phi_3(x)}
// plus \Phi_2(x)\Phi_3(x) and the same thing multiplied by \Phi_1(x)
// concretely, P_n(x) = sum over all subsets a_i of the divisors of n s.t. lcm(a_i) = n of \prod_{k\in a_i}\Phi_k(x)
// (P_n)'(x) = \sum_{d|n}{P_d(x)} = \prod_{d|n}(1+\Phi_d(x))
// \Phi_d(2) can be computed using a sieve, and so can (P_n)'(2)
// use mobius inversion to get P_n(2)
pub fn main() {}
