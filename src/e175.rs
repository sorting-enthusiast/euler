use crate::e169::aviv as f;

// f is the function from e169
pub fn main() {
    for i in 1..=256 {
        println!("{i}, {}/{}", f(i), f(i - 1));
    }
    todo!()
}
