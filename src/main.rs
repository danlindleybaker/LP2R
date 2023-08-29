use lp2r::get_input::get_input;

fn main() {
    let input_file: &str = "inp.dat";
    println!("LP2R Rust...");
    let test = get_input(input_file).unwrap();
    println!("Input Struct:");
    println!("{:?}", test);
}
