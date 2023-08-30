use lp2r::get_input::get_input;

fn main() {
    let input_file: &str = "inp.dat";
    println!("LP2R Rust...");

    // parse input file
    let test = get_input(input_file).expect("Problem parsing input file");
    println!("Input Struct:");
    println!("{:?}", test);

    // 
}
