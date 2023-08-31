use lp2r::prep::init_parameters;

fn main() {
    let input_file: &str = "inp.dat";
    println!("LP2R Rust...");

    // parse input file
    let test = init_parameters(input_file).expect("Problem parsing input file");
    println!("Input Struct:");
    println!("{:?}", test);

    // 
}
