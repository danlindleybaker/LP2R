use lp2r::prep::init_parameters;

fn main() {
    let input_file: &str = "inp.dat";
    println!("LP2R Rust...");

    // parse input file
    let (parameters, lpoly, data_arrays) = init_parameters(input_file).expect("Problem parsing input file");
    println!("#####################");
    println!("Input Struct:");
    println!("{:?}", parameters);
    println!("#####################");
    println!("Lpoly[0]");
    println!("{:?}", lpoly[0]);
    println!("#####################");
    println!("DataArrays: ");
    println!("{:?}", data_arrays);


    // 
}
