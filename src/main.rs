use lp2r::prep::init_parameters;
use lp2r::relax::time_step;
use lp2r::rheology::lin_rheology;

fn main() {
    let input_file: &str = "inp.dat";
    println!("LP2R Rust...");

    // parse input file
    let (mut parameters, mut lpoly, mut data_arrays) =
        init_parameters(input_file).expect("Problem parsing input file");
    println!("#####################");
    println!("Input Struct:");
    println!("{:?}", parameters);
    println!("#####################");
    println!("Lpoly[0]");
    println!("{:?}", lpoly[0]);
    println!("#####################");
    println!("DataArrays: ");
    println!("{:?}", data_arrays);

    println!("#####################");
    println!("Entangled dynamics = {}", parameters.entangled_dynamics);

    let mut nalive: i32;
    if parameters.entangled_dynamics {
        nalive = time_step(0, &mut parameters, &mut lpoly, &mut data_arrays);

        while nalive > 0 {
            nalive = time_step(1, &mut parameters, &mut lpoly, &mut data_arrays);

        }
    }

    println!("#####################");
    let results = lin_rheology(&mut parameters, &mut lpoly, &mut data_arrays);
    println!("Results: {:?}", results); 
    //
}
