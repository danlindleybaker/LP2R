use lp2r::prep::init_parameters;
use lp2r::relax::time_step;
use lp2r::rheology::lin_rheology;

fn main() {
    let input_file: &str = "inp.dat";
    println!("LP2R Rust...");

    // parse input file
    let (mut parameters, mut lpoly, mut data_arrays) =
        init_parameters(input_file).expect("Problem parsing input file");

    let mut nalive: i32;
    let mut count = 1;
    if parameters.entangled_dynamics {
        nalive = time_step(0, &mut parameters, &mut lpoly, &mut data_arrays);

        while nalive > 0 {
            nalive = time_step(1, &mut parameters, &mut lpoly, &mut data_arrays);
            if count == 609 {
         //       println!("phi_st_ar[69] = {}", data_arrays.phi_st_ar[count]);
            }
//            println!("Count = {}", count);
            count += 1;
        }
    }

    let _results = lin_rheology(&mut parameters, &mut lpoly, &mut data_arrays);

    //
}
