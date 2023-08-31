use crate::{prep::get_input::get_input, Parameters};

pub fn init_parameters(input_file: &str) -> Result<Parameters, String> {
    let parameters = get_input(input_file).expect("problem reading input file");
    Ok(parameters)
}

pub mod gen_poly_lin;
pub mod get_input;
