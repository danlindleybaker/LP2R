use std::fs::read_to_string;
use scan_fmt::scan_fmt;
pub fn get_input(input_file: &str) -> Result<String, String> {
    let file_string = read_file_to_string(input_file)?;
    Ok(file_string)
}


fn read_file_to_string(input_file: &str) -> Result<String, String> {
    Ok(read_to_string(input_file).unwrap())
}
