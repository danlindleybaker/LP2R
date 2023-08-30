use scan_fmt::scan_fmt;
use std::fs::read_to_string;
use crate::Parameters;

pub fn get_input(input_file: &str) -> Result<Parameters, String> {
    let file_string = read_file_to_string(input_file).expect("Can't open file");
    let mut first_line: bool = true;
    let mut inp_par = Parameters::new();
    for (i, line) in file_string.lines().enumerate() {
        if i == 1 {
            (inp_par.freq_min, inp_par.freq_max, inp_par.freq_ratio) =
                scan_fmt!(line, "{} {} {}", f32, f32, f32).expect("Line doesn't look right!");
        } else if i == 2 {
            (inp_par.m_kuhn, inp_par.m_e, inp_par.g_0, inp_par.tau_e) =
                scan_fmt!(line, "{} {} {} {}", f32, f32, f32, f32)
                    .expect("Line doesn't look right!");
        } else if i == 3 {
            (inp_par.g_glass, inp_par.tau_glass, inp_par.beta_glass) =
                scan_fmt!(line, "{} {} {} {}", f32, f32, f32)
                    .expect("Line doesn't look right!");
        } else if i == 4 {
            inp_par.number_of_components = scan_fmt!(line, "{}", i32).expect("Line doesn't look right!");
        } else if i > 4 {
            if first_line {
            let (ptype_tmp, w_comp_tmp) = scan_fmt!(line, "{} {}", i32, f32).expect("Line doesn't look right!");
            inp_par.p_type.push(ptype_tmp);
            inp_par.wt_comp.push(w_comp_tmp);
                first_line = false;
            }
            else {
                let (npoly_tmp, mw_tmp, pdi_tmp) = scan_fmt!(line, "{} {} {}", i32, f32, f32).expect("Line doesn't look right!");
                inp_par.n_poly.push(npoly_tmp);
                inp_par.m_w.push(mw_tmp);
                inp_par.pdi.push(pdi_tmp);
                first_line = true;

            }

        }
    }
    inp_par.initialise();
    Ok(inp_par)
}

fn read_file_to_string(input_file: &str) -> Result<String, String> {
    Ok(read_to_string(input_file).unwrap())
}
