use scan_fmt::scan_fmt;
use std::fs::read_to_string;

#[derive(Debug)]
pub struct InputParameters {
    pub freq_min: f32,
    pub freq_max: f32,
    pub freq_ratio: f32,
    pub m_kuhn: f32,
    pub m_e: f32,
    pub g_n: f32,
    pub tau_e: f32,
    pub g_inf: f32,
    pub tau_g: f32,
    pub beta_g: f32,
    pub number_of_components: i32,
    pub p_type: Vec<i32>,
    pub w_comp: Vec<f32>, 
    pub n_poly: Vec<i32>,
    pub m_w: Vec<f32>,
    pub pdi: Vec<f32>,

}

impl InputParameters {
    fn new() -> InputParameters {
        InputParameters {
            freq_min: 0.0,
            freq_max: 0.0,
            freq_ratio: 0.0,
            m_kuhn: 0.0,
            m_e: 0.0,
            g_n: 0.0,
            tau_e: 0.0,
            g_inf: 0.0,
            tau_g: 0.0,
            beta_g: 0.0,
            number_of_components: 0,
            p_type: Vec::new(),
            w_comp: Vec::new(),
            n_poly: Vec::new(),
            m_w: Vec::new(),
            pdi: Vec::new(),

        }
    }
}

pub fn get_input(input_file: &str) -> Result<InputParameters, String> {
    let file_string = read_file_to_string(input_file).expect("Can't open file");
    let mut first_line: bool = true;
    let mut inp_par = InputParameters::new();
    for (i, line) in file_string.lines().enumerate() {
        if i == 1 {
            (inp_par.freq_min, inp_par.freq_max, inp_par.freq_ratio) =
                scan_fmt!(line, "{} {} {}", f32, f32, f32).expect("Line doesn't look right!");
        } else if i == 2 {
            (inp_par.m_kuhn, inp_par.m_e, inp_par.g_n, inp_par.tau_e) =
                scan_fmt!(line, "{} {} {} {}", f32, f32, f32, f32)
                    .expect("Line doesn't look right!");
        } else if i == 3 {
            (inp_par.g_inf, inp_par.tau_g, inp_par.beta_g) =
                scan_fmt!(line, "{} {} {} {}", f32, f32, f32)
                    .expect("Line doesn't look right!");
        } else if i == 4 {
            inp_par.number_of_components = scan_fmt!(line, "{}", i32).expect("Line doesn't look right!");
        } else if i > 4 {
            if first_line {
            let (ptype_tmp, w_comp_tmp) = scan_fmt!(line, "{} {}", i32, f32).expect("Line doesn't look right!");
            inp_par.p_type.push(ptype_tmp);
            inp_par.w_comp.push(w_comp_tmp);
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
    Ok(inp_par)
}

fn read_file_to_string(input_file: &str) -> Result<String, String> {
    Ok(read_to_string(input_file).unwrap())
}
