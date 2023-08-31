use crate::{CLPoly, DataArrays};
use crate::{prep::get_input::get_input, Parameters};

use self::gen_poly_lin::gen_linlog_normal;
pub fn init_parameters(input_file: &str) -> Result<(Parameters, Vec<CLPoly>, DataArrays), String> {
    let mut parameters = get_input(input_file).expect("problem reading input file");
    let mut lpoly: Vec<CLPoly> = Vec::new();
    for i in 0..parameters.number_of_components as usize {
        gen_linlog_normal(
            parameters.n_poly[i],
            parameters.m_w[i],
            parameters.pdi[i],
            parameters.wt_comp[i],
            &mut lpoly,
            &mut parameters,
        );
    }
    parameters.wt_tot = 0.0;
    for i in 0..parameters.number_of_polymers as usize {
        parameters.wt_tot += lpoly[i].wt;
    }

    for i in 0..parameters.number_of_polymers as usize {
        lpoly[i].wt /= parameters.wt_tot;
    }

    let mut num_entangled: i32 = 0;

    for i in 0..parameters.number_of_polymers as usize {
        if lpoly[i].z_chain < parameters.rouse_switch_factor {
            lpoly[i].relax_free_rouse = true;
            lpoly[i].alive = false;
            parameters.rouse_wt += lpoly[i].wt;
        } else {
            num_entangled += 1;
        }
        parameters.sys_mn += lpoly[i].wt / lpoly[i].mass;
        parameters.sys_mw += lpoly[i].wt * lpoly[i].mass;
    }

    parameters.sys_mn = 1.0 / parameters.sys_mn;
    parameters.sys_pdi = parameters.sys_mw / parameters.sys_mn;

    if num_entangled == 0 {
        parameters.entangled_dynamics = false;
    }

    let data_arrays = DataArrays::new();

    Ok((parameters, lpoly, data_arrays))
}

pub mod gen_poly_lin;
pub mod get_input;
