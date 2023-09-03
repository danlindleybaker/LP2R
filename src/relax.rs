use crate::relax::arm_retraction::arm_retraction;
use crate::relax::frac_unrelaxed::frac_unrelaxed;
use crate::relax::try_reptate::try_reptate;
use crate::{CLPoly, DataArrays, Parameters};

pub mod arm_retraction;
pub mod frac_unrelaxed;
pub mod try_reptate;
pub fn time_step(
    index: i32,
    parameters: &mut Parameters,
    lpoly: &mut Vec<CLPoly>,
    data_arrays: &mut DataArrays,
) -> i32 {
    if index != 0 {
        parameters.cur_time *= parameters.dt_mult;
    }

    for i in 0..parameters.number_of_polymers as usize {
        if lpoly[i].alive {
            if !lpoly[i].rept_set {
                arm_retraction(i, index, parameters, lpoly);
            }

            if lpoly[i].alive && parameters.cur_time > 1.0 {
                try_reptate(i, parameters, lpoly);
            }
        }

        if lpoly[i].alive {
            if lpoly[i].z_chain * parameters.phi_st.powf(parameters.alpha)
                < parameters.disentanglement_switch
            {
                lpoly[i].z = 0.50 * lpoly[i].z_chain;
                lpoly[i].alive = false;
                if i == 3 {
                //    println!("z (polymer 4) = {}", lpoly[i].z);
                }
            }
        }
    }

    frac_unrelaxed(parameters, lpoly, data_arrays);

    let mut nalive: i32 = 0;
    for i in 0..parameters.number_of_polymers as usize {
        if lpoly[i].alive {
            nalive += 1;
        }
    }
    return nalive;
}
