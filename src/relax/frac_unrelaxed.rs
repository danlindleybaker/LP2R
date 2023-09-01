use crate::{CLPoly, Parameters, DataArrays};


pub fn frac_unrelaxed(parameters: &mut Parameters, lpoly: &mut Vec<CLPoly>, data_arrays: &mut DataArrays) {

  let deltaphi: f64;
  let dphi_max:f64;
  let phit_sv = parameters.phi_true;

  parameters.phi_true = 0.0;
  for i in 0..parameters.number_of_polymers as usize
    {
    if lpoly[i].alive {
      parameters.phi_true += lpoly[i].wt * (1.0 - 2.0 * lpoly[i].z / lpoly[i].z_chain);
    }
  }

  if parameters.cur_time > parameters.t_cr_start { // only thin tube available before this time
    let a_zeta_inv = (1.0 / ((1.0 - parameters.delta_cr) * (1.0 - parameters.delta_cr))) - 1.0;
    let mut delta_cr_now =
        1.0 - (1.0 - parameters.delta_cr) * (1.0 + a_zeta_inv / parameters.cur_time).sqrt();
    if delta_cr_now < 0.0 {
      delta_cr_now = 0.0;
    }

    if parameters.above_tau_e_first {
      parameters.above_tau_e_first= false;
      deltaphi = 1.0 - parameters.phi_true; // tube relaxation starts now
    } else {
      deltaphi = phit_sv - parameters.phi_true;
    }
    let st_max_drop_now = (1.0 - parameters.st_max_drop) / (1.0 - delta_cr_now * parameters.st_max_drop);
    dphi_max = parameters.phi_st * st_max_drop_now;
    // ********************************************
    if !parameters.supertube_activated {
      if deltaphi <= dphi_max {
        parameters.phi_st = parameters.phi_true;
      } else { // start of new CR Rouse region
        parameters.supertube_activated = true;
        parameters.phi_st_0 = parameters.phi_st - delta_cr_now * deltaphi;
        parameters.st_activ_time = parameters.cur_time;
        parameters.phi_st -= dphi_max; // this much can relax in one time step
      }
    } else { // already in CR Rouse
      let tv = parameters.phi_st_0 * ((parameters.st_activ_time / parameters.cur_time).ln() / (2.0 * parameters.alpha)).exp();
      if tv < parameters.phi_true {
        parameters.phi_st = parameters.phi_true;
        parameters.supertube_activated = false;
      } else {                      // continue in CR Rouse
        if deltaphi <= dphi_max { // normal power-law
          parameters.phi_st = tv;
        } else { // additional CR; add extra drop to parameters.phi_st_0
          parameters.phi_st_0 -= delta_cr_now * deltaphi;

          parameters.phi_st =
              parameters.phi_st_0 * ((parameters.st_activ_time / parameters.cur_time).ln() / (2.0 * parameters.alpha)).exp();
        }
      }
    }
    // ********************************************
  } // parameters.cur_time > t_CR_START condition

  // store data and calculate t_eq, a_eq
  let t_equil_cur_tube = parameters.a_eq * parameters.cur_time * (1.0 + parameters.b_eq / parameters.cur_time.sqrt());
  data_arrays.t_ar.push(parameters.cur_time);
  data_arrays.phi_ar.push(parameters.phi_true);
  data_arrays.phi_st_ar.push(parameters.phi_st);
  data_arrays.t_eq_ar.push(t_equil_cur_tube);
  parameters.phi_eq = get_phi_eq(parameters, data_arrays);
  let mut tpast: f64;
  if !parameters.supertube_activated {
    if parameters.cur_time > 1.0 {
      tpast = -1.0 * parameters.b_eq + (parameters.b_eq * parameters.b_eq + 4.0 * parameters.cur_time / parameters.a_eq).sqrt();
      tpast = 0.25 * tpast * tpast;
      if parameters.phi_eq < 0.999999 {
        let tv = parameters.b_zeta * parameters.phi_eq.powf(3.0 * parameters.alpha) * tpast;
        if tv < parameters.psi_rept {
          parameters.psi_rept = tv;
          parameters.phi_rept = parameters.phi_eq;
        }
      }
    }
  }
}


fn get_phi_eq(parameters: &mut Parameters, data_arrays: &mut DataArrays) -> f64
{
let mut phi_eq1=1.0;
let mut n1=parameters.phi_eq_index;
if parameters.cur_time < 1.0{ // no interpolation
 while  n1 < data_arrays.t_eq_ar.len() && parameters.cur_time > data_arrays.t_eq_ar[n1]{ n1+=1};
 if n1 > 0 {n1=n1-1;}
 parameters.phi_eq_index=n1;
                }
else{
 while  n1 < data_arrays.t_eq_ar.len() && parameters.cur_time > data_arrays.t_eq_ar[n1]{ n1+=1};
 if n1 > 0 {n1=n1-1;}
 parameters.phi_eq_index=n1;
if n1 > 5{ // Avoid problem with having zero as time at the beginning
 phi_eq1=data_arrays.phi_st_ar[n1];
 let deltat=parameters.cur_time - data_arrays.t_eq_ar[n1];
 if deltat > 1.0e-6{ //guard against having negative in log from finite precision
phi_eq1 += (data_arrays.phi_st_ar[n1+1] - data_arrays.phi_st_ar[n1])*(parameters.cur_time/data_arrays.t_eq_ar[n1]).ln()/parameters.log_dt_mult;
                    }
          }
  else{phi_eq1=1.0;}
    }
return phi_eq1;
}
