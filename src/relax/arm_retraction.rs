use crate::{Parameters, CLPoly};

pub fn arm_retraction(np: usize, index: i32, parameters: &mut Parameters, lpoly: &mut Vec<CLPoly>) {

  let clfpref =
      parameters.ret_pref_0 + (parameters.ret_pref - parameters.ret_pref_0) /
                       (1.0 + (-1.0 * parameters.ret_switch_exponent * (parameters.cur_time.ln())).exp());
  let mut dz;
  if index == 0 { // first call; z=0
    dz = (2.0 * clfpref * parameters.cur_time.sqrt()).sqrt();
    lpoly[np].z = dz;
  } else {
    let z0 = lpoly[np].z;
    if parameters.alpha ==
        1.0 { // Dan Comment - let's make this pow() call conditional as well...
      dz = 0.50 * clfpref * (parameters.cur_time / parameters.phi_eq).sqrt() * parameters.log_dt_mult /
           (z0 * parameters.psi_rept.sqrt());
    } else {
      dz = 0.50 * clfpref * (parameters.cur_time / (parameters.phi_eq.powf(parameters.alpha))).sqrt() * parameters.log_dt_mult /
           (z0 * parameters.psi_rept.sqrt());
    }
    if dz > (0.50 * lpoly[np].z_chain -
              lpoly[np].z) { // Some chains may end up here at the interface
                               // of reptation and disentanglement
      dz = 0.50 * lpoly[np].z_chain - lpoly[np].z;
      lpoly[np].alive = false;
    }
    lpoly[np].z = z0 + dz;

     
  }
}
