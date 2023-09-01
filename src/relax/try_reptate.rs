use crate::{Parameters, CLPoly};

pub fn try_reptate(np: usize, parameters: &mut Parameters, lpoly: &mut Vec<CLPoly>) {
  if !lpoly[np].rept_set {
    let tcomp = parameters.rept_switch_factor * 3.0 * lpoly[np].z_chain *
                   lpoly[np].z * lpoly[np].z * parameters.psi_rept;
    if tcomp < parameters.cur_time {
      lpoly[np].rept_set = true;
      let len_to_rept = lpoly[np].z_chain - 2.0 * lpoly[np].z;
      lpoly[np].tau_d_0 =
          3.0 * lpoly[np].z_chain * len_to_rept * len_to_rept * parameters.psi_rept;

      // For only linear polymer, the following ensures that reptation time for
      // a given chain cannot be lower than that of some shorter chain. May need
      // rethinking when branched polymers are added in future.
      if lpoly[np].tau_d_0 < parameters.last_reptation_time {
        if lpoly[np].z_chain < parameters.last_rept_z {
          lpoly[np].tau_d_0 = parameters.last_reptation_time;
        }
      } else {
        parameters.last_reptation_time = lpoly[np].tau_d_0;
        parameters.last_rept_z = lpoly[np].z_chain;
      }
      lpoly[np].z_rept = len_to_rept;

      if lpoly[np].tau_d_0 < parameters.cur_time {
        lpoly[np].p_max = 1;
      } else {
        lpoly[np].p_max = (lpoly[np].tau_d_0 / parameters.cur_time).sqrt().floor() as i32;
        if (lpoly[np].p_max) % 2 == 0 {
          lpoly[np].p_max = lpoly[np].p_max - 1;
        }
      }
      if lpoly[np].p_max < 1 {
        lpoly[np].p_max = 1;
      }
      lpoly[np].p_next = lpoly[np].p_max;
      lpoly[np].rept_wt = 0.0;
      for i in 1..lpoly[np].p_max {
        lpoly[np].rept_wt += 1.0 / ((i * i) as f64);
      }

    } // Set reptation
  }   // Check for reptation

  if lpoly[np].rept_set {
    let psq = (lpoly[np].p_next * lpoly[np].p_next) as f64;
    if parameters.cur_time > (lpoly[np].tau_d_0 / psq) {
      lpoly[np].z += 0.50 * lpoly[np].z_rept / (lpoly[np].rept_wt * psq);
      if lpoly[np].p_next == 1 {
        lpoly[np].z = 0.50 * lpoly[np].z_chain;
        lpoly[np].alive = false;
      } else {
        lpoly[np].p_next = lpoly[np].p_next - 2;
      }
    }
  }
}
