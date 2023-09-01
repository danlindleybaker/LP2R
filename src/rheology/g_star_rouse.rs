use crate::{Parameters, CLPoly};
pub fn g_star_rouse(freq: f64, grs: &mut f64,  g2rs: &mut f64, ers: &mut f64,
                e2rs: &mut f64, parameters: &mut Parameters, lpoly: &mut Vec<CLPoly>) {
  *grs = 0.0;
  *g2rs = 0.0;
  *ers = 0.0;
  *e2rs = 0.0;

  let mut g_r: f64;
  let mut g2_r: f64;
  let mut e_r: f64;
  let mut e2_r: f64;
  let mut taup: f64;
  let mut tv: f64;
  let mut tv2: f64;
  for  i1 in 0..parameters.number_of_polymers as usize {
    if lpoly[i1].relax_free_rouse {
      (g_r , g2_r , e_r , e2_r) =( 0.0, 0.0 ,0.0, 0.0);
      let tau1 = lpoly[i1].t_f_rouse; // Z^2 * tau_e
      let pmax = (lpoly[i1].z_chain * parameters.n_e).ceil() as i32;

      for p in 1..=pmax {
        let psq = (p * p) as f64;
        taup = tau1 / (2.0 * psq); // stress relaxation time
        tv = freq * taup;
        tv2 = tv * tv;
        g_r += tv2 / (1.0 + tv2);
        g2_r += tv / (1.0 + tv2);
        if p % 2 != 0 {
          taup = tau1 / psq;
          tv = freq * taup;
          tv2 = tv * tv;
          e_r += tv2 / (psq * (1.0 + tv2));
          e2_r += tv / (psq * (1.0 + tv2));
        }
      }
      g_r = g_r * 5.0 * lpoly[i1].wt / (4.0 * lpoly[i1].z_chain);
      g2_r = g2_r * 5.0 * lpoly[i1].wt / (4.0 * lpoly[i1].z_chain);
      *grs += g_r;
      *g2rs += g2_r;
      *ers += lpoly[i1].wt * e_r;
      *e2rs += lpoly[i1].wt * e2_r;
    }
  }
}
