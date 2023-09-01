use crate::{Parameters, kwws, kwwc};

pub fn g_star_glass(w: f64, gp: &mut f64, g2p: &mut f64, parameters: &mut Parameters) {
    let omega = w*parameters.tau_glass;
    unsafe {
        *gp=parameters.g_glass*omega*kwws(omega,parameters.beta_glass);
        *g2p=parameters.g_glass*omega*kwwc(omega,parameters.beta_glass);

    }
}
