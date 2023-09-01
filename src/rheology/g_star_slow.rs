use crate::{DataArrays, Parameters};

fn symbint(tk: f64, td: f64, w: f64, rint1: &mut f64, rint2: &mut f64) {
    let a = tk / td;
    let b = w * w * tk * tk;
    let alpha = (1.0 + b).sqrt();
    let beta = (1.0 + alpha).sqrt();
    let gamma = (a * (alpha - 1.0)).sqrt();
    let delta = (a * (alpha + 1.0)).sqrt();
    let rt2 = 2.0f64.sqrt();

    let mut t1 = ((rt2 * gamma + a + alpha) / (a + alpha - rt2 * gamma)).ln();
    let t2 =
        (2.0 * rt2 * alpha * delta / (delta * delta + gamma * gamma - 2.0 * alpha * alpha)).atan();

    *rint1 = -gamma * beta * t1 - 2.0 * a.sqrt() * (1.0 + alpha) * t2;
    *rint2 = beta * gamma * (1.0 + alpha) * t1 / b - 2.0 * a.sqrt() * t2;
    t1 = 1.0 / (2.0 * rt2 * alpha * beta * a);
    *rint1 *= t1;
    *rint2 *= t1;
}

/**
 * \brief Calculate tube escape contrubution to relaxation moduli
 *
 *
 * \param[in] w : frequency at which result is required
 * \param[out] gp : elastic modulus G'
 * \param[out] g2p : viscous modulus G"
 * \param[out] ep : dielectric storage modulus epsilon'
 * \param[out] e2p : dielectric loss modulus epsilon"
 */
pub fn g_star_slow(
    w: f64,
    gp: &mut f64,
    g2p: &mut f64,
    ep: &mut f64,
    e2p: &mut f64,
    data_arrays: &mut DataArrays,
    parameters: &mut Parameters,
) {
    *gp = 0.0;
    *g2p = 0.0;
    *ep = 0.0;
    *e2p = 0.0;
    let mut dphi: f64;
    let mut dphi_st: f64;
    let mut tk: f64;
    let mut tm: f64;
    let mut tkm: f64;
    let mut tv: f64;

    let wsq = w * w;
    let n = data_arrays.t_ar.len() as usize;

    for k in 1..n {
        // loop over phi
        dphi = data_arrays.phi_ar[k - 1] - data_arrays.phi_ar[k];
        tk = data_arrays.t_ar[k];
        for m in 1..n {
            // loop over phi_ST

            // Dan Comment: pow is slow. If Alpha == 1, then pow() is redundant. Let's
            // make this pow() invocation conditional.
            if parameters.alpha == 1.0 {
                dphi_st = data_arrays.phi_st_ar[m - 1] - data_arrays.phi_st_ar[m];
            } else {
                dphi_st = data_arrays.phi_st_ar[m - 1].powf(parameters.alpha)
                    - data_arrays.phi_st_ar[m].powf(parameters.alpha);
            }
            tm = data_arrays.t_ar[m];
            tkm = tk * tm / (tk + tm);
            tv = tkm / (1.0 + wsq * tkm * tkm);
            *gp += tv * tkm * dphi * dphi_st;
            *g2p += tv * dphi * dphi_st;
        }
        let mut rint1 = 0.0;
        let mut rint2 = 0.0;
        symbint(tk, data_arrays.t_ar[n - 1], w, &mut rint1, &mut rint2);

        if parameters.alpha == 1.0 {
            rint1 *= 0.50 * data_arrays.phi_st_ar[n - 1];
            rint2 *= 0.50 * data_arrays.phi_st_ar[n - 1];
        } else {
            rint1 *= 0.50 * data_arrays.phi_st_ar[n - 1].powf(parameters.alpha);
            rint2 *= 0.50 * data_arrays.phi_st_ar[n - 1].powf(parameters.alpha);
        }
        *gp += rint2 * tk * tk * dphi;
        *g2p += rint1 * tk * dphi;
        tv = tk / (1.0 + wsq * tk * tk);
        *ep += tv * tk * dphi;
        *e2p += tv * dphi;
    }
    *gp *= wsq;
    *g2p *= w;
    *ep *= wsq;
    *e2p *= w;
}
