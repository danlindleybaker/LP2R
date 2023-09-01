use crate::rheology::g_star_fast_rouse::g_star_fast_rouse;
use crate::rheology::g_star_glass::g_star_glass;
use crate::rheology::g_star_rouse::g_star_rouse;
use crate::rheology::g_star_slow::g_star_slow;
use crate::{CLPoly, DataArrays, Parameters, Results};
use std::io::Write;
use std::fs::File;
pub fn calc_g_star(
    parameters: &mut Parameters,
    lpoly: &mut Vec<CLPoly>,
    data_arrays: &mut DataArrays,
) -> Results {

    let mut output_file = File::create("outp.dat").unwrap();  // automatically closed when out of scope
    let mut results = Results::new();
    parameters.freq_min = parameters.freq_min * parameters.tau_e;
    parameters.freq_max = parameters.freq_max * parameters.tau_e;
    let mut freq = parameters.freq_min / parameters.freq_ratio;
    while freq < parameters.freq_max {
        freq *= parameters.freq_ratio;
        let mut gp: f64;
        let mut g2p: f64;
        let mut ep: f64 = 0.0;
        let mut e2p: f64= 0.0;

        let mut gptmp = 0.0;
        let mut g2ptmp = 0.0;
        let mut eptmp = 0.0;
        let mut e2ptmp = 0.0;

        g_star_glass(freq, &mut gptmp, &mut g2ptmp, parameters); // glassy response
        gp = gptmp;
        g2p = g2ptmp;

        g_star_rouse(
            freq,
            &mut gptmp,
            &mut g2ptmp,
            &mut eptmp,
            &mut e2ptmp,
            parameters,
            lpoly,
        ); // unentangled Rouse
        gp += gptmp;
        g2p += g2ptmp;
        ep += eptmp;
        e2p += e2ptmp;

        if parameters.entangled_dynamics {
            g_star_fast_rouse(freq, &mut gptmp, &mut g2ptmp, parameters, lpoly); // In tube Rouse modes
            gp += gptmp;
            g2p += g2ptmp;

            g_star_slow(
                freq,
                &mut gptmp,
                &mut g2ptmp,
                &mut eptmp,
                &mut e2ptmp,
                data_arrays,
                parameters,
            );
            gp += (1.0 - parameters.rouse_wt) * gptmp;
            g2p += (1.0 - parameters.rouse_wt) * g2ptmp;
            ep += (1.0 - parameters.rouse_wt) * eptmp;
            e2p += (1.0 - parameters.rouse_wt) * e2ptmp;
        }
        
        let visc = parameters.g_0 * parameters.tau_e * (gp * gp + g2p * g2p).sqrt() / freq;
        results.freq.push(freq/parameters.tau_e);
        results.gp_ar.push(gp*parameters.g_0);
        results.g2p_ar.push(g2p*parameters.g_0);
        results.ep_ar.push(ep);
        results.e2p_ar.push(e2p);
        results.viscosity_ar.push(visc);

        writeln!(&mut output_file, "{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}", freq/parameters.tau_e, gp*parameters.g_0, g2p*parameters.g_0, ep, e2p, visc).unwrap();    
    }
  results
}
