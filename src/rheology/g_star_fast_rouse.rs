use crate::{CLPoly, Parameters};

pub fn g_star_fast_rouse(
    w: f64,
    gpf: &mut f64,
    g2pf: &mut f64,
    parameters: &mut Parameters,
    lpoly: &mut Vec<CLPoly>,
) {
    *gpf = 0.0;
    *g2pf = 0.0;

    let mut gp: f64;
    let mut g2p: f64;
    let mut tv: f64;
    let mut tv2: f64;
    let mut tv4: f64;
    let w2 = w * w;
    let mut zi: i64;

    for i1 in 0..parameters.number_of_polymers as usize {
        if !lpoly[i1].relax_free_rouse {
            let zz = lpoly[i1].z_chain;
            zi = zz.ceil() as i64;
            gp = 0.0;
            g2p = 0.0;

            for i in 1..zi as usize {
                // longitudinal contribution
                tv = i as f64 / zz;
                tv2 = tv * tv;
                tv4 = tv2 * tv2;
                gp += 1.0 / (w2 + tv4); // multiply by w^2/(4*Z) later
                g2p += tv2 / (w2 + tv4);
            }

            let max_term = (parameters.n_e * zz).ceil() as i64;
            for i in zi..max_term {
                // internal Rouse modes
                tv = i as f64 / zz;
                tv2 = tv * tv;
                tv4 = tv2 * tv2;
                tv = 1.0 / (w2 + 4.0 * tv4);
                gp += 5.0 * tv;
                g2p += 10.0 * tv2 * tv;
            }
            *gpf += w2 * gp * lpoly[i1].wt / (4.0 * zz);
            *g2pf += w * g2p * lpoly[i1].wt / (4.0 * zz);
        }
    }
}
