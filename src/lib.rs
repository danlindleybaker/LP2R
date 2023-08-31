pub mod prep;


#[derive(Debug)]
pub struct CLPoly {
    pub mass: f32,
    pub wt: f32,
    pub z_chain: f32,
    pub z: f32,
    pub alive: bool,
    pub relax_free_rouse: bool,
    pub rept_set: bool,
    pub tau_d_0: f32,
    pub z_rept: f32,
    pub rept_wt: f32,
    pub p_max: i32,
    pub p_next: i32,
    pub t_f_rouse: f32,
}

impl CLPoly {
    fn new() -> CLPoly {
        CLPoly {
            mass: 0.0,
            wt: 0.0,
            z_chain: 0.0,
            z: 0.0,
            alive: true,
            relax_free_rouse: false,
            rept_set: false,
            tau_d_0: 1.0e22,
            z_rept: 0.0,
            rept_wt: 0.0,
            p_max: 0,
            p_next: 0,
            t_f_rouse: 0.0,
        }
    }
}

#[derive(Debug)]
pub struct Parameters {
    // parameters from input file
    pub freq_min: f32,
    pub freq_max: f32,
    pub freq_ratio: f32,
    pub m_kuhn: f32,
    pub m_e: f32,
    pub g_0: f32,
    pub tau_e: f32,
    pub g_glass: f32,
    pub tau_glass: f32,
    pub beta_glass: f32,
    pub number_of_components: i32,
    pub p_type: Vec<i32>,
    pub wt_comp: Vec<f32>,
    pub n_poly: Vec<i32>,
    pub m_w: Vec<f32>,
    pub pdi: Vec<f32>,
    pub wt_tot: f32,

    // rc params
    pub alpha: f32,
    pub t_cr_start: f32,
    pub delta_cr: f32,
    pub b_zeta: f32,
    pub a_eq: f32,
    pub b_eq: f32,
    pub ret_pref: f32,
    pub rept_switch_factor: f32,
    pub rouse_switch_factor: f32,
    pub disentanglement_switch: f32,
    pub ret_pref_0: f32,
    pub ret_switch_exponent: f32,

    // discrete time evolution control
    pub cur_time: f32,
    pub dt_mult: f32,
    pub log_dt_mult: f32,

    // polymer parameters
    pub n_e: f32,
    pub st_max_drop: f32,
}

impl Parameters {
    fn new() -> Parameters {
        Parameters {
            freq_min: 0.0,
            freq_max: 0.0,
            freq_ratio: 0.0,
            m_kuhn: 0.0,
            m_e: 0.0,
            g_0: 0.0,
            tau_e: 0.0,
            g_glass: 0.0,
            tau_glass: 0.0,
            beta_glass: 0.0,
            number_of_components: 0,
            p_type: Vec::new(),
            wt_comp: Vec::new(),
            n_poly: Vec::new(),
            m_w: Vec::new(),
            pdi: Vec::new(),
            wt_tot: 0.0,
            alpha: 1.0,
            t_cr_start: 1.0,
            delta_cr: 0.3,
            b_zeta: 2.0,
            a_eq: 2.0,
            b_eq: 10.0,
            ret_pref: 0.189,
            rept_switch_factor: 1.664,
            rouse_switch_factor: 1.5,
            disentanglement_switch: 1.0,
            ret_pref_0: 0.020,
            ret_switch_exponent: 0.42,
            cur_time: 1.0e-3,
            dt_mult: 1.02,
            log_dt_mult: 0.0,
            n_e: 0.0,
            st_max_drop: 0.0,
        }
    }

    fn initialise(&mut self) {
        self.log_dt_mult = self.dt_mult.ln();
        self.st_max_drop = ((-1.0 * self.log_dt_mult) / (2.0 * self.alpha)).exp();
        self.tau_glass /= self.tau_e;
        self.g_glass = (self.g_glass - 1.250 * self.g_0) / self.g_0;
        self.n_e = self.m_e / self.m_kuhn;

        for i in 0..self.number_of_components as usize {
            self.wt_tot += self.wt_comp[i]
        }
        for i in 0..self.number_of_components as usize {
            self.wt_comp[i] /= self.wt_tot;
        }
    }
}
