pub mod prep;
pub mod relax;
pub mod rheology;

#[derive(Debug)]
pub struct CLPoly {
    pub mass: f64,
    pub wt: f64,
    pub z_chain: f64,
    pub z: f64,
    pub alive: bool,
    pub relax_free_rouse: bool,
    pub rept_set: bool,
    pub tau_d_0: f64,
    pub z_rept: f64,
    pub rept_wt: f64,
    pub p_max: i32,
    pub p_next: i32,
    pub t_f_rouse: f64,
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

    fn initialise(&mut self, m_e: f64) {
        self.z_chain = self.mass / m_e;
        self.t_f_rouse = self.z_chain * self.z_chain;
    }
}

#[derive(Debug)]
pub struct Parameters {
    // parameters from input file
    pub freq_min: f64,
    pub freq_max: f64,
    pub freq_ratio: f64,
    pub m_kuhn: f64,
    pub m_e: f64,
    pub g_0: f64,
    pub tau_e: f64,
    pub g_glass: f64,
    pub tau_glass: f64,
    pub beta_glass: f64,
    pub number_of_components: i32,
    pub p_type: Vec<i32>,
    pub wt_comp: Vec<f64>,
    pub n_poly: Vec<i32>,
    pub m_w: Vec<f64>,
    pub pdi: Vec<f64>,
    pub wt_tot: f64,

    // rc params
    pub alpha: f64,
    pub t_cr_start: f64,
    pub delta_cr: f64,
    pub b_zeta: f64,
    pub a_eq: f64,
    pub b_eq: f64,
    pub ret_pref: f64,
    pub rept_switch_factor: f64,
    pub rouse_switch_factor: f64,
    pub disentanglement_switch: f64,
    pub ret_pref_0: f64,
    pub ret_switch_exponent: f64,

    // discrete time evolution control
    pub cur_time: f64,
    pub dt_mult: f64,
    pub log_dt_mult: f64,

    // polymer parameters
    pub n_e: f64,
    pub st_max_drop: f64,
    pub rouse_wt: f64,
    pub sys_mn: f64,
    pub sys_mw: f64,
    pub sys_pdi: f64,
    pub entangled_dynamics: bool,

    // time dependent relaxation variables
    pub phi_true: f64,
    pub phi_st: f64,
    pub phi_rept: f64,
    pub phi_eq: f64,
    pub psi_rept: f64,
    pub last_reptation_time: f64,
    pub last_rept_z: f64,
    pub supertube_activated: bool,
    pub above_tau_e_first: bool,
    pub phi_st_0: f64,
    pub st_activ_time: f64,

    // count of total number of polymers
    pub number_of_polymers: i32,
    pub phi_eq_index: usize,
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
            rouse_wt: 0.0,
            sys_mn: 0.0,
            sys_mw: 0.0,
            sys_pdi: 1.0,
            entangled_dynamics: true,
            phi_true: 1.0,
            phi_st: 1.0,
            phi_rept: 1.0,
            phi_eq: 1.0,
            psi_rept: 1.0,
            last_reptation_time: 1.0,
            last_rept_z: 1.0,
            supertube_activated: false,
            above_tau_e_first: false,
            phi_st_0: 1.0,
            st_activ_time: 1.0,
            number_of_polymers: 0,
            phi_eq_index: 0,
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
#[derive(Debug)]
pub struct DataArrays {
    pub t_ar: Vec<f64>,
    pub phi_ar: Vec<f64>,
    pub phi_st_ar: Vec<f64>,
    pub t_eq_ar: Vec<f64>,
}

impl DataArrays {
    fn new() -> DataArrays {
        let mut x = DataArrays {
            t_ar: Vec::new(),
            phi_ar: Vec::new(),
            phi_st_ar: Vec::new(),
            t_eq_ar: Vec::new(),
        };
        x.t_ar.push(0.0);
        x.phi_ar.push(1.0);
        x.phi_st_ar.push(1.0);
        x.t_eq_ar.push(0.0);
        return x;
    }
}

#[derive(Debug)]
pub struct Results {
    pub freq: Vec<f64>,
    pub gp_ar: Vec<f64>,
    pub g2p_ar: Vec<f64>,
    pub ep_ar: Vec<f64>,
    pub e2p_ar: Vec<f64>,
    pub viscosity_ar: Vec<f64>,
}

impl Results {
    fn new() -> Results {
        Results {
            freq: Vec::new(),
            gp_ar: Vec::new(),
            g2p_ar: Vec::new(),
            ep_ar: Vec::new(),
            e2p_ar: Vec::new(),
            viscosity_ar: Vec::new()
        }
    }
}
extern "C" {
    pub fn kwws(win: cty::c_double, beta: cty::c_double) -> cty::c_double;
    pub fn kwwc(win: cty::c_double, beta: cty::c_double) -> cty::c_double;
}
