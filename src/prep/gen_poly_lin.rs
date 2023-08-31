use crate::{CLPoly, Parameters};

fn aaerfcc(x: f64) -> f64
{
    let t: f64; 
    let z: f64;
    let ans: f64;
    
  z = x.abs();
  t = 1.0 / (1.0 + 0.5 * z);
  ans = t * (-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277))))))))).exp();
  if x >= 0.0 {return ans} else {return 2.0 - ans};
}


fn log_normal_wt(mw: f64, pdi: f64, m1: f64, m2: f64 , mw_bin: &mut f64) -> f64
{
  let mu = mw.ln() - 1.50 * pdi.ln();
  let sigma = (pdi.ln()).sqrt();
  let wt_bin: f64; 
  if m1 < 0.0
  { // (0, m2)
    wt_bin = 0.50 * aaerfcc((mu + sigma * sigma - m2.ln()) / (2.0f64.sqrt() * sigma));
    let t2 = aaerfcc((mu + 2.0 * sigma * sigma - m2.ln()) / (2.0f64.sqrt() * sigma));
    *mw_bin = 0.50 * mw * t2 / wt_bin; // t1=erfc(infinity)=0
  }
  else
  {
    if m2 < 0.0
    { // (m1, infinity)
      wt_bin = 1.0 - 0.5 * aaerfcc((mu + sigma * sigma - m1.ln()) / (2.0f64.sqrt() * sigma));
      let t1 = aaerfcc((mu + 2.0 * sigma * sigma - m1.ln()) / (2.0f64.sqrt() * sigma));
      *mw_bin = 0.50 * mw * (2.0 - t1) / wt_bin;
    }
    else
    { // (m1, m2)
      let w1 = 0.5 * aaerfcc((mu + sigma * sigma - m1.ln()) / (2.0f64.sqrt() * sigma));
      let w2 = 0.5 * aaerfcc((mu + sigma * sigma - m2.ln()) / (2.0f64.sqrt() * sigma));
      wt_bin = w2 - w1;
      let t1 = aaerfcc((mu + 2.0 * sigma * sigma - m1.ln()) / (2.0f64.sqrt() * sigma));
      let t2 = aaerfcc((mu + 2.0 * sigma * sigma - m2.ln()) / (2.0f64.sqrt() * sigma));
      *mw_bin = 0.50 * mw * (t2 - t1) / wt_bin;
    }
  }
  return wt_bin;
}

/**
 * Generate polymers from a logNormal distribution characterized by molar mass mw and pdi pdi.\n
 * Special case: if either the number of discrete molar mass is one or pdi < 1, create
 * a single polymer with the molar mass supplied.
 * \param[in] n number of discrete molar mass to represent the distribution
 * \param[in] mw weight averaged molar mass
 * \param[in] pdi Polydispersity index
 * \param[in] wtcomp weight fraction of this polymer componennt
 */
pub fn gen_linlog_normal(n: i32, mw: f64, pdi: f64, wtcomp: f64, lpoly: &mut Vec<CLPoly>, parameters: &mut Parameters)
{
  if (n <= 1) || (pdi <= 1.0)
  {
    let mut ptmp = CLPoly::new();
    ptmp.mass = mw;
    ptmp.wt = wtcomp;
    ptmp.initialise(parameters.m_e);
    
    lpoly.push(ptmp);
    //npoly+=1;
    parameters.number_of_polymers += 1;
    
  }

  
  let mu = mw.ln() - 1.50 * pdi.ln();
  let sigma = (pdi.ln()).sqrt();
  let ln_m_low = mu + sigma * sigma - 2.63 * 2.0f64.sqrt() * sigma;
  let ln_m_high = mu + sigma * sigma + 2.63 * 2.0f64.sqrt() * sigma;
  let mut wt_bin: f64;
    let mut mw_bin: f64= 0.0;
    let mut m1:f64; let mut m2: f64;
  m2 = ln_m_low.exp();
  wt_bin = log_normal_wt(mw, pdi, -1.0, m2, &mut mw_bin);
  let mut ptmp = CLPoly::new();
    ptmp.mass = mw_bin; 
    ptmp.wt = wt_bin*wtcomp;
    ptmp.initialise(parameters.m_e);
  lpoly.push(ptmp);
  //npoly++;
    parameters.number_of_polymers += 1;

  let delta_ln_m = (ln_m_high - ln_m_low) / ((n as f64 - 2.0));
  for i in 0..n-2
  {
    m1 = (ln_m_low + (i as f64) * delta_ln_m).exp();
    m2 = (ln_m_low + (i as f64 + 1.0) * delta_ln_m).exp();
    wt_bin = log_normal_wt(mw, pdi, m1, m2, &mut mw_bin);
let mut ptmp = CLPoly::new();
    ptmp.mass = mw_bin; 
    ptmp.wt = wt_bin*wtcomp;
    ptmp.initialise(parameters.m_e);
  lpoly.push(ptmp);
//    npoly++;
    parameters.number_of_polymers+=1;
  }
  m1 = (ln_m_high).exp();
  wt_bin = log_normal_wt(mw, pdi, m1, -1.0, &mut mw_bin );

 
let mut ptmp = CLPoly::new();
    ptmp.mass = mw_bin; 
    ptmp.wt = wt_bin*wtcomp;
    ptmp.initialise(parameters.m_e);
  lpoly.push(ptmp);
//  npoly++;
    parameters.number_of_polymers +=1;
}
