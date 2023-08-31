use crate::CLPoly;

// This is an old conversion... haven't really checked it so if output is wrong, good to start
// here... 

#[allow(non_snake_case)]
pub fn gen_lin_log_normal(n1: &mut i32, n: i32, mw: f32, pdi: f32, wtcomp: f32, LPoly: &mut Vec<CLPoly> ) {
    let mut m_ar: Vec<f32> = vec![0.0; n as usize];
    let mut w_ar: Vec<f32> = vec![0.0; n as usize];
    let mu = mw.ln() - 1.50*pdi.ln();
    let sigma = (pdi.ln()).sqrt(); 
    let lnMlow = mu + sigma*sigma-2.63*2.0_f32.sqrt()*sigma;
    let lnMhigh =mu + sigma*sigma+2.63*2.0_f32.sqrt()*sigma;
    let M1 = lnMhigh.exp();
    let M2 = lnMlow.exp();
    
    let (wtBin,MwBin)=LogNormalWt(mw, pdi, -1., M2);
    
    m_ar[0]=MwBin; w_ar[0]=wtBin;
    let delta_lnM=(lnMhigh-lnMlow)/((n as f32-2.0));

    for i in 0..n-2 {

        let M1=(lnMlow + i as f32 *delta_lnM).exp();
        let M2=(lnMlow + (i as f32+1.0)*delta_lnM).exp();
        let (wtBin,MwBin)=LogNormalWt(mw, pdi, M1, M2);
        m_ar[i as usize +1]=MwBin; w_ar[i as usize + 1]=wtBin;  
        
    }
    let (wtBin,MwBin)=LogNormalWt(mw, pdi, M1, -1.0);
    
    m_ar[n as usize -1]=MwBin; w_ar[n as usize - 1]=wtBin;

    let mut wtot: f32=0.0;
    for i  in 0..n {
        wtot += w_ar[i as usize];
    }
    for i in 0..n{
        w_ar[i as usize] /= wtot;
    }
    

    for i in 0..n {
        LPoly.push(CLPoly {mass: m_ar[i as usize],
        wt: wtcomp*w_ar[i as usize],z: 0.,z_chain: 0.,
         alive: true, rept_set: false, relax_free_rouse: false,
         tau_d_0: 1.0e22, z_rept: 0.,rept_wt: 0.,p_max: 0,p_next: 0, t_f_rouse: 0.});
        *n1+=1;
    }


}
#[allow(non_snake_case)]
fn LogNormalWt(Mw: f32, PDI: f32, M1: f32, M2: f32) -> (f32, f32) 
{
    let mu = Mw.ln() - 1.50*PDI.ln();
    let sigma = (PDI.ln()).sqrt();
    let WtBin: f32;
    let MwBin: f32;
    if M1 < 0.0 { 
        WtBin = 0.5*aaerfcc((mu + sigma*sigma - M2.ln())/(2.0_f32.sqrt()*sigma));
        let t2=aaerfcc( (mu + 2.0*sigma*sigma - M2.ln())/(2.0_f32.sqrt()*sigma) );
        MwBin=0.50*Mw*t2/WtBin;  // t1=erfc(infinity)=0;
    }
    else if M2 < 0.0 {
        WtBin=1.0 - 0.5*aaerfcc( (mu + sigma*sigma - M1.ln())/(2.0_f32.sqrt()*sigma) );
        let t1=aaerfcc( (mu + 2.0*sigma*sigma - M1.ln())/(2.0_f32.sqrt()*sigma) );
        MwBin=0.50*Mw*(2.0 - t1)/WtBin;
    } else {
        let w1=0.5*aaerfcc( (mu + sigma*sigma - M1.ln())/(2.0_f32.sqrt()*sigma) );
        let w2=0.5*aaerfcc( (mu + sigma*sigma - M2.ln())/(2.0_f32.sqrt()*sigma) );
         WtBin=w2-w1;
        let t1=aaerfcc( (mu + 2.0*sigma*sigma - M1.ln())/(2.0_f32.sqrt()*sigma) );
        let t2=aaerfcc( (mu + 2.0*sigma*sigma - M2.ln())/(2.0_f32.sqrt()*sigma) );
        MwBin=0.50*Mw*(t2 - t1)/WtBin;
    }
    
    return (WtBin, MwBin)
}

fn aaerfcc(x: f32) -> f32 {
    let z = x.abs();
    let t=1.0/(1.0+0.5*z);
    let ans=t*(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
            t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
            t*(-0.82215223+t*0.17087277))))))))).exp();
    if x >= 0.0 {return ans} else {return 2.0-ans} 
}
