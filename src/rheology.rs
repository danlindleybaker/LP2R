use crate::{Parameters, CLPoly, DataArrays, Results};
use crate::rheology::calc_g_star::calc_g_star;
pub mod calc_g_star;
pub mod g_star_glass;
pub mod g_star_rouse;
pub mod g_star_fast_rouse;
pub mod g_star_slow;
pub fn lin_rheology(parameters: &mut Parameters, lpoly: &mut Vec<CLPoly>, data_arrays: &mut DataArrays) -> Results{
     calc_g_star(parameters, lpoly, data_arrays)
}
