#include "../include/LP2R.h"
#include "../include/LP2R_global.h"
#include "../include/tclap/CmdLine.h"
#include <pybind11/pybind11.h>

namespace py = pybind11;

void reset_all_the_stuff()
{
  using namespace LP2R_NS;
  Alpha = 1.0;
  t_CR_START = 1.0;
  deltaCR = 0.30;
  B_zeta = 2.0;
  A_eq = 2.0;
  B_eq = 10.0;
  ret_pref = 0.189;
  Rept_Switch_Factor = 1.664;
  Rouse_Switch_Factor = 1.5;
  Disentanglement_Switch = 1.0;
  ret_pref_0 = 0.020;
  ret_switch_exponent = 0.42;

  // discrete time evolution control
  cur_time = 1.0e-3;
  DtMult = 1.02;

  // Output control
  FreqMin = 1.0e-3;
  FreqMax = 1.0e3;
  FreqRatio = 1.1;
  CalcDielectric = false;
  OutTermRelaxPathways = false;
  OutPhiPhiST = false;
  Output_G_of_t = false;
  OutputFormat = "Default";
  CSVdelimiter = ",";
  Add_header = true;
  has_temp = false;
  has_origin = false;
  has_label = false;
  has_chem = false;
  reptate_temp = 0.0;
  GenLogFL = false;

  // Polymer collection
  npoly = 0;
  // std::vector<C_LPoly *> LPoly;
  Rouse_wt = 0.0;
  Sys_MN = 0.0;
  Sys_MW = 0.0;
  Sys_PDI = 1.0;
  Entangled_Dynamics = true;

  // time dependent relaxation variables
  phi_true = 1.0;
  phi_ST = 1.0;
  phi_rept = 1.0;
  phi_eq = 1.0;
  Psi_rept = 1.0;
  LastReptationTime = 1.0;
  LastReptZ = 1.0;

  supertube_activated = false;
  AboveTauEFirst = false;
  phi_ST_0 = 1.0;
  ST_activ_time = 1.0;
  STmaxDrop = 1.0;
  phi_eq_indx = 0;
  t_ar.clear();
  phi_ar.clear();
  phi_ST_ar.clear();
  t_eq_ar.clear();
  LPoly.erase(LPoly.begin(), LPoly.end());
}

int LP2R_run(std::string filename, std::string out_file)
{
  using namespace LP2R_NS;
  int rtval = 0;
  // rtval = parse_arg(argc, argv);
  // if (rtval != 0) {
  //   return rtval;
  // }
  // std::cout << filename << std::endl; // uncomment to print filename
  inpFNM = filename;
  reset_all_the_stuff();
  
  rtval = ReadInput();
  if (rtval != 0)
  {
    if (GenLogFL)
    {
      f_Log << "Program ended encounting error during setup" << std::endl;
      f_Log.close();
    }
    return rtval;
  }
  int count = 0;
  if (Entangled_Dynamics)
  {
    int nalive = 0;
    nalive = time_step(0);
    while (nalive > 0)
    {
      nalive = time_step(1);
      count++;
    }
  }

  MechRelSpecFNM = out_file;
  LinRheology();

  if (OutPhiPhiST)
  {
    f_phi.close();
  }
  if (OutTermRelaxPathways)
  {
    f_trelax.close();
  }
  if (GenLogFL)
  {
    f_Log << "Program ended normally." << std::endl;
    f_Log.close();
  }
  // LPoly.clear();
  return rtval;
}

PYBIND11_MODULE(LP2R, m)
{
  m.doc() = "pybind11 example plugin";
  m.def("LP2R_run", &LP2R_run, "Function that adds");
}
