#include <pybind11/pybind11.h>
#include "../include/LP2R.h"
#include "../include/LP2R_global.h"   
#include "../include/tclap/CmdLine.h"

namespace py = pybind11;

int LP2R() {
  using namespace LP2R_NS;
  int rtval = 0;
  // rtval = parse_arg(argc, argv);
  // if (rtval != 0) {
  //   return rtval;
  // }

  rtval = ReadInput();

  if (rtval != 0) {
    if (GenLogFL) {
      f_Log << "Program ended encounting error during setup" << std::endl;
      f_Log.close();
    }
    return rtval;
  }

  if (Entangled_Dynamics) {
    int nalive = 0;
    nalive = time_step(0);
    while (nalive > 0) {
      nalive = time_step(1);
    }
  }

  LinRheology();

  if (OutPhiPhiST) {
    f_phi.close();
  }
  if (OutTermRelaxPathways) {
    f_trelax.close();
  }
  if (GenLogFL) {
    f_Log << "Program ended normally." << std::endl;
    f_Log.close();
  }

  return rtval;
}


PYBIND11_MODULE(LP2R, m) {
  m.doc() = "pybind11 example plugin";
  m.def("LP2R", &LP2R, "Function that adds");
}
