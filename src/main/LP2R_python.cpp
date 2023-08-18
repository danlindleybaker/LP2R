/** \file
 * \brief \c main for LP2R
 * \param[in] argc Number of command line arguments
 * \param[in] argv argument list
 * \return nonzero integer is run fails
 */
// need to rewrite this such that it can be used as an entry
// point for python module rather than CLI app
//
//

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

//int main(int argc, char *argv[]) { LP2R(); }
