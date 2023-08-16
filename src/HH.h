#ifdef ___USEGLOBAL___
  #ifdef ___USEFIXED___
    #include "auto/globalfixed_init_ext.h"
  #else
    #include "auto/global_init_ext.h"
  #endif
#endif


void HHall(const int starttime, const bool calibinit, const double scaleA, const int nthrds ___DataOLG_dataOLG___);
void HH(const int sage0, const int z0, const int maxiter, const double stol, const double atol ___DataOLG_dataOLG___);
double HH_root(const double lambdain, const int sage0, const int z0 ___DataOLG_dataOLG___);
