#ifdef ___USEGLOBAL___
  #ifdef ___USEFIXED___
    #include "auto/globalfixed_init_ext.h"
  #else
    #include "auto/global_init_ext.h"
  #endif
#endif

arma::colvec lifeexpect(arma::colvec gamv);
void compdemo(___DataOLG_dataOLG0___);
List solveOLG_(int starttime, int maxiter, double tol, double damping_budget, double damping_assets, double damping_ab, double damping_r, double damping_new_assets, int nthrds, List dataOLGin);

