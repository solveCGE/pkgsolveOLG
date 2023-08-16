#ifdef ___USEGLOBAL___
  #ifdef ___USEFIXED___
    #include "auto/globalfixed_init_ext.h"
  #else
    #include "auto/global_init_ext.h"
  #endif
#endif

arma::rowvec fY(const arma::rowvec& K, const arma::rowvec& L, const arma::rowvec& TFP, double alpha);
arma::rowvec MPK(const arma::rowvec& K, const arma::rowvec& L, const arma::rowvec& TFP, double alpha);
arma::rowvec MPKinv(const arma::rowvec& MPKin, const arma::rowvec& L, const arma::rowvec& TFP, double alpha);
arma::rowvec MPL(const arma::rowvec& K, const arma::rowvec& L, const arma::rowvec& TFP, double alpha);

arma::rowvec rdemand(const arma::rowvec& assetsupply ___DataOLG_dataOLG___, int maxiter = 20, double tol = 1e-6, bool verbose = false);
