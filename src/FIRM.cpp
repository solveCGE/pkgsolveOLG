#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "support.h"
#include "FIRM.h"

// production function and marginal products
arma::rowvec fY(const arma::rowvec& K, const arma::rowvec& L, const arma::rowvec& TFP, double alpha) {return TFP%pow(K,alpha)%pow(L,1-alpha);}
arma::rowvec MPK(const arma::rowvec& K, const arma::rowvec& L, const arma::rowvec& TFP, double alpha) {return alpha*TFP%pow(K/L,alpha-1);}
arma::rowvec MPKinv(const arma::rowvec& MPKin, const arma::rowvec& L, const arma::rowvec& TFP, double alpha) {return L%pow(MPKin/(TFP*alpha),1/(alpha-1));}
arma::rowvec MPL(const arma::rowvec& K, const arma::rowvec& L, const arma::rowvec& TFP, double alpha) {return (1-alpha)*TFP%pow(K/L,alpha);}

// invert the firm problem: finds r for given V (and LS, TFP and taxes)
arma::rowvec rdemand(const arma::rowvec& assetsupply ___DataOLG_dataOLG___, int maxiter, double tol, bool verbose) {
  
  #include "unpack.h"
  
  arma::rowvec K2      = K;
  arma::rowvec uck2    = uck;
  arma::rowvec r2      = r;
  arma::rowvec qTob2   = qTob;
  
  double error = arma::datum::inf;
  int iter     = 0;
  
  double error_old;
  arma::rowvec qTob_old;
  
  while(true) {
    
    iter++;
    error_old      = error;
    qTob_old       = qTob2;
    
    K2             = assetsupply/qTob2;
    //K2(0)        = K(0); //predetermined
    uck2           = MPK(K2,LD,TFP,alpha);
    qTob2          = (1-tauprof)%uck2 + tauprof*delta + (1-delta);
    
    error = arma::sum(abs(qTob2-qTob_old));
    
    if (verbose) Rcout << "Iteration:\t" << iter << "\t\tError:\t" << error << "\n";
    
    if (iter > maxiter)    {if (verbose) Rcout << "No convergence!!\n"; break;}
    if (error < tol)       {if (verbose) Rcout << "Convergence!!\n"; break;}
    if (error > error_old) {if (verbose) Rcout << "Increasing error: stop at previous step!\n"; qTob2 = qTob_old; break;}
    
  }
  
  r2.subvec(0,tend-2) = qTob2.subvec(1,tend-1)-1;
  r2(tend-1) = r2(tend-2);
  
  return r2;
}

