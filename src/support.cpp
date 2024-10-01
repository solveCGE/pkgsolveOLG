#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "support.h"

#ifdef ___USEGLOBAL___
#ifdef ___USEFIXED___
#include "auto/globalfixed_init.h"
#else
#include "auto/global_init.h"
#endif
#include "auto/global_export.h"
#include "auto/global_import.h"
#else
#ifdef ___USEFIXED___
#include "auto/classfixed_export.h"
void deleteDataOLG(DataOLG* dataOLG) { delete dataOLG; }
#else
#include "auto/class_export.h"
#endif
#include "auto/class_import.h"
#endif

arma::colvec onescol(int dim) {
  return arma::ones<arma::colvec>(dim);
}

arma::rowvec onesrow(int dim) {
  return arma::ones<arma::rowvec>(dim);
}

arma::mat onesmat(int dim1, int dim2) {
  return arma::ones<arma::mat>(dim1, dim2);
}

arma::colvec zeroscol(int dim) {
  return arma::zeros<arma::colvec>(dim);
}

arma::rowvec zerosrow(int dim) {
  return arma::zeros<arma::rowvec>(dim);
}

arma::mat zerosmat(int dim1, int dim2) {
  return arma::zeros<arma::mat>(dim1, dim2);
}

// convert variable from cohort-view to period-view
// [[Rcpp::export]]
arma::mat coh2per(const arma::mat& inmat) {
  int maxage = inmat.n_rows;
  int numcoh = inmat.n_cols;
  int numper = numcoh - (maxage - 1);

  if (numper <= 0) stop("coh2per: insufficient number of columns in input matrix");

  arma::mat outmat = zerosmat(maxage, numper);

  for (int a = 0; a < maxage; ++a) {
    outmat.row(a) = arma::rowvec(inmat.submat(a, maxage - 1 - a, a, numcoh - 1 - a));
  }

  return outmat;
}

// [[Rcpp::export]]
arma::rowvec aggcoh2per(const arma::mat& inmat) {
  return arma::rowvec(arma::sum(coh2per(inmat), 0));
}

// [[Rcpp::export]]
arma::mat per2coh(const arma::mat& inmat) {
  int maxage = inmat.n_rows;
  int numper = inmat.n_cols;
  int numcoh = numper + (maxage - 1);

  arma::colvec calibvec = inmat.col(0);

  arma::mat outmat = zerosmat(maxage, numcoh);

  for (int a = 0; a < maxage; ++a) {
    if (a < (maxage - 1)) outmat.submat(a, 0, a, maxage - 2 - a) = onesrow(maxage - a - 1) * calibvec(a);  // CHECK: one less than in R code
    outmat.submat(a, maxage - 1 - a, a, numcoh - 1 - a) = arma::rowvec(inmat.row(a));
    if (a > 0) outmat.submat(a, numcoh - a, a, numcoh - 1) = onesrow(a) * inmat(a, numper - 1);  // CHECK: one more than in R code
  }

  return outmat;
}

// [[Rcpp::export]]
arma::ivec fixeddim() {
  arma::ivec out = {___TEND___, ___NAG___};
  return out;
}