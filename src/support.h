using namespace Rcpp;
using namespace std;
#include <sys/time.h> // Not needed in LINUX as already loaded in armadillo, but needed in WIN
#include <chrono>

// some macro switches here
#define ___USEFIXED___
//#define ___USEGLOBAL___

#ifdef ___USEGLOBAL___
  List export2R();
  void import2cpp(const List& datain);
  #define ___DataOLG_dataOLG___
  #define ___DataOLG_dataOLG0___
  #define ___dataOLG___
  #define ___dataOLG0___
#else
  #ifdef ___USEFIXED___
    #include "auto/classfixed_init.h"
    List export2R(DataOLG* dataOLG);
    #define ___DataOLG_dataOLG___ , DataOLG* dataOLG
    #define ___DataOLG_dataOLG0___ DataOLG* dataOLG
    #define ___dataOLG___ , dataOLG
    #define ___dataOLG0___ dataOLG
  #else
    #include "auto/class_init.h"
    List export2R(DataOLG dataOLG);
    #define ___DataOLG_dataOLG___ , DataOLG& dataOLG
    #define ___DataOLG_dataOLG0___ DataOLG& dataOLG
    #define ___dataOLG___ , dataOLG
    #define ___dataOLG0___ dataOLG
  #endif
#endif

// some debugging macros
#define OK Rcout << "\n>> [" << __FILE__ << "] line " << __LINE__ << ": OK" << endl;
#define PRINT(X) Rcout << (#X) << " = " << (X) << endl;

arma::colvec onescol(int dim);
arma::rowvec onesrow(int dim);
arma::mat    onesmat(int dim1, int dim2);
arma::colvec zeroscol(int dim);
arma::rowvec zerosrow(int dim);
arma::mat    zerosmat(int dim1, int dim2);

arma::mat coh2per(const arma::mat& inmat);
arma::rowvec aggcoh2per(const arma::mat& inmat);
arma::mat per2coh(const arma::mat& inmat);

arma::ivec fixeddim();