#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)))

#include <omp.h>
// [[Rcpp::plugins(openmp)))

#include "support.h"
#include "HH.h"

void HHall(const int starttime, const bool calibinit, const double scaleA, const int nthrds ___DataOLG_dataOLG___) {

#include "unpack.h"
  
#pragma omp parallel for num_threads(nthrds) schedule(dynamic)
  for (int z0 = starttime; z0 <= ncoh; ++z0) {
    if (z0 <= nag-fag+starttime-1) {
      if (calibinit == true) {
        Az.col(z0-1) = Av0;
      }
      Az(nag-1-(z0-starttime),z0-1) *= scaleA;
      HH(nag-(z0-starttime), z0, 30, 1e-10, 0.1 ___dataOLG___);
    } else {
      HH(fag, z0, 30, 1e-10, 0.1 ___dataOLG___);
    }
  }
}

void HH(const int sage0, const int z0, const int maxiter, const double stol, const double atol ___DataOLG_dataOLG___) {
  
  #include "unpack.h"
  
  double err     = arma::datum::inf;
  int iter       = 0;
  int trys       = 0;
  const double stepsize = 1e-6; // for numerical gradient
  
  const arma::vec lambdatrys = {1.0,0.5,1.5,0.25,1.25,0.1,1.0};
  const int maxtrys          = lambdatrys.n_elem;
  bool while_continue        = true;
  
  double lambdazsave, lambdaz0, lambdaz1, lambdaz2, f0, f1, f2;
  bool breakwhile;
  int iterpertry;
  
  while (while_continue) {
    
    while_continue = false;
    lambdazsave    = lambdaz(sage0-1,z0-1);
    
    while (((err > stol)||(abs(Savz(nag-1,z0-1)) > atol)) && (trys < maxtrys)) {
      
      trys++;
      iterpertry = 0;
      lambdaz1 = lambdazsave*lambdatrys(trys-1);
      
      breakwhile = false;
      while ((err > stol) && (iterpertry < maxiter) && (breakwhile == FALSE)) {
        if (iterpertry == 0) { // Newton step for first iteration
        f2 = HH_root(lambdaz1+stepsize,sage0,z0 ___dataOLG___); iter++;
        
        if (!std::isfinite(f2)) {breakwhile = true; break;}
        f1 = HH_root(lambdaz1,sage0,z0 ___dataOLG___); iter++;
        
        if (!std::isfinite(f1)) {breakwhile = true; break;}
        lambdaz2 = lambdaz1 - f1*stepsize/(f2-f1);
        if (!std::isfinite(lambdaz2)||(lambdaz2<0)) {breakwhile = true; break;}
        } else { // Secant method
          f1 = HH_root(lambdaz1,sage0,z0 ___dataOLG___); iter++;
          
          if (!std::isfinite(f1)) {breakwhile = true; break;}
          lambdaz2 = lambdaz1 - f1*(lambdaz1-lambdaz0)/(f1-f0);
          if (!std::isfinite(lambdaz2)||(lambdaz2<0)) {breakwhile = true; break;}
        }
        err = abs(lambdaz2-lambdaz1);
        lambdaz0 = lambdaz1;
        lambdaz1 = lambdaz2;
        f0       = f1;
        iterpertry++;
      }
    }
  }
  
  if (abs(Savz(nag-1,z0-1)) > atol) {
    HH_nonconvz(z0-1) = 1; // counter
  }
  
}

double HH_root(const double lambdain, const int sage0, const int z0 ___DataOLG_dataOLG___) {
  
  const int sage = sage0 - 1;
  const int z    = z0 - 1;
  
  #include "unpack.h"
  
  // EULER EQUATION: solve forward in age
  lambdaz(sage,z) = lambdain;
  if (sage < nag-1) {
    for (int a = sage; a < nag-1; ++a) {
      lambdaz(a+1,z) = lambdaz(a,z)/((1/(1+rho))*gamz(a,z)*(1+rz(a,z)));
    }
  }
  
  // CONSUMPTION
  pcz.submat(sage,z,nag-1,z)      = 1+tauCz.submat(sage,z,nag-1,z);
  Consz.submat(sage,z,nag-1,z)   = pow(pcz.submat(sage,z,nag-1,z)%lambdaz.submat(sage,z,nag-1,z),-sigma);
  
  // HOURS SUPPLY
  ellz.submat(sage,z,nag-1,z)     = pow((wz.submat(sage,z,nag-1,z)%(1-tauWz.submat(sage,z,nag-1,z))%thetaz.submat(sage,z,nag-1,z)/pcz.submat(sage,z,nag-1,z)%pow(Consz.submat(sage,z,nag-1,z),-1/sigma))/parlv0.subvec(sage,nag-1),sigL);
  dis_totz.submat(sage,z,nag-1,z) = (sigL/(1+sigL))*parlv0.subvec(sage,nag-1)%pow(ellz.submat(sage,z,nag-1,z),(1+sigL)/sigL)-parlv1.subvec(sage,nag-1);
  
  // CONSUMPTION AND SAVINGS
  yz.submat(sage,z,nag-1,z)       = notretz.submat(sage,z,nag-1,z)%(wz.submat(sage,z,nag-1,z)%(1-tauWz.submat(sage,z,nag-1,z))%ellz.submat(sage,z,nag-1,z)%thetaz.submat(sage,z,nag-1,z))+(1-notretz.submat(sage,z,nag-1,z))%(1-tauWz.submat(sage,z,nag-1,z))%pz.submat(sage,z,nag-1,z)-taulz.submat(sage,z,nag-1,z);
  
  // ASSETS: solve forward in age
  Az(0,z) = 0;
  
  if (sage < nag) {
    for (int a = sage; a < nag-1; ++a) {
      Az(a+1,z)   = (1+rz(a,z))*(Az(a,z)+yz(a,z)+ivz(a,z)+abz(a,z)-pcz(a,z)*Consz(a,z)); // if sage > 1 take previous age entry in Az as starting value! (i.e. has to be given globally not passed in function)
    }
  }
  
  Savz.submat(sage,z,nag-1,z)  = Az.submat(sage,z,nag-1,z)+yz.submat(sage,z,nag-1,z)+ivz.submat(sage,z,nag-1,z)+abz.submat(sage,z,nag-1,z)-pcz.submat(sage,z,nag-1,z)%Consz.submat(sage,z,nag-1,z);
  
  return Savz(nag-1,z);
}