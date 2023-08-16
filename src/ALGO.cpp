#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

#include <omp.h>
// [[Rcpp::plugins(openmp)]]

#include "support.h"
#include "ALGO.h"
#include "FIRM.h"
#include "HH.h"

// computes demographic transition (and updates intervivo transfers accordingly)
void compdemo(___DataOLG_dataOLG0___) {

  #include "unpack.h"

  // compute demography transition
  for (int t = 1; t < tend; ++t) {
    Nv(0,t)  = NB(t);
    for (int i = 1; i < nag; ++i) {
      Nv(i,t) = Nv(i-1,t-1)*gamv(i-1,t-1);
    }
  }

  Nz = per2coh(Nv);
  N  = aggcoh2per(Nz);
  Nc = arma::sum(arma::mat(Nv.submat(0,0,fag-2,tend-1)),0);

  // Compute neutral intervivo-transfers by rescaling received transfers
  arma::Row<int> ivverror(tend); ivverror.zeros();
  
  #pragma omp parallel for schedule(static)
  for (int t = 0; t < tend; ++t) {
    double ivgiven          = -arma::accu(Nv.col(t)%ivv.col(t)%(ivv.col(t)<0));
    double ivreceived       = arma::accu(Nv.col(t)%ivv.col(t)%(ivv.col(t)>0));
    arma::colvec ivvtemp    = ivv.col(t);
    ivvtemp.elem(find(ivvtemp>0)) = ivvtemp.elem(find(ivvtemp>0))*(ivgiven/ivreceived);
    ivv.col(t)              = ivvtemp;
    
    if (abs(arma::accu(ivv.col(t)%Nv.col(t)))>1e-10) ivverror(t) = 1;
  }
  
  if (arma::sum(ivverror) > 0) stop("ERROR IN UPDATEDEMO: Unbalanced intervivo transfers!");

  ivz = per2coh(ivv);

}

// computes the further life expectancy
// [[Rcpp::export]]
arma::colvec lifeexpect(const arma::colvec& gamv) {
  
  const int nag = gamv.n_elem;
  
  arma::colvec lifeexpectageN = zeroscol(nag);
  for (int a = nag-2; a >= 0; a--) {
    lifeexpectageN(a) = (lifeexpectageN(a+1)+1)*gamv(a)+gamv(a+1)*(1-gamv(a));
  }
  
  return lifeexpectageN;
}


// [[Rcpp::export]]
List solveOLG_(int starttime, int maxiter, double tol, double damping_budget, double damping_assets, double damping_ab, double damping_r, double damping_new_assets, int nthrds, List dataOLGin) {
  
  // converts Rcpp-List dataOLGin to individual variables
  #include "import.h"
  #include "unpack.h"
  
  double scaleA            = 1.0;           // initialize
  arma::rowvec scaleab     = onesrow(tend); // initialize
  
  std::ios_base::fmtflags myRcoutflags( Rcout.flags() );
  std::chrono::high_resolution_clock::time_point tstart_algo, tend_algo, tstart_loop, tend_loop;
  long long time_algo, time_loop;
  
  //===== demography ======//
  compdemo(___dataOLG0___); // recomputes demographic transition
  
  tstart_algo = std::chrono::high_resolution_clock::now();
  
  for (int iter = 1; iter <= maxiter; ++iter) {
    
    tstart_loop = std::chrono::high_resolution_clock::now();
    
    //===== solve the firm problem for given labor demand ======//
    uck.subvec(starttime,tend-1)   = (r.subvec(starttime-1,tend-2)+delta*(1-tauprof.subvec(starttime,tend-1)))/(1-tauprof.subvec(starttime,tend-1));
    K.subvec(starttime,tend-1)     = MPKinv(uck.subvec(starttime,tend-1),LD.subvec(starttime,tend-1),TFP.subvec(starttime,tend-1),alpha);
    Inv.subvec(starttime-1,tend-2) = K.subvec(starttime,tend-1) - (1-delta)*K.subvec(starttime-1,tend-2);
    
    Inv(tend-1)     = delta*K(tend-1);
    qTob            = (1-tauprof)%MPK(K,LD,TFP,alpha) + tauprof*delta + (1-delta);
    
    Y               = fY(K,LD,TFP,alpha);
    w               = MPL(K,LD,TFP,alpha)/(1+tauF);
    wv              = arma::kron(w,onescol(nag));
    wz              = per2coh(wv);
    //wz              = per2coh(w,nag); // <- IS THIS IMPLEMENTED???!!! rather convert first to wv
    V               = qTob%K;
    TaxF            = tauprof%(Y-(1+tauF)%w%LD-delta*K);
    Div             = Y-(1+tauF)%w%LD-Inv-TaxF;
    
    //===== solve the households' problem for given prices and tax rates ======//
    HHall(starttime, (iter == 1), scaleA, nthrds ___dataOLG___);
    
    //===== aggregation ======//
    Cons      = aggcoh2per(Consz%Nz);
    LS        = aggcoh2per(notretz%ellz%thetaz%Nz);
    A         = aggcoh2per(Az%Nz);
    ab        = aggcoh2per(abz%Nz);
    iv        = aggcoh2per(ivz%Nz); // should be 0 by construction
    Nw        = aggcoh2per(notretz%Nz);
    Nr        = aggcoh2per((1-notretz)%Nz);
    
    // government budget
    P         = aggcoh2per((1-notretz)%pz%Nz);
    tauW      = aggcoh2per(tauWz%notretz%ellz%thetaz%Nz)/LS;
    arma::rowvec TaxP = aggcoh2per((1-notretz)%tauWz%pz%Nz);
    arma::rowvec Taxl = aggcoh2per(taulz%Nz);
    Rev       = TaxF+(tauF%LD+tauW%LS)%w+Taxl+tauC%Cons+TaxP;
    CG        = arma::sum(cGv%Nv,0);
    Exp       = CG+P;
    
    // follow given debt-path
    PB.subvec(starttime-1,tend-2)  = DG.subvec(starttime-1,tend-2)-DG.subvec(starttime,tend-1)/(1+r.subvec(starttime-1,tend-2));
    PB(tend-1)                     = r(tend-1)*DG(tend-1)/(1+r(tend-1));
    
    //===== excess demands ======//
    edy       = Inv+Cons+CG-Y;
    edg       = Rev-Exp-PB;
    edl       = LD-LS;
    eda       = DG+V-A;
    ediv      = -iv;
    edab      = aggcoh2per((1-gamz)%Savz%Nz)-ab;
    
    arma::rowvec eda_lag     = eda;
    eda_lag.subvec(0,tend-2) = eda.subvec(1,tend-1);
    edw       = 1*edy + w%edl + ediv + edab + edg + eda - eda_lag/(1+r); // Walras' Law
    
    // check Walras' Law: this always has to hold (even out of equilibrium)! If not there is something wrong with accounting in the model
    if (max(abs(edw.subvec(starttime-1,tend-2)))> 1e-10) stop("Error: Walras Law does not hold!");
    
    tend_loop = std::chrono::high_resolution_clock::now();
    time_loop = chrono::duration_cast<chrono::nanoseconds>(tend_loop-tstart_loop).count();
    
    //===== checking error and breaking loop ======//
    double err             = arma::accu(abs(edy.subvec(starttime-1,tend-1)))+sum(abs(edg.subvec(starttime-1,tend-1)))+sum(abs(edl.subvec(starttime-1,tend-1)))+sum(abs(eda.subvec(starttime-1,tend-1)))+sum(abs(ediv.subvec(starttime-1,tend-1)))+sum(abs(edab.subvec(starttime-1,tend-1)));
    double err2            = log(err/tol);

    Rcout << "Iteration:  "       << std::fixed << std::setw( 3 ) << std::setprecision( 0 ) << std::setfill( ' ' ) << iter;
    Rcout << "  scaleA: "         << std::fixed << std::setw( 7 ) << std::setprecision( 6 ) << std::setfill( ' ' ) << scaleA;
    Rcout << "  scaleab: "        << std::fixed << std::setw( 7 ) << std::setprecision( 6 ) << std::setfill( ' ' ) << arma::mean(scaleab); 
    Rcout << "  non-conv.HH: "    << std::fixed << std::setw( 4 ) << std::setprecision( 0 ) << std::setfill( ' ' ) << arma::accu(HH_nonconvz); 
    Rcout << "  Time: "           << std::fixed << std::setw( 7 ) << std::setprecision( 2 ) << std::setfill( ' ' ) << double(time_loop)/1000000 << " ms"; 
    Rcout << "  log of err/tol: " << std::fixed << std::setw( 11 ) << std::setprecision( 8 ) << std::setfill( ' ' ) << err2 << "\n";

    if (err2 < 0.0) {    
      Rcout << "                                                                                                            Convergence!\n\n";
      break;
    }
    if (iter == maxiter) {    
      Rcout << "                                                                                                            No Convergence!\n\n";
      break;
    }
    
    //======= updating for next iteration =======//
    // budget rules
    arma::rowvec budget_surplus  = edg*damping_budget;
    
    if (budget_bal == 1) {
      tauWv     = tauWv - arma::kron(budget_surplus/(w%LS),onescol(nag));
      tauWz     = per2coh(tauWv);
    }
    if (budget_bal == 2) {
      tauF       = tauF - budget_surplus/(w%LD); 
    }
    if (budget_bal == 3) {
      tauC       = tauC - budget_surplus/Cons; 
      tauCv      = arma::kron(tauC,onescol(nag));
      tauCz      = per2coh(tauCv);
    }
    if (budget_bal == 4) {
      taul            = taul - budget_surplus/(N-Nc);
      taulv.submat(fag-1,0,nag-1,tend-1) = arma::kron(taul,onescol(nag-fag+1));
      taulz           = per2coh(taulv);
    } 
    if (budget_bal == 5) {
      tauprof    = tauprof - budget_surplus/(Y-(1+tauF)%w%LD-delta*K); 
    }
    if (budget_bal == 6) {
      cGv        = cGv + kron(budget_surplus/N,onescol(nag));
      CG         = arma::sum(cGv%Nv,0);
    }
    
    // price updating
    arma::rowvec newassets       = damping_new_assets*(A-DG) + (1-damping_new_assets)*V;
    arma::rowvec r_new           = rdemand(newassets ___dataOLG___);
    r               = damping_r*r_new + (1-damping_r)*r;
    rv              = arma::kron(r,onescol(nag));
    //rz              = per2coh(r,nag);
    rz              = per2coh(rv);
    
    scaleab         = 1+(aggcoh2per((1-gamz)%Savz%Nz)/ab-1)*damping_ab;
    abv             = abv%arma::kron(scaleab,onescol(nag));
    abz             = per2coh(abv);
    LD              = LS;
    scaleA          = 1+((DG(starttime-1)+V(starttime-1))/A(starttime-1)-1)*damping_assets;

  }
  
  tend_algo = std::chrono::high_resolution_clock::now();
  time_algo = chrono::duration_cast<chrono::nanoseconds>(tend_algo-tstart_algo).count();
  
  Rcout << "Computation time:     " << std::fixed << std::setw( 8 ) << std::setprecision( 4 ) << std::setfill( ' ' ) << double(time_algo)/1000000000 << " sec" << endl;
  
  // convert cohort-view variables back to period-view variables
  // (those where only cohort-view variables were altered in solveOLG)
  Av       = coh2per(Az);
  Consv    = coh2per(Consz);
  lambdav  = coh2per(lambdaz);
  Savv     = coh2per(Savz);
  dis_totv = coh2per(dis_totz);
  ellv     = coh2per(ellz);
  pcv      = coh2per(pcz);
  yv       = coh2per(yz);
 
  // collects individual variables again in Rcpp-List dataOLGout
  #include "export.h"
  
  return dataOLGout;
}
