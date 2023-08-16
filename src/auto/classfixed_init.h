// copy this code to auto/classfixed_init.h and rebuild package.
// Time stamp: 2023-08-16 11:36:41

class DataOLG {

private:
~DataOLG(){} // private destructor: only allows dynamically allocated instances of DataOLG

public:
DataOLG(const List&);
friend void deleteDataOLG(DataOLG*);

int                         tend;
int                         nag;
int                         budget_bal;
int                         ncoh;
double                      alpha;
double                      delta;
double                      sigL;
double                      rho;
double                      sigma;
int                         fag;
int                         tt;
arma::colvec::fixed<100>    parlv0;
arma::colvec::fixed<100>    parlv1;
arma::rowvec::fixed<300>    A;
double                      A0;
arma::rowvec::fixed<300>    Cons;
double                      Cons0;
arma::rowvec::fixed<300>    CG;
double                      CG0;
arma::rowvec::fixed<300>    DG;
double                      DG0;
arma::rowvec::fixed<300>    Div;
double                      Div0;
arma::rowvec::fixed<300>    Exp;
double                      Exp0;
arma::rowvec::fixed<300>    Inv;
double                      Inv0;
arma::rowvec::fixed<300>    K;
double                      K0;
arma::rowvec::fixed<300>    LD;
double                      LD0;
arma::rowvec::fixed<300>    LS;
double                      LS0;
arma::rowvec::fixed<300>    N;
double                      N0;
arma::rowvec::fixed<300>    NB;
double                      NB0;
arma::rowvec::fixed<300>    Nc;
double                      Nc0;
arma::rowvec::fixed<300>    Nr;
double                      Nr0;
arma::rowvec::fixed<300>    Nw;
double                      Nw0;
arma::rowvec::fixed<300>    P;
double                      P0;
arma::rowvec::fixed<300>    PB;
double                      PB0;
arma::rowvec::fixed<300>    Rev;
double                      Rev0;
arma::rowvec::fixed<300>    TaxF;
double                      TaxF0;
arma::rowvec::fixed<300>    TFP;
double                      TFP0;
arma::rowvec::fixed<300>    V;
double                      V0;
arma::rowvec::fixed<300>    Y;
double                      Y0;
arma::rowvec::fixed<300>    ab;
double                      ab0;
arma::rowvec::fixed<300>    eda;
double                      eda0;
arma::rowvec::fixed<300>    edab;
double                      edab0;
arma::rowvec::fixed<300>    edg;
double                      edg0;
arma::rowvec::fixed<300>    ediv;
double                      ediv0;
arma::rowvec::fixed<300>    edl;
double                      edl0;
arma::rowvec::fixed<300>    edw;
double                      edw0;
arma::rowvec::fixed<300>    edy;
double                      edy0;
arma::rowvec::fixed<300>    iv;
double                      iv0;
arma::rowvec::fixed<300>    pc;
double                      pc0;
arma::rowvec::fixed<300>    qTob;
double                      qTob0;
arma::rowvec::fixed<300>    r;
double                      r0;
arma::rowvec::fixed<300>    rag;
double                      rag0;
arma::rowvec::fixed<300>    tauC;
double                      tauC0;
arma::rowvec::fixed<300>    tauF;
double                      tauF0;
arma::rowvec::fixed<300>    tauW;
double                      tauW0;
arma::rowvec::fixed<300>    taul;
double                      taul0;
arma::rowvec::fixed<300>    tauprof;
double                      tauprof0;
arma::rowvec::fixed<300>    uck;
double                      uck0;
arma::rowvec::fixed<300>    w;
double                      w0;
arma::mat::fixed<100,300>   Av;
arma::mat::fixed<100,399>   Az;
arma::colvec::fixed<100>    Av0;
arma::mat::fixed<100,300>   Consv;
arma::mat::fixed<100,399>   Consz;
arma::colvec::fixed<100>    Consv0;
arma::mat::fixed<100,300>   HH_nonconvv;
arma::mat::fixed<100,399>   HH_nonconvz;
arma::colvec::fixed<100>    HH_nonconvv0;
arma::mat::fixed<100,300>   Nv;
arma::mat::fixed<100,399>   Nz;
arma::colvec::fixed<100>    Nv0;
arma::mat::fixed<100,300>   Savv;
arma::mat::fixed<100,399>   Savz;
arma::colvec::fixed<100>    Savv0;
arma::mat::fixed<100,300>   abv;
arma::mat::fixed<100,399>   abz;
arma::colvec::fixed<100>    abv0;
arma::mat::fixed<100,300>   dis_totv;
arma::mat::fixed<100,399>   dis_totz;
arma::colvec::fixed<100>    dis_totv0;
arma::mat::fixed<100,300>   cGv;
arma::mat::fixed<100,399>   cGz;
arma::colvec::fixed<100>    cGv0;
arma::mat::fixed<100,300>   ellv;
arma::mat::fixed<100,399>   ellz;
arma::colvec::fixed<100>    ellv0;
arma::mat::fixed<100,300>   gamv;
arma::mat::fixed<100,399>   gamz;
arma::colvec::fixed<100>    gamv0;
arma::mat::fixed<100,300>   ivv;
arma::mat::fixed<100,399>   ivz;
arma::colvec::fixed<100>    ivv0;
arma::mat::fixed<100,300>   lambdav;
arma::mat::fixed<100,399>   lambdaz;
arma::colvec::fixed<100>    lambdav0;
arma::mat::fixed<100,300>   notretv;
arma::mat::fixed<100,399>   notretz;
arma::colvec::fixed<100>    notretv0;
arma::mat::fixed<100,300>   pv;
arma::mat::fixed<100,399>   pz;
arma::colvec::fixed<100>    pv0;
arma::mat::fixed<100,300>   pcv;
arma::mat::fixed<100,399>   pcz;
arma::colvec::fixed<100>    pcv0;
arma::mat::fixed<100,300>   rv;
arma::mat::fixed<100,399>   rz;
arma::colvec::fixed<100>    rv0;
arma::mat::fixed<100,300>   tauCv;
arma::mat::fixed<100,399>   tauCz;
arma::colvec::fixed<100>    tauCv0;
arma::mat::fixed<100,300>   tauWv;
arma::mat::fixed<100,399>   tauWz;
arma::colvec::fixed<100>    tauWv0;
arma::mat::fixed<100,300>   taulv;
arma::mat::fixed<100,399>   taulz;
arma::colvec::fixed<100>    taulv0;
arma::mat::fixed<100,300>   thetav;
arma::mat::fixed<100,399>   thetaz;
arma::colvec::fixed<100>    thetav0;
arma::mat::fixed<100,300>   wv;
arma::mat::fixed<100,399>   wz;
arma::colvec::fixed<100>    wv0;
arma::mat::fixed<100,300>   yv;
arma::mat::fixed<100,399>   yz;
arma::colvec::fixed<100>    yv0;

};

void deleteDataOLG(DataOLG* dataOLG);


#define ___COLVECFIXED___ arma::colvec::fixed<100>
#define ___TEND___ 300
#define ___NAG___ 100
