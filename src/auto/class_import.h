// copy this code to auto/class_import.h and rebuild package.
// Time stamp: 2023-08-16 11:36:41

DataOLG::DataOLG(const List& datain) {

tend           = as<int>(datain["tend"]);
nag            = as<int>(datain["nag"]);
budget_bal     = as<int>(datain["budget_bal"]);
ncoh           = as<int>(datain["ncoh"]);
alpha          = as<double>(datain["alpha"]);
delta          = as<double>(datain["delta"]);
sigL           = as<double>(datain["sigL"]);
rho            = as<double>(datain["rho"]);
sigma          = as<double>(datain["sigma"]);
fag            = as<int>(datain["fag"]);
tt             = as<int>(datain["tt"]);
parlv0         = as<arma::colvec>(datain["parlv0"]);
parlv1         = as<arma::colvec>(datain["parlv1"]);
A              = as<arma::rowvec>(datain["A"]);
A0             = as<double>(datain["A0"]);
Cons           = as<arma::rowvec>(datain["Cons"]);
Cons0          = as<double>(datain["Cons0"]);
CG             = as<arma::rowvec>(datain["CG"]);
CG0            = as<double>(datain["CG0"]);
DG             = as<arma::rowvec>(datain["DG"]);
DG0            = as<double>(datain["DG0"]);
Div            = as<arma::rowvec>(datain["Div"]);
Div0           = as<double>(datain["Div0"]);
Exp            = as<arma::rowvec>(datain["Exp"]);
Exp0           = as<double>(datain["Exp0"]);
Inv            = as<arma::rowvec>(datain["Inv"]);
Inv0           = as<double>(datain["Inv0"]);
K              = as<arma::rowvec>(datain["K"]);
K0             = as<double>(datain["K0"]);
LD             = as<arma::rowvec>(datain["LD"]);
LD0            = as<double>(datain["LD0"]);
LS             = as<arma::rowvec>(datain["LS"]);
LS0            = as<double>(datain["LS0"]);
N              = as<arma::rowvec>(datain["N"]);
N0             = as<double>(datain["N0"]);
NB             = as<arma::rowvec>(datain["NB"]);
NB0            = as<double>(datain["NB0"]);
Nc             = as<arma::rowvec>(datain["Nc"]);
Nc0            = as<double>(datain["Nc0"]);
Nr             = as<arma::rowvec>(datain["Nr"]);
Nr0            = as<double>(datain["Nr0"]);
Nw             = as<arma::rowvec>(datain["Nw"]);
Nw0            = as<double>(datain["Nw0"]);
P              = as<arma::rowvec>(datain["P"]);
P0             = as<double>(datain["P0"]);
PB             = as<arma::rowvec>(datain["PB"]);
PB0            = as<double>(datain["PB0"]);
Rev            = as<arma::rowvec>(datain["Rev"]);
Rev0           = as<double>(datain["Rev0"]);
TaxF           = as<arma::rowvec>(datain["TaxF"]);
TaxF0          = as<double>(datain["TaxF0"]);
TFP            = as<arma::rowvec>(datain["TFP"]);
TFP0           = as<double>(datain["TFP0"]);
V              = as<arma::rowvec>(datain["V"]);
V0             = as<double>(datain["V0"]);
Y              = as<arma::rowvec>(datain["Y"]);
Y0             = as<double>(datain["Y0"]);
ab             = as<arma::rowvec>(datain["ab"]);
ab0            = as<double>(datain["ab0"]);
eda            = as<arma::rowvec>(datain["eda"]);
eda0           = as<double>(datain["eda0"]);
edab           = as<arma::rowvec>(datain["edab"]);
edab0          = as<double>(datain["edab0"]);
edg            = as<arma::rowvec>(datain["edg"]);
edg0           = as<double>(datain["edg0"]);
ediv           = as<arma::rowvec>(datain["ediv"]);
ediv0          = as<double>(datain["ediv0"]);
edl            = as<arma::rowvec>(datain["edl"]);
edl0           = as<double>(datain["edl0"]);
edw            = as<arma::rowvec>(datain["edw"]);
edw0           = as<double>(datain["edw0"]);
edy            = as<arma::rowvec>(datain["edy"]);
edy0           = as<double>(datain["edy0"]);
iv             = as<arma::rowvec>(datain["iv"]);
iv0            = as<double>(datain["iv0"]);
pc             = as<arma::rowvec>(datain["pc"]);
pc0            = as<double>(datain["pc0"]);
qTob           = as<arma::rowvec>(datain["qTob"]);
qTob0          = as<double>(datain["qTob0"]);
r              = as<arma::rowvec>(datain["r"]);
r0             = as<double>(datain["r0"]);
rag            = as<arma::rowvec>(datain["rag"]);
rag0           = as<double>(datain["rag0"]);
tauC           = as<arma::rowvec>(datain["tauC"]);
tauC0          = as<double>(datain["tauC0"]);
tauF           = as<arma::rowvec>(datain["tauF"]);
tauF0          = as<double>(datain["tauF0"]);
tauW           = as<arma::rowvec>(datain["tauW"]);
tauW0          = as<double>(datain["tauW0"]);
taul           = as<arma::rowvec>(datain["taul"]);
taul0          = as<double>(datain["taul0"]);
tauprof        = as<arma::rowvec>(datain["tauprof"]);
tauprof0       = as<double>(datain["tauprof0"]);
uck            = as<arma::rowvec>(datain["uck"]);
uck0           = as<double>(datain["uck0"]);
w              = as<arma::rowvec>(datain["w"]);
w0             = as<double>(datain["w0"]);
Av             = as<arma::mat>(datain["Av"]);
Az             = as<arma::mat>(datain["Az"]);
Av0            = as<arma::colvec>(datain["Av0"]);
Consv          = as<arma::mat>(datain["Consv"]);
Consz          = as<arma::mat>(datain["Consz"]);
Consv0         = as<arma::colvec>(datain["Consv0"]);
HH_nonconvv    = as<arma::mat>(datain["HH_nonconvv"]);
HH_nonconvz    = as<arma::mat>(datain["HH_nonconvz"]);
HH_nonconvv0   = as<arma::colvec>(datain["HH_nonconvv0"]);
Nv             = as<arma::mat>(datain["Nv"]);
Nz             = as<arma::mat>(datain["Nz"]);
Nv0            = as<arma::colvec>(datain["Nv0"]);
Savv           = as<arma::mat>(datain["Savv"]);
Savz           = as<arma::mat>(datain["Savz"]);
Savv0          = as<arma::colvec>(datain["Savv0"]);
abv            = as<arma::mat>(datain["abv"]);
abz            = as<arma::mat>(datain["abz"]);
abv0           = as<arma::colvec>(datain["abv0"]);
dis_totv       = as<arma::mat>(datain["dis_totv"]);
dis_totz       = as<arma::mat>(datain["dis_totz"]);
dis_totv0      = as<arma::colvec>(datain["dis_totv0"]);
cGv            = as<arma::mat>(datain["cGv"]);
cGz            = as<arma::mat>(datain["cGz"]);
cGv0           = as<arma::colvec>(datain["cGv0"]);
ellv           = as<arma::mat>(datain["ellv"]);
ellz           = as<arma::mat>(datain["ellz"]);
ellv0          = as<arma::colvec>(datain["ellv0"]);
gamv           = as<arma::mat>(datain["gamv"]);
gamz           = as<arma::mat>(datain["gamz"]);
gamv0          = as<arma::colvec>(datain["gamv0"]);
ivv            = as<arma::mat>(datain["ivv"]);
ivz            = as<arma::mat>(datain["ivz"]);
ivv0           = as<arma::colvec>(datain["ivv0"]);
lambdav        = as<arma::mat>(datain["lambdav"]);
lambdaz        = as<arma::mat>(datain["lambdaz"]);
lambdav0       = as<arma::colvec>(datain["lambdav0"]);
notretv        = as<arma::mat>(datain["notretv"]);
notretz        = as<arma::mat>(datain["notretz"]);
notretv0       = as<arma::colvec>(datain["notretv0"]);
pv             = as<arma::mat>(datain["pv"]);
pz             = as<arma::mat>(datain["pz"]);
pv0            = as<arma::colvec>(datain["pv0"]);
pcv            = as<arma::mat>(datain["pcv"]);
pcz            = as<arma::mat>(datain["pcz"]);
pcv0           = as<arma::colvec>(datain["pcv0"]);
rv             = as<arma::mat>(datain["rv"]);
rz             = as<arma::mat>(datain["rz"]);
rv0            = as<arma::colvec>(datain["rv0"]);
tauCv          = as<arma::mat>(datain["tauCv"]);
tauCz          = as<arma::mat>(datain["tauCz"]);
tauCv0         = as<arma::colvec>(datain["tauCv0"]);
tauWv          = as<arma::mat>(datain["tauWv"]);
tauWz          = as<arma::mat>(datain["tauWz"]);
tauWv0         = as<arma::colvec>(datain["tauWv0"]);
taulv          = as<arma::mat>(datain["taulv"]);
taulz          = as<arma::mat>(datain["taulz"]);
taulv0         = as<arma::colvec>(datain["taulv0"]);
thetav         = as<arma::mat>(datain["thetav"]);
thetaz         = as<arma::mat>(datain["thetaz"]);
thetav0        = as<arma::colvec>(datain["thetav0"]);
wv             = as<arma::mat>(datain["wv"]);
wz             = as<arma::mat>(datain["wz"]);
wv0            = as<arma::colvec>(datain["wv0"]);
yv             = as<arma::mat>(datain["yv"]);
yz             = as<arma::mat>(datain["yz"]);
yv0            = as<arma::colvec>(datain["yv0"]);

}
