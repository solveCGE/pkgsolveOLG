# Rcpp package for solveOLG_closed_Rcpp
An Rcpp package containing routines to solve a simple AK-OLG-model for closed economy

## About
This package is a prerequisite for running <https://github.com/solveCGE/solveOLG_closed_Rcpp>. The package provides the same routines as in <https://github.com/solveCGE/solveOLG_closed_R> but rewritten in C++ using the linear algebra library Armadillo (`Rcpparmadillo`). Solving the household problem is parallalized with OpenMP. In comparison to the pure R-implementation this results in a significant speed increase. Many initialization header files can be automatically generated using the included R function `gencppcode()`. The package experiments with different ways of memory management: using fixed dimension `arma` objects vs. non-fixed and treating all model data as global variables vs. collecting them in a struct and passing a reference to the struct between functions (including struct unpacking). The different ways can be tested by defining the macros `___USEFIXED___` and `___USEGLOBAL___` in `src/support.h`.

A model description can be found here: <https://github.com/solveCGE/solveOLG_doc>.

## How to install

- simply run `devtools::install_github("solveCGE/pkgsolveOLG")` or
- clone the repository, then open the R-project file in Rstudio and click Build -> Install.

## Author
Philip Schuster
