#' @useDynLib fdaPDE2
#' @import methods Rcpp
#' @import RcppEigen
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
NULL

## load required modules


## utils
Rcpp::loadModule("cpp_network_domain", TRUE)
Rcpp::loadModule("cpp_1d_domain", TRUE)
Rcpp::loadModule("cpp_2d_domain", TRUE)
Rcpp::loadModule("cpp_surface_domain", TRUE)
Rcpp::loadModule("cpp_3d_domain", TRUE)
Rcpp::loadModule("cpp_pde_1d_fe1", TRUE)
Rcpp::loadModule("cpp_pde_2d_fe1", TRUE)
Rcpp::loadModule("cpp_pde_2d_fe2", TRUE)
Rcpp::loadModule("cpp_pde_surface_fe1", TRUE)
Rcpp::loadModule("cpp_pde_3d_fe1", TRUE)
Rcpp::loadModule("cpp_lagrange_basis_1d_fe1", TRUE)
Rcpp::loadModule("cpp_lagrange_basis_2d_fe1", TRUE)
Rcpp::loadModule("cpp_lagrange_basis_2d_fe2", TRUE)
Rcpp::loadModule("cpp_lagrange_basis_surface_fe1", TRUE)
Rcpp::loadModule("cpp_lagrange_basis_3d_fe1", TRUE)


## spatial models

## RegressionModels
Rcpp::loadModule("cpp_regression_model", TRUE)

## functional models

## centering
Rcpp::loadModule("cpp_center", TRUE)

## fPCA
Rcpp::loadModule("cpp_fpca_spaceonly", TRUE)
Rcpp::loadModule("cpp_fpca_spacetimeseparable", TRUE)

## fPLS
Rcpp::loadModule("cpp_fpls_r_spaceonly", TRUE)
Rcpp::loadModule("cpp_fpls_a_spaceonly", TRUE)
Rcpp::loadModule("cpp_fpls_sb_spaceonly", TRUE)
