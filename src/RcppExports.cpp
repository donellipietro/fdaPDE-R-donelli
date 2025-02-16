// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif


RcppExport SEXP _rcpp_module_boot_cpp_center();
RcppExport SEXP _rcpp_module_boot_cpp_fpca_spaceonly();
RcppExport SEXP _rcpp_module_boot_cpp_fpca_spacetimeseparable();
RcppExport SEXP _rcpp_module_boot_cpp_fpls_r_spaceonly();
RcppExport SEXP _rcpp_module_boot_cpp_fpls_a_spaceonly();
RcppExport SEXP _rcpp_module_boot_cpp_fpls_sb_spaceonly();
RcppExport SEXP _rcpp_module_boot_cpp_lagrange_basis_1d_fe1();
RcppExport SEXP _rcpp_module_boot_cpp_lagrange_basis_2d_fe1();
RcppExport SEXP _rcpp_module_boot_cpp_lagrange_basis_2d_fe2();
RcppExport SEXP _rcpp_module_boot_cpp_lagrange_basis_surface_fe1();
RcppExport SEXP _rcpp_module_boot_cpp_lagrange_basis_3d_fe1();
RcppExport SEXP _rcpp_module_boot_cpp_network_domain();
RcppExport SEXP _rcpp_module_boot_cpp_1d_domain();
RcppExport SEXP _rcpp_module_boot_cpp_2d_domain();
RcppExport SEXP _rcpp_module_boot_cpp_surface_domain();
RcppExport SEXP _rcpp_module_boot_cpp_3d_domain();
RcppExport SEXP _rcpp_module_boot_cpp_pde_1d_fe1();
RcppExport SEXP _rcpp_module_boot_cpp_pde_2d_fe1();
RcppExport SEXP _rcpp_module_boot_cpp_pde_2d_fe2();
RcppExport SEXP _rcpp_module_boot_cpp_pde_surface_fe1();
RcppExport SEXP _rcpp_module_boot_cpp_pde_3d_fe1();
RcppExport SEXP _rcpp_module_boot_cpp_regression_model();

static const R_CallMethodDef CallEntries[] = {
    {"_rcpp_module_boot_cpp_center", (DL_FUNC) &_rcpp_module_boot_cpp_center, 0},
    {"_rcpp_module_boot_cpp_fpca_spaceonly", (DL_FUNC) &_rcpp_module_boot_cpp_fpca_spaceonly, 0},
    {"_rcpp_module_boot_cpp_fpca_spacetimeseparable", (DL_FUNC) &_rcpp_module_boot_cpp_fpca_spacetimeseparable, 0},
    {"_rcpp_module_boot_cpp_fpls_r_spaceonly", (DL_FUNC) &_rcpp_module_boot_cpp_fpls_r_spaceonly, 0},
    {"_rcpp_module_boot_cpp_fpls_a_spaceonly", (DL_FUNC) &_rcpp_module_boot_cpp_fpls_a_spaceonly, 0},
    {"_rcpp_module_boot_cpp_fpls_sb_spaceonly", (DL_FUNC) &_rcpp_module_boot_cpp_fpls_sb_spaceonly, 0},
    {"_rcpp_module_boot_cpp_lagrange_basis_1d_fe1", (DL_FUNC) &_rcpp_module_boot_cpp_lagrange_basis_1d_fe1, 0},
    {"_rcpp_module_boot_cpp_lagrange_basis_2d_fe1", (DL_FUNC) &_rcpp_module_boot_cpp_lagrange_basis_2d_fe1, 0},
    {"_rcpp_module_boot_cpp_lagrange_basis_2d_fe2", (DL_FUNC) &_rcpp_module_boot_cpp_lagrange_basis_2d_fe2, 0},
    {"_rcpp_module_boot_cpp_lagrange_basis_surface_fe1", (DL_FUNC) &_rcpp_module_boot_cpp_lagrange_basis_surface_fe1, 0},
    {"_rcpp_module_boot_cpp_lagrange_basis_3d_fe1", (DL_FUNC) &_rcpp_module_boot_cpp_lagrange_basis_3d_fe1, 0},
    {"_rcpp_module_boot_cpp_network_domain", (DL_FUNC) &_rcpp_module_boot_cpp_network_domain, 0},
    {"_rcpp_module_boot_cpp_1d_domain", (DL_FUNC) &_rcpp_module_boot_cpp_1d_domain, 0},
    {"_rcpp_module_boot_cpp_2d_domain", (DL_FUNC) &_rcpp_module_boot_cpp_2d_domain, 0},
    {"_rcpp_module_boot_cpp_surface_domain", (DL_FUNC) &_rcpp_module_boot_cpp_surface_domain, 0},
    {"_rcpp_module_boot_cpp_3d_domain", (DL_FUNC) &_rcpp_module_boot_cpp_3d_domain, 0},
    {"_rcpp_module_boot_cpp_pde_1d_fe1", (DL_FUNC) &_rcpp_module_boot_cpp_pde_1d_fe1, 0},
    {"_rcpp_module_boot_cpp_pde_2d_fe1", (DL_FUNC) &_rcpp_module_boot_cpp_pde_2d_fe1, 0},
    {"_rcpp_module_boot_cpp_pde_2d_fe2", (DL_FUNC) &_rcpp_module_boot_cpp_pde_2d_fe2, 0},
    {"_rcpp_module_boot_cpp_pde_surface_fe1", (DL_FUNC) &_rcpp_module_boot_cpp_pde_surface_fe1, 0},
    {"_rcpp_module_boot_cpp_pde_3d_fe1", (DL_FUNC) &_rcpp_module_boot_cpp_pde_3d_fe1, 0},
    {"_rcpp_module_boot_cpp_regression_model", (DL_FUNC) &_rcpp_module_boot_cpp_regression_model, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_fdaPDE2(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
