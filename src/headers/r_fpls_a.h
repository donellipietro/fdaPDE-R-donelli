// This file is part of fdaPDE, a C++ library for physics-informed
// spatial and functional data analysis.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __R_FPLS_A_H__
#define __R_FPLS_A_H__

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

// utils
#include <fdaPDE/utils/symbols.h>
#include <fdaPDE/models.h>
using fdapde::models::SpaceOnly;
using fdapde::models::SpaceTimeSeparable;

// pde
#include "r_pde.h"

// fpls
#include <fdaPDE/models/functional/fpls_base.h>
#include <fdaPDE/models/functional/fpls_a.h>
using fdapde::models::FPLS_BASE;
using fdapde::models::FPLS_A;

// rsvd
#include <fdaPDE/models/functional/regularized_svd.h>
#include "r_rsvd.h"
using fdapde::models::RSVDType;

// generic fPLS-R model wrapper signature
template<typename RegularizationType> class R_FPLS_A {
  public:
    using ModelType = std::decay_t<FPLS_A<RegularizationType>>;
    using ModelBaseType = std::decay_t<FPLS_BASE<RegularizationType, ModelType>>;
    using SolverType = RSVDType<ModelBaseType>;
    // constructor
    R_FPLS_A() : calibration_strategy_(Calibration::off) {}
    // getters
    const SpMatrix<double>&  Psi() const { return model_.Psi(); }
    const SpMatrix<double>&  R0() const { return model_.R0(); }
    DMatrix<double> fitted(std::size_t h){ return model_.fitted(h); }
    DMatrix<double> reconstructed(std::size_t h){ return model_.reconstructed(h); }
    DMatrix<double> Y_space_directions(){ return model_.Y_space_directions(); }
    DMatrix<double> Y_loadings(){ return model_.Y_loadings(); }
    DMatrix<double> X_space_directions(){ return model_.X_space_directions(); }
    DMatrix<double> X_loadings(){ return model_.X_loadings(); }
    DMatrix<double> X_latent_scores(){ return model_.X_latent_scores(); }
    DMatrix<double> Y_latent_scores(){ return model_.Y_latent_scores(); }
    // setters
    void set_data(const Rcpp::List & data) {
      BlockFrame<double, int> df;
      df.template insert<double>(OBSERVATIONS_BLK, Rcpp::as<DMatrix<double>>(data["Y"]));
      df.template insert<double>(DESIGN_MATRIX_BLK, Rcpp::as<DMatrix<double>>(data["X"]));
      model_.set_data(df);
    }
    void set_spatial_locations(const DMatrix<double>& locs) { model_.set_spatial_locations(locs); }
    void set_ncomp(std::size_t n_comp) { model_.set_ncomp(n_comp); }
    void set_lambda(Rcpp::List R_lambda){
      // parse the R input
      Rcpp::NumericVector lambda_D_grid = R_lambda["space"];
      // Rcpp::NumericVector lambda_T_grid = R_lambda["time"];
      // move Rcpp::NumericVector into something understandable by cpp layer
      lambda_grid_.resize(lambda_D_grid.size());
      for(std::size_t i = 0; i < lambda_grid_.size(); ++i){
        lambda_grid_[i].resize(1); // 2
        lambda_grid_[i][0] = lambda_D_grid[i];
        // lambda_grid_[i][1] = lambda_T_grid[i];
      }
    }
    void set_solver(int policy, Rcpp::List rsvd_params,
                    int calibration_strategy, Rcpp::List calibrator_params, Rcpp::List lambda) {
      rsvd_ = R_RSVD<ModelBaseType>(RSVDSolutionPolicy(policy), rsvd_params, Calibration(calibration_strategy), calibrator_params);
      calibration_strategy_ = rsvd_.calibration();
      set_lambda(lambda);
    }
    // utilities
    void init() { 
      load_solver();
      model_.init();
    }
    void solve() { model_.solve(); }
  protected:
    ModelType model_;
    R_RSVD<ModelBaseType> rsvd_;
    Calibration calibration_strategy_;
    std::vector<DVector<double>> lambda_grid_;
    // utilities
    void load_solver() {
      if(calibration_strategy_ == Calibration::off){
        model_.set_rsvd(rsvd_());
        model_.set_lambda_D(lambda_grid_.front()[0]);
      } else {
        model_.set_rsvd(rsvd_(lambda_grid_));
      }
    }
};

// implementation of the fPLS model wrapper for SpaceOnly regulatization
class R_FPLS_A_SpaceOnly : public R_FPLS_A<SpaceOnly> {
  public:
    R_FPLS_A_SpaceOnly(Rcpp::Environment pde,
                      int sampling_type,
                      Rcpp::List fPLS_params) {
        // recover pointer to penalty
        SEXP pdeptr = pde[".pointer"];
        PDEWrapper* ptr = reinterpret_cast<PDEWrapper*>(R_ExternalPtrAddr(pdeptr));
        // set model instance
        model_ = ModelType(ptr->get_pde(), Sampling(sampling_type));
        // model configuration
        model_.set_ncomp(fPLS_params["n_comp"]);
    }
};

#endif // __R_FPLS_A_H__