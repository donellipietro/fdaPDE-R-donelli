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

#ifndef __R_REGRESSION_MODEL_H__
#define __R_REGRESSION_MODEL_H__

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

// pde utilities
#include "r_pde.h"

// regression models
#include <fdaPDE/models/regression/regression_base.h>
#include <fdaPDE/models/regression/srpde.h>
#include <fdaPDE/models/regression/strpde.h>
#include <fdaPDE/models/regression/gsrpde.h>
using fdapde::models::RegressionView;
using fdapde::models::SRPDE;
using fdapde::models::STRPDE;
using fdapde::models::GSRPDE;
enum RegressionModelEnum{ srpde };

// calibrators
#include "r_calibrator.h"

class R_REGRESSION_MODEL {
  private:

    // model to wrap
    RegressionModelEnum regression_model_;
    fdapde::models::SRPDE model_srpde_;
    RegressionView<void> model_view_;

    // calibrator
    R_CALIBRATOR calibrator_;

    // data
    BlockFrame<double, int> data_;

  public:

    // constructor
    R_REGRESSION_MODEL(int regression_model,
                       const Rcpp::Environment & pde, int sampling_type, const Rcpp::List & smoother_params)
                       : regression_model_(RegressionModelEnum(regression_model)) {
      switch (regression_model_) {
        case RegressionModelEnum::srpde : {
          // recover pointer to penalty
          SEXP pdeptr = pde[".pointer"];
          PDEWrapper* ptr = reinterpret_cast<PDEWrapper*>(R_ExternalPtrAddr(pdeptr));
          // set model instance
          model_srpde_ = SRPDE(ptr->get_pde(), Sampling(sampling_type));
          // smoother customization
          // ... SRPDE does not have any smoother parameter yet
          // create a model view
          model_view_ = model_srpde_;
          break;
        }
        // other regression models ...
      }
    }

    // getters
    RegressionView<void> get_model_view() const { return model_view_; }
    Calibrator<RegressionView<void>> get_configured_calibrator() { return calibrator_.get_configured_calibrator(); }
    const DVector<double>& f() const { return model_view_.f(); }
    const DVector<double>& beta() const { return model_view_.beta(); }
    int get_calibration_strategy() { return calibrator_.get_calibration_strategy();}

    // setters
    void set_observations(const DMatrix<double>& y) { data_.template insert<double>(OBSERVATIONS_BLK, y); }
    void set_covariates(const DMatrix<double>& X) { data_.template insert<double>(DESIGN_MATRIX_BLK, X); }
    void set_spatial_locations(const DMatrix<double>& locs) { model_view_.set_spatial_locations(locs); }
    void set_calibrator(int calibration_strategy, const Rcpp::List & calibrator_params, const Rcpp::List & R_lambda) {
      calibrator_ = R_CALIBRATOR{Calibration(calibration_strategy), calibrator_params};
      calibrator_.configure_calibrator(R_lambda);
    }
    void set_lambda(Rcpp::List lambda) {
      model_view_.set_lambda_D(lambda["space"]);
      // if(lambda["time"] != R_NilValue)
      //  model_view_.set_lambda_T(lambda["time"]);
    }

    // utils
    void init() {
      model_view_.set_data(data_);
      model_view_.init();
    }
    void configure_calibrator(const Rcpp::List & R_lambda) {
      calibrator_.configure_calibrator(R_lambda);
    }
    DVector<double> calibrate() {
      DVector<double> lambda_opt = calibrator_.fit(model_view_);
      model_view_.set_lambda(lambda_opt);
      return lambda_opt;
    }
    DVector<double> optimum() { return calibrator_.optimum(); }
    void solve() { model_view_.solve(); }

  public:
    void set_lambda(DVector<double> lambda) {
      model_view_.set_lambda_D(lambda[0]);
      // if(lambda.size() > 1)
      //  model_view_.set_lambda_T(lambda[1]);
    }

    // ** model specific methods ** //


    // ** calibrator specific methods ** //

    // gcv getters and setters
    std::vector<double> gcvs() const { 
      // fdapde_assert(calibration_strategy_ == Calibration::gcv);
      return calibrator_.gcvs();
    }
    std::vector<double> edfs() const {
      // fdapde_assert(calibration_strategy_ == Calibration::gcv);
      return calibrator_.edfs();
    }
    void set_step(double step) {
      // fdapde_assert(calibration_strategy_ == Calibration::gcv);
      calibrator_.set_step(step);
    }

    // kcv getters and setters
    const DVector<double>& avg_scores() const {
      // fdapde_assert(calibration_strategy_ == Calibration::kcv);
      return calibrator_.avg_scores();
    }
    const DVector<double>& std_scores() const {
      // fdapde_assert(calibration_strategy_ == Calibration::kcv);
      return calibrator_.std_scores();
    }
    const DMatrix<double>& scores() const {
      // fdapde_assert(calibration_strategy_ == Calibration::kcv);
      return calibrator_.scores();
    }
    void set_n_folds(std::size_t n_folds) {
      // fdapde_assert(calibration_strategy_ == Calibration::kcv);
      calibrator_.set_n_folds(n_folds);
    }
};

#endif // __R_REGRESSION_MODEL_H__