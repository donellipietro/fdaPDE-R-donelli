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

#ifndef __R_CENTER_H__
#define __R_CENTER_H__

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

// utils
#include <fdaPDE/utils/symbols.h>
#include "r_pde.h"

// models
#include <fdaPDE/models.h>
#include "r_regression_model.h"
using fdapde::models::Sampling;
using fdapde::models::RegressionView;

// calibrators
#include "r_calibrator.h"
#include <fdaPDE/calibration/calibration_base.h>
using fdapde::calibration::Calibrator;

// center method
#include <fdaPDE/models/functional/center.h>
using fdapde::models::CenterReturnType;


// implementation of CENTER wrapper
class R_CENTER {

private:
    // regression model
    R_REGRESSION_MODEL regression_model_;
    // data
    DMatrix<double> X_;
    DMatrix<double> w_;
    bool weighted_ = FALSE;
    // results
    CenterReturnType results_;

public:

    // constructor
    R_CENTER(int regression_model,
             const Rcpp::Environment & pde, int sampling_type, const Rcpp::List & smoother_params) : 
             regression_model_(regression_model, pde, sampling_type, smoother_params) {}

    // getters
    DMatrix<double> centered() {return results_.fitted; }
    DMatrix<double> mean() {return results_.mean; }

    // setters
    void set_data(const DMatrix<double>& X) { X_ = X; }
    void set_weights(const DVector<double>& w) {
      w_ = w;
      weighted_ = TRUE;
    }
    void set_spatial_locations(const DMatrix<double>& locs) { regression_model_.set_spatial_locations(locs); }
    void set_calibrator(int calibration_strategy, const Rcpp::List & calibrator_params, const Rcpp::List & R_lambda) {
      regression_model_.set_calibrator(Calibration(calibration_strategy), calibrator_params, R_lambda);
    }

    // utilities
    void init() { return; }
    void configure_calibrator(const Rcpp::List & R_lambda) {
      regression_model_.configure_calibrator(R_lambda);
    }
    void solve() {
      if(weighted_)
        results_ = center(X_, w_, regression_model_.get_model_view(), regression_model_.get_configured_calibrator());
      else 
        results_ = center(X_, regression_model_.get_model_view(), regression_model_.get_configured_calibrator());
    }

};

#endif // __R_CENTER_H__
