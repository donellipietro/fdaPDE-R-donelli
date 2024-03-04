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


#ifndef __R_CALIBRATOR_H__
#define __R_CALIBRATOR_H__

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

// calibration strategy enumerator
#include <fdaPDE/calibration/symbols.h>
using CalibrationStrategy = fdapde::calibration::Calibration;

// calibrators classes
#include <fdaPDE/calibration/off.h> 
#include <fdaPDE/calibration/gcv.h> 
#include <fdaPDE/models/regression/gcv.h> 
#include <fdaPDE/calibration/kfold_cv.h>
#include <fdaPDE/calibration/rmse.h>
using fdapde::calibration::RMSE;
using fdapde::models::GCV;
using fdapde::calibration::Calibrator;
using fdapde::calibration::Calibration;

// optimization
#include <fdaPDE/optimization.h>
using fdapde::models::EDFStrategy;
using fdapde::models::StochasticEDF;
using fdapde::models::ExactEDF;

// regression models
#include <fdaPDE/models/regression/regression_type_erasure.h>
using fdapde::models::RegressionView;

/* TODOs:
// TODO gcv: generalize the GCV to work with other optimizers, by now the optimizer is forced to be opt_grid.
// TODO gcv: generalize the GCV to work with other EDF strategies, by now the EDF strategies is forced to be stochastic.
// TODO kcv: the fit method needs to be generalized, by now the ScoreType is set ot RMSE
// TODO gcv/kcv: the set_lambda method needs to be generalized to space-time
*/

class R_CALIBRATOR {
  private:

    Calibration calibration_strategy_;
    fdapde::calibration::Off calibrator_off_;
    fdapde::calibration::GCV<void> calibrator_gcv_;
    fdapde::calibration::KCV calibrator_kcv_;
    Calibrator<RegressionView<void>> configured_calibrator_;
    std::vector<DVector<double>> lambda_grid_;
    
  public:
    // constructor
    R_CALIBRATOR() : calibration_strategy_(Calibration::off) {};
    R_CALIBRATOR(Calibration calibration_strategy, const Rcpp::List & calibrator_params) : calibration_strategy_(calibration_strategy) {
      switch (calibration_strategy_) {
        case Calibration::off : {
          break;
        }
        case Calibration::gcv : {
          // grid optimizer (fixed to grid)
          fdapde::core::Grid<fdapde::Dynamic> opt;
          // EDF strategy (set to stochastic)
          StochasticEDF edf;
          edf.set_seed(Rcpp::as<int>(calibrator_params["seed"]));
          edf.set_n_mc_samples(Rcpp::as<int>(calibrator_params["mc_samples"]));
          // calibrator
          calibrator_gcv_ = fdapde::calibration::GCV {opt, edf};
          break;
        }
        case Calibration::kcv : {
          // number of folds
          if (calibrator_params["n_folds"] == R_NilValue){
            calibrator_kcv_ = fdapde::calibration::KCV {};
            return;
          }
          std::size_t n_folds = calibrator_params["n_folds"];
          // folds randomization
          bool shuffle = calibrator_params["shuffle"];
          if (calibrator_params["seed"] != R_NilValue){
            calibrator_kcv_ = fdapde::calibration::KCV{n_folds, calibrator_params["seed"], shuffle};
            return;
          } else{
            calibrator_kcv_ = fdapde::calibration::KCV {n_folds, shuffle};
            return;
          }
          break;
        }
      }
    };

    // getters
    int get_calibration_strategy() { return calibration_strategy_;}
    DVector<double> optimum() { return configured_calibrator_.optimum(); }
    Calibrator<RegressionView<void>> & get_configured_calibrator() {
      return configured_calibrator_;
    }

    // utils
    void parse_R_lambda(Rcpp::List R_lambda){
      Rcpp::NumericVector lambda_D_grid = R_lambda["space"];
      // move Rcpp::NumericVector into something understandable by cpp layer
      lambda_grid_.resize(lambda_D_grid.size());
      for(std::size_t i = 0; i < lambda_grid_.size(); ++i){
        lambda_grid_[i].resize(1);
        lambda_grid_[i][0] = lambda_D_grid[i];
      }
    }
    Calibrator<RegressionView<void>> configure_calibrator(const Rcpp::List & R_lambda){
      switch (calibration_strategy_) {
        case Calibration::off : 
          return configure_calibrator_off(R_lambda);
        case Calibration::gcv :
          return configure_calibrator_gcv(R_lambda);
        case Calibration::kcv :
          return configure_calibrator_kcv(R_lambda);
      }
    }
    DVector<double> fit(RegressionView<void> & model_view) {
      switch (calibration_strategy_) {
        case Calibration::off : 
          return fit_off(model_view);
        case Calibration::gcv :
          return fit_gcv(model_view);
        case Calibration::kcv :
          return fit_kcv(model_view);
      }
    }

    // ** calibrator specific methods ** //

    // gcv getters and setters
    std::vector<double> gcvs() const { 
      fdapde_assert(calibration_strategy_ == Calibration::gcv);
      return calibrator_gcv_.gcvs();
    }
    std::vector<double> edfs() const {
      fdapde_assert(calibration_strategy_ == Calibration::gcv);
      return calibrator_gcv_.edfs();
    }
    void set_step(double step) {
      fdapde_assert(calibration_strategy_ == Calibration::gcv);
      calibrator_gcv_.set_step(step);
    }

    // kcv getters and setters
    const DVector<double>& avg_scores() const {
      fdapde_assert(calibration_strategy_ == Calibration::kcv);
      return calibrator_kcv_.avg_scores();
    }
    const DVector<double>& std_scores() const {
      fdapde_assert(calibration_strategy_ == Calibration::kcv);
      return calibrator_kcv_.std_scores();
    }
    const DMatrix<double>& scores() const {
      fdapde_assert(calibration_strategy_ == Calibration::kcv);
      return calibrator_kcv_.scores();
    }
    void set_n_folds(std::size_t n_folds) {
      fdapde_assert(calibration_strategy_ == Calibration::kcv);
      calibrator_kcv_.set_n_folds(n_folds);
    }

  private:

    // off utils
    Calibrator<RegressionView<void>> configure_calibrator_off(const Rcpp::List & R_lambda){
      parse_R_lambda(R_lambda); // it initializes lambda_grid_
      DVector<double> lambda = lambda_grid_.front();
      configured_calibrator_ = calibrator_off_(lambda);
      return get_configured_calibrator();
    }
    DVector<double> fit_off(RegressionView<void> model_view){
      return calibrator_off_.fit(model_view);
    }

    // gcv utils 
    Calibrator<RegressionView<void>> configure_calibrator_gcv(const Rcpp::List & R_lambda){
      parse_R_lambda(R_lambda); // it initializes lambda_grid_
      configured_calibrator_ = calibrator_gcv_(lambda_grid_);
      return get_configured_calibrator();;
    }
    DVector<double> fit_gcv(RegressionView<void> & model_view){
      return calibrator_gcv_.fit(model_view, lambda_grid_);
    }

    // kcv utils
    Calibrator<RegressionView<void>> configure_calibrator_kcv(const Rcpp::List & R_lambda){
      parse_R_lambda(R_lambda); // it initializes lambda_grid_
      configured_calibrator_ = calibrator_kcv_(lambda_grid_, RMSE());
      return get_configured_calibrator();;
    }
    DVector<double> fit_kcv(RegressionView<void> & model_view){
      return calibrator_kcv_.fit(model_view, lambda_grid_, RMSE(model_view));
    }

};

#endif   // __R_CALIBRATOR_H__