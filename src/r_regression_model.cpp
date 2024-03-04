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

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include "headers/r_regression_model.h"

// Rcpp modules definition
using cpp_rm = R_REGRESSION_MODEL;
using cpp_regression_model = R_REGRESSION_MODEL;
RCPP_MODULE(cpp_regression_model) {
    Rcpp::class_<cpp_regression_model>("cpp_regression_model") 
      .constructor<int, Rcpp::Environment, int, Rcpp::List>() 
      .method("get_calibration_strategy"   , &cpp_regression_model::get_calibration_strategy            )
      .method("f"                          , &cpp_regression_model::f                                   )
      .method("beta"                       , &cpp_regression_model::beta                                ) 
      .method("optimum"                    , &cpp_regression_model::optimum                             ) 
      .method("gcvs"                       , &cpp_regression_model::gcvs                                )
      .method("edfs"                       , &cpp_regression_model::edfs                                ) 
      .method("avg_scores"                 , &cpp_regression_model::avg_scores                          ) 
      .method("set_lambda"                 , (void (cpp_rm::*)(Rcpp::List))(&cpp_rm::set_lambda)        ) 
      .method("set_observations"           , &cpp_regression_model::set_observations                    )
      .method("set_covariates"             , &cpp_regression_model::set_covariates                      )
      .method("set_spatial_locations"      , &cpp_regression_model::set_spatial_locations               ) 
      .method("set_calibrator"             , &cpp_regression_model::set_calibrator                      ) 
      .method("init"                       , &cpp_regression_model::init                                )
      .method("calibrate"                  , &cpp_regression_model::calibrate                           )
      .method("solve"                      , &cpp_regression_model::solve                               );  
}
