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
#include "headers/r_fpls_r.h"
#include "headers/r_fpls_a.h"
#include "headers/r_fpls_sb.h"


using cpp_fpls_r_spaceonly_base = R_FPLS_R<SpaceOnly>;
using cpp_fpls_r_spaceonly = R_FPLS_R_SpaceOnly;
RCPP_MODULE(cpp_fpls_r_spaceonly) {
    // inherited methods
    Rcpp::class_<cpp_fpls_r_spaceonly_base>("cpp_fpls_r_spaceonly_base")
      // getters
      .method("Psi",                      &cpp_fpls_r_spaceonly_base::Psi                    )
      .method("R0",                       &cpp_fpls_r_spaceonly_base::R0                     )
      // setters
      .method("set_data",                 &cpp_fpls_r_spaceonly_base::set_data               )
      .method("set_spatial_locations",    &cpp_fpls_r_spaceonly_base::set_spatial_locations  )
      .method("set_ncomp",                &cpp_fpls_r_spaceonly_base::set_ncomp              )
      .method("set_lambda",               &cpp_fpls_r_spaceonly_base::set_lambda             )
      .method("set_solver",               &cpp_fpls_r_spaceonly_base::set_solver             )
      // getters
      .method("fitted",                   &cpp_fpls_r_spaceonly_base::fitted                 )
      .method("reconstructed",            &cpp_fpls_r_spaceonly_base::reconstructed          )
      .method("Y_space_directions",       &cpp_fpls_r_spaceonly_base::Y_space_directions     )
      .method("Y_loadings",               &cpp_fpls_r_spaceonly_base::Y_loadings             )
      .method("X_space_directions",       &cpp_fpls_r_spaceonly_base::X_space_directions     )
      .method("X_loadings",               &cpp_fpls_r_spaceonly_base::X_loadings             )
      .method("X_latent_scores",          &cpp_fpls_r_spaceonly_base::X_latent_scores        )
      .method("Y_latent_scores",          &cpp_fpls_r_spaceonly_base::Y_latent_scores        )
      .method("Beta",                     &cpp_fpls_r_spaceonly_base::Beta                   )
      // utilities
      .method("init",                     &cpp_fpls_r_spaceonly_base::init                   )
      .method("solve",                    &cpp_fpls_r_spaceonly_base::solve                  );
    // specific methods
    Rcpp::class_<R_FPLS_R_SpaceOnly>("cpp_fpls_r_spaceonly")
      .derives<cpp_fpls_r_spaceonly_base>("cpp_fpls_r_spaceonly_base")
      .constructor<Rcpp::Environment,
                   int,
                   Rcpp::List>();
}

using cpp_fpls_a_spaceonly_base = R_FPLS_A<SpaceOnly>;
using cpp_fpls_a_spaceonly = R_FPLS_A_SpaceOnly;
RCPP_MODULE(cpp_fpls_a_spaceonly) {
    // inherited methods
    Rcpp::class_<cpp_fpls_a_spaceonly_base>("cpp_fpls_a_spaceonly_base")
      // getters
      .method("Psi",                      &cpp_fpls_a_spaceonly_base::Psi                    )
      .method("R0",                       &cpp_fpls_a_spaceonly_base::R0                     )
      // setters
      .method("set_data",                 &cpp_fpls_a_spaceonly_base::set_data               )
      .method("set_spatial_locations",    &cpp_fpls_a_spaceonly_base::set_spatial_locations  )
      .method("set_ncomp",                &cpp_fpls_a_spaceonly_base::set_ncomp              )
      .method("set_lambda",               &cpp_fpls_a_spaceonly_base::set_lambda             )
      .method("set_solver",               &cpp_fpls_a_spaceonly_base::set_solver             )
      // getters
      .method("fitted",                   &cpp_fpls_a_spaceonly_base::fitted                 )
      .method("reconstructed",            &cpp_fpls_a_spaceonly_base::reconstructed          )
      .method("Y_space_directions",       &cpp_fpls_a_spaceonly_base::Y_space_directions     )
      .method("Y_loadings",               &cpp_fpls_a_spaceonly_base::Y_loadings             )
      .method("X_space_directions",       &cpp_fpls_a_spaceonly_base::X_space_directions     )
      .method("X_loadings",               &cpp_fpls_a_spaceonly_base::X_loadings             )
      .method("X_latent_scores",          &cpp_fpls_a_spaceonly_base::X_latent_scores        )
      .method("Y_latent_scores",          &cpp_fpls_a_spaceonly_base::Y_latent_scores        )
      // utilities
      .method("init",                     &cpp_fpls_a_spaceonly_base::init                   )
      .method("solve",                    &cpp_fpls_a_spaceonly_base::solve                  );
    // specific methods
    Rcpp::class_<R_FPLS_A_SpaceOnly>("cpp_fpls_a_spaceonly")
      .derives<cpp_fpls_a_spaceonly_base>("cpp_fpls_a_spaceonly_base")
      .constructor<Rcpp::Environment,
                   int,
                   Rcpp::List>();
}


using cpp_fpls_sb_spaceonly_base = R_FPLS_SB<SpaceOnly>;
using cpp_fpls_sb_spaceonly = R_FPLS_SB_SpaceOnly;
RCPP_MODULE(cpp_fpls_sb_spaceonly) {
    // inherited methods
    Rcpp::class_<cpp_fpls_sb_spaceonly_base>("cpp_fpls_sb_spaceonly_base")
      // getters
      .method("Psi",                      &cpp_fpls_sb_spaceonly_base::Psi                    )
      .method("R0",                       &cpp_fpls_sb_spaceonly_base::R0                     )
      // setters
      .method("set_data",                 &cpp_fpls_sb_spaceonly_base::set_data               )
      .method("set_spatial_locations",    &cpp_fpls_sb_spaceonly_base::set_spatial_locations  )
      .method("set_ncomp",                &cpp_fpls_sb_spaceonly_base::set_ncomp              )
      .method("set_lambda",               &cpp_fpls_sb_spaceonly_base::set_lambda             )
      .method("set_solver",               &cpp_fpls_sb_spaceonly_base::set_solver             )
      // getters
      .method("fitted",                   &cpp_fpls_sb_spaceonly_base::fitted                 )
      .method("reconstructed",            &cpp_fpls_sb_spaceonly_base::reconstructed          )
      .method("Y_space_directions",       &cpp_fpls_sb_spaceonly_base::Y_space_directions     )
      .method("Y_loadings",               &cpp_fpls_sb_spaceonly_base::Y_loadings             )
      .method("X_space_directions",       &cpp_fpls_sb_spaceonly_base::X_space_directions     )
      .method("X_loadings",               &cpp_fpls_sb_spaceonly_base::X_loadings             )
      .method("X_latent_scores",          &cpp_fpls_sb_spaceonly_base::X_latent_scores        )
      .method("Y_latent_scores",          &cpp_fpls_sb_spaceonly_base::Y_latent_scores        )
      // utilities
      .method("init",                     &cpp_fpls_sb_spaceonly_base::init                   )
      .method("solve",                    &cpp_fpls_sb_spaceonly_base::solve                  );
    // specific methods
    Rcpp::class_<R_FPLS_SB_SpaceOnly>("cpp_fpls_sb_spaceonly")
      .derives<cpp_fpls_sb_spaceonly_base>("cpp_fpls_sb_spaceonly_base")
      .constructor<Rcpp::Environment,
                   int,
                   Rcpp::List>();
}