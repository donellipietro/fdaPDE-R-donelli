## TODOs:
## TODO regression_model_factory: let this factory return also a calibrator

# models ----

# model_traits
# - regularization_type
# - sampling_type

# smoother_init_list
# - name
# - penalty
#   - pde_callable
#   - parameters
# - parameters
# ...

regression_models_factory <- function(domain, model_traits, smoother_init_list) {
  ## pde initialization
  pde_callable <- smoother_init_list$penalty$pde_callable
  pde_parameters <- smoother_init_list$penalty$parameters
  pde <- pde_callable(domain, pde_parameters)
  ## model initialization
  cpp_model <- new(
    ## cpp module
    cpp_regression_model,
    ## contructor's arguments
    RegressionModel(smoother_init_list$name),
    pde,
    Sampling(model_traits$sampling_type),
    smoother_init_list$parameters
  )
  return(cpp_model)
}


# model_traits
# - regularization_type
# - sampling_type

# fm_init_list
# - name
# - penalty
#   - pde_callable
#   - parameters
# - parameters
# ...

functional_models_factory <- function(domain, model_traits, fm_init_list) {
  ## pde initialization
  pde_callable <- fm_init_list$penalty$pde_callable
  parameters <- fm_init_list$penalty$parameters
  pde <- pde_callable(domain, parameters)
  ## model initialization
  functional_model <- switch(fm_init_list$name,
    "fCentering" <- {
      cpp_fCentering <- new(
        ## cpp module
        cpp_center,
        ## contructor's arguments
        RegressionModel(fm_init_list$name),
        pde,
        Sampling(model_traits$sampling_type),
        fm_init_list$parameters
      )
      return(cpp_fCentering)
    },
    "fPCA" = {
      ## model initialization
      switch(model_traits$regularization_type,
        "SpaceOnly" = {
          ## statistical model initialization
          cpp_fPCA <- new(
            cpp_fpca_spaceonly,
            ## contructor's arguments
            pde,
            Sampling(model_traits$sampling_type),
            fm_init_list$parameters
          )
          return(cpp_fPCA)
        }
      )
    },
    "fPLS" = {
      ## model initialization
      switch(model_traits$regularization_type,
        "SpaceOnly" = {
          cpp_fpls_mode <- switch(fm_init_list$MODE,
            "fPLS-R" = {
              cpp_fpls_r_spaceonly
            },
            "fPLS-A" = {
              cpp_fpls_a_spaceonly
            },
            "fPLS-SB" = {
              cpp_fpls_sb_spaceonly
            },
            {
              stop("This fPLS mode is not implemented.")
            }
          )
          ## statistical model initialization
          cpp_fPLS <- new(
            ## cpp module
            cpp_fpls_mode,
            ## contructor's arguments
            pde,
            Sampling(model_traits$sampling_type),
            fm_init_list$parameters
          )
          return(cpp_fPLS)
        }
      )
    }
  )
  return(functional_model)
}
