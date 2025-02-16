## TODOs:
## TODO fdaPDE_Regression_Model: implement regression data sanity checks:
##      - check that data of type spatial_data have been passed
##      - coherence between locations and observations dimensions
##      - formula parser and check that everything has been passed available
## ???? does the penalty need any sanity check?
## TODO fdaPDE_Regression_Model: possible fault when passing smoother_params
##                               if the user does not use the function
##                               SmootherName_params the parameters are not
##                               checked and passed directly to the constructor
## TODO: add sanity checks on objects provided to the fit method,
##       now they are assumed to be generated using the suitable
##       formatted_lists functions.

## fdaPDE_Regression_Model models ----

## RMK: spatial models == regression models (by now)
## This might need to be changed when we will include
## other methods as density estimation, ...

## interface for a generic regression model
fdaPDE_Regression_Model <- R6::R6Class(
  classname = "fdaPDE_Regression_Model",
  inherit = fdaPDE_Base_Model,
  public = list(
    ## regression data
    observations = NULL,
    covariates = NULL
  ),
  private = list(
    ## smoother_init_list
    regression_model = NULL,
    get_regression_model = function() {
      return(private$regression_model)
    },
    set_regression_model = function(regression_model) {
      private$regression_model <- regression_model
    },
    ## initialize utils
    sanity_check_data = function(formula = NULL, data) {
      ## domain informations
      super$sanity_check_domain(data)
      ## regression data
      private$sanity_check_regression_data(formula, data)
    },
    sanity_check_regression_data = function(formula, data) {
      self$observations <- as.matrix(data$observations)
      if (!is.null(data$covariates)) {
        self$covariates <- as.matrix(data$covariates)
      }
      ## TODO: check coherence between observations and locations dimensions
    },
    init_model_and_calibrator = function(smoother_init_list) {
      ## init cpp_model & cpp_calibrator
      super$cpp_model <- regression_models_factory(
        self$domain,
        self$model_traits,
        smoother_init_list
      )
      ## set locations
      super$cpp_model$set_spatial_locations(self$locations)
      super$display("  Locations have been set.")
      ## set data
      super$cpp_model$set_observations(self$observations)
      super$display("  Observations have been set.")
      if (!is.null(self$covariates)) {
        super$cpp_model$set_covariates(self$covariates)
        super$display("  Covariates have been set.")
      }
      ## set calibrator
      private$set_calibrator(smoother_init_list$calibrator)
    },
    ## fit utils
    calibrate = function() {
      ## model initialization (necessary for GCV and KCV calibration strategies)
      super$cpp_model$init()
      ## calibration
      lambda_cpp <- super$cpp_model$calibrate()
      lambda_opt <- hyperparameters(lambda_cpp[1], lambda_cpp[2])
      ## save calibrator's results
      self$results$calibrator <- export_calibrator_results(super$cpp_model)
    },
    ## cpp_model interfaces
    set_lambda = function(lambda) {
      super$cpp_model$set_lambda(lambda)
    },
    set_calibrator = function(calibrator) {
      super$cpp_model$set_calibrator(
        Calibration(calibrator$name),
        calibrator$parameters,
        calibrator$lambda
      )
    },
    get_f = function() {
      return(as.matrix(super$cpp_model$f()))
    },
    get_beta = function() {
      return(as.matrix(super$cpp_model$beta()))
    },
    fit = function() {
      super$cpp_model$init()
      super$cpp_model$solve()
    }
  )
)


### SRPDE model implementation ----

SRPDE_Class <- R6::R6Class(
  classname = "SRPDE_Class",
  inherit = fdaPDE_Regression_Model,
  public = list(
    ## constructor
    initialize = function(formula, data, smoother, VERBOSE) {
      super$set_verbosity(VERBOSE)
      super$display("\n\nSRPDE model\n")
      ## inputs sanity check
      super$display("- Performing inputs sanity check")
      super$sanity_check_data(formula, data)
      ## model and calibrator initialization
      super$display("- Initializing the model and the calibrator")
      super$set_regression_model(smoother)
      super$init_model_and_calibrator(super$get_regression_model())
    },
    ## utils
    fit = function(calibrator = NULL) {
      regression_model <- super$get_regression_model()
      ## clear previous results
      self$results <- list()
      ## calibrator update
      if (!is.null(calibrator)) {
        regression_model$calibrator <- parse_calibrator(calibrator) ## to keep the init list updated
        ## save updates
        super$set_regression_model(regression_model)
        ## initialize a new calibrator
        super$set_calibrator(super$get_regression_model()$calibrator)
      }
      ## calibration
      super$display("- Calibrating the model")
      super$calibrate()
      ## fit
      super$display("- Fitting the model")
      super$fit()
      ## save results
      super$display("- Saving the results")
      self$results$f <- super$get_f()
      self$results$beta <- super$get_beta()
    }
  ),
  private = list()
)

SRPDE_pro <- function(formula = NULL, data,
                      smoother = smoothing(),
                      VERBOSE = FALSE) {
  model <- SRPDE_Class$new(
    formula = formula,
    data = data,
    smoother = smoother,
    VERBOSE = VERBOSE
  )
  return(model)
}

#' @export
SRPDE <- function(formula = NULL, data,
                  penalty = simple_laplacian_penalty(),
                  smoother_options = SRPDE_params(),
                  calibrator = off(),
                  VERBOSE = FALSE) {
  ## overwriting smoother defaults with the required ones
  ## necessary to guarantee the consistency with the other statistical models
  ## while providing a clear interface to the final user
  smoother <- smoothing()
  smoother$penalty <- penalty
  smoother$calibrator <- parse_calibrator(calibrator)
  smoother$smoother_params <- smoother_options
  ## return wrapped model
  return(SRPDE_pro(
    formula = formula,
    data = data,
    smoother = smoother,
    VERBOSE = VERBOSE
  ))
}

# #' @export
# GSRPDE <- function(formula = NULL, data,
#                    penalty = simple_laplacian_penalty(),
#                    smoother_params = GSRPDE_params(),
#                    calibrator = off(),
#                    VERBOSE = FALSE) {
#   model <- SRPDE_Class$new(
#     formula = formula,
#     data = data,
#     penalty = penalty,
#     smoother_params = NULL,
#     calibrator = parse_calibrator(calibrator)),
#     VERBOSE = VERBOSE
#   )
#   return(model)
# }
