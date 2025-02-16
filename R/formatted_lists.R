## TODOs:
## TODO: possible fault when setting a centering device in functional models
##       - NULL: default values are used
##       - logic: default values are used if true
##       - other: used as centering device
##       In this lasta case we have to check that it is a proper centering
##       device initialization list


# enumerators ----

match.list <- function(options_list, x) {
  if (is.numeric(x)) {
    return(names(options_list)[options_list == x])
  } else {
    return(options_list[[x]])
  }
}

Calibration <- function(x) {
  calibration_options <- list(
    "off" = 0,
    "gcv" = 1,
    "kcv" = 2
  )
  match.list(calibration_options, x)
}

Sampling <- function(x) {
  sampling_options <- list(
    "mesh_nodes" = 0,
    "pointwise" = 1,
    "areal" = 2
  )
  match.list(sampling_options, x)
}

Regularization <- function(x) {
  regularization_options <- list(
    "SpaceOnly" = 0,
    "SpaceTimeSeparable" = 1,
    "SpaceTimeParabolic" = 2
  )
  match.list(regularization_options, x)
}

SolutionPolicy <- function(x) {
  solution_policy_options <- list(
    "sequential" = 0,
    "monolithic" = 1
  )
  match.list(solution_policy_options, x)
}

RegressionModel <- function(x) {
  regression_model_options <- list(
    "SRPDE" = 0
  )
  match.list(regression_model_options, x)
}


# hyperparameters ----

#' @export
hyperparameters <- function(space, time = 0) {
  lambda <- list(
    type = "fdaPDE_hyperparameters",
    space = space,
    time = time
  )
  return(lambda)
}

# data ----

#' @export
spatial_data <- function(domain, observations,
                         locations = NULL,
                         covariates = NULL) {
  data <- list(
    type = "fdaPDE_spatialData",
    domain = domain,
    locations = locations,
    observations = observations,
    covariates = covariates
  )
  return(data)
}

#' @export
functional_data <- function(domain, X,
                            locations = NULL,
                            w = NULL) {
  data <- list(
    type = "fdaPDE_functionalData",
    domain = domain,
    locations = locations,
    X = X,
    w = w
  )
  return(data)
}

#' @export
functional_regression_data <- function(domain, Y, X,
                                       locations = NULL) {
  data <- list(
    type = "fdaPDE_functionalRegressionData",
    domain = domain,
    locations = locations,
    Y = Y,
    X = X
  )
  return(data)
}

# calibrators ----

off <- function(lambda = hyperparameters(space = 1e-4)) {
  ## calibration strategy specific parameters
  off_params <- NULL
  ## init list assembly
  off_init_list <- list(
    type = "fdaPDE_calibrator",
    ## calibration strategy name
    name = "off",
    ## lambda
    lambda = lambda,
    ## parameters
    parameters = off_params
  )
  return(off_init_list)
}

#' @export
gcv <- function(lambda = hyperparameters(space = 10^seq(-9, 0, by = 0.5)),
                optimizer = c("grid"),
                edf_computation = c("stochastic", "exact"),
                seed = -1, mc_samples = 100,
                max_iter = NULL, step = NULL,
                tolerance = NULL) {
  ## calibration strategy specific parameters
  gcv_params <- list(
    ## properties related to the computation of the expected degrees of freedom
    edf_computation = match.arg(edf_computation),
    seed = seed, ## meaningfull only for stochastic edf_computation
    mc_samples = mc_samples, ## meaningfull only for stochastic edf_computation
    ## properties related to the optimization algorithm
    optimizer = match.arg(optimizer),
    max_iter = max_iter,
    step = step,
    tolerance = tolerance
  )
  ## init list assembly
  gcv_init_list <- list(
    type = "fdaPDE_calibrator",
    ## calibration strategy name
    name = "gcv",
    ## lambda grid
    lambda = lambda,
    ## parameters
    parameters = gcv_params
  )
  return(gcv_init_list)
}

#' @export
kcv <- function(lambda = hyperparameters(space = 10^seq(-9, 0, by = 1)),
                n_folds = 10, shuffle = TRUE, seed = NULL) {
  ## calibration strategy specific parameters
  kcv_params <- list(
    n_folds = n_folds,
    shuffle = shuffle,
    seed = seed
  )
  ## init list assembly
  kcv_init_list <- list(
    type = "fdaPDE_calibrator",
    ## calibration strategy name
    name = "kcv",
    ## lambda grid
    lambda = lambda,
    ## parameters
    parameters = kcv_params
  )
  return(kcv_init_list)
}

parse_calibrator <- function(calibrator) {
  if (calibrator$type == "fdaPDE_hyperparameters") {
    return(off(lambda = calibrator))
  } else {
    ## TODO: check if this actually is a calibrator
    return(calibrator)
  }
}

#' @export
calibration <- function(strategy = c("off", "gcv", "kcv"), ...) {
  calibrator_params <- switch(match.arg(strategy),
    "off" = {
      off(...)
    },
    "gcv" = {
      gcv(...)
    },
    "kcv" = {
      kcv(...)
    },
    stop("The selected calibrator does not exist.")
  )
  return(calibrator_params)
}


# models ----


## Regularized SVD ----

### sequential ----

sequential_params <- function(tolerance = 1e-6, max_iter = 20) {
  parameters <- list(
    tolerance = tolerance,
    max_iter = max_iter
  )
  return(parameters)
}

#' @export
sequential <- function(calibrator = off(), ...) {
  sequential_init_list <- list(
    policy = "sequential",
    parameters = sequential_params(...),
    calibrator = parse_calibrator(calibrator)
  )
  return(sequential_init_list)
}

### monolithic ----

monolithic_params <- function() {
  return(NULL)
}

#' @export
monolithic <- function(calibrator = hyperparameters(space = 1e-4), ...) {
  ## check calibrator availability
  if (parse_calibrator(calibrator)$name != "off") {
    stop("Only off calibrator is available for this solution policy")
  }
  ##
  monolithic_init_list <- list(
    policy = "monolithic",
    parameters = monolithic_params(...),
    calibrator = parse_calibrator(calibrator)
  )
  return(monolithic_init_list)
}

### generic ---

#' @export
RSVD <- function(policy = c("sequential", "monolithic"), ...) {
  RSVD_params <- switch(match.arg(policy),
    "sequential" = {
      sequential(...)
    },
    "monolithic" = {
      monolithic(...)
    },
    stop("The selected policy does not exist.")
  )
  return(RSVD_params)
}


## spatial ----


### regression ----

#### SRPDE ----

SRPDE_params <- function() {
  return(NULL)
}

#### generic ----

#' @export
smoothing <- function(name = c("SRPDE"),
                      penalty = simple_laplacian_penalty(),
                      calibrator = gcv(),
                      ...) {
  ## method specific parameters initialization
  parameters <- switch(match.arg(name),
    "SRPDE" = {
      SRPDE_params(...)
    },
    stop("The selected smoother does not exist.")
  )
  ## init list assembly
  smoother_init_list <- list(
    name = match.arg(name),
    ## regularization
    penalty = penalty,
    ## method specific parameters
    parameters = parameters,
    ## calibration
    calibrator = parse_calibrator(calibrator)
  )
  return(smoother_init_list)
}


## functional ----

parse_center <- function(center, as.option = FALSE) {
  if (!as.option) {
    if (is.null(center) || is.logical(center)) {
      return(centering())
    } else {
      return(center)
    }
  } else {
    if (is.null(center)) {
      return(TRUE)
    } ## fPCA defaults center to TRUE
    else if (is.logical(center)) {
      return(center)
    } else {
      return(TRUE)
    }
  }
}

### fCentering ----

centering <- function(...) {
  fCentering_init(...)
}
fCentering_init <- function(calibrator = gcv(),
                            penalty = simple_laplacian_penalty(),
                            smoother = smoothing()) {
  smoother$calibrator <- parse_calibrator(calibrator)
  fCentering_init_list <- list(
    name = "fCentering",
    ## regularization
    penalty = penalty,
    ## analysis
    smoother = smoother
  )
  return(fCentering_init_list)
}

### fPCA ----

fPCA_params <- function(n_pc = 3) {
  parameters <- list(
    n_pc = n_pc
  )
  return(parameters)
}

fPCA_init <- function(penalty = simple_laplacian_penalty(),
                      center = NULL, ## also defaults centering lambda
                      solver = sequential(), ## also defaults solver's lambda
                      ...) {
  fPCA_init_list <- list(
    name = "fPCA",
    ## regularization
    penalty = penalty,
    ## centering
    center = parse_center(center),
    ## analsysis
    solver = solver,
    ## method parameters
    parameters = fPCA_params(...),
    ## options
    CENTER = parse_center(center, as.option = TRUE)
  )
  return(fPCA_init_list)
}

### fPLS ----

fPLS_params <- function(n_comp = 3) {
  parameters <- list(
    n_comp = n_comp
  )
  return(parameters)
}

fPLS_init <- function(penalty = simple_laplacian_penalty(),
                      center = NULL,
                      solver = sequential(),
                      smoother = smoothing(),
                      mode = "fPLS-R",
                      ...) {
  fPLS_init_list <- list(
    name = "fPLS",
    ## regularization
    penalty = penalty,
    ## centering
    center = parse_center(center),
    ## analysis
    solver = solver,
    smoother = smoother,
    ## method parameters
    parameters = fPLS_params(...),
    ## options
    MODE = mode,
    CENTER = parse_center(center, as.option = TRUE)
  )
  return(fPLS_init_list)
}


### generic ----

functional_model <- function(name = c("fPCA", "fPLS"),
                             penalty = simple_laplacian_penalty(),
                             center = NULL,
                             solver = sequential(),
                             ...) {
  functional_model_init_list <- switch(match.arg(name),
    "fCentering" = {
      fPCA_init(
        penalty = penalty,
        ...
      )
    },
    "fPCA" = {
      fPCA_init(
        penalty = penalty,
        solver = solver,
        center = center,
        ...
      )
    },
    "fPLS" = {
      fPLS_init(
        penalty = penalty,
        solver = solver,
        center = center,
        ...
      )
    },
    stop("The selected functional model does not exist.")
  )
  return(functional_model_init_list)
}
