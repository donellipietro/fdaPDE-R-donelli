## functional centering models ----

fCentering_class <- R6::R6Class(
  classname = "fCentering",
  inherit = fdaPDE_Base_Model,
  public = list(
    ## data
    X = NULL, ## data to center
    w = NULL, ## weights
    ## constructor
    initialize = function(data, center, VERBOSE) {
      super$set_verbosity(VERBOSE)
      super$display("\n\nCentering model\n")
      ## inputs sanity check
      super$display("- Performing inputs sanity check")
      private$sanity_check_data(data)
      ## centering initialization
      super$display("- Initializing the centering model")
      private$center_model <- center
      private$init_centering(private$center_model)
    },
    ## utilities
    fit = function(calibrator = NULL, smoother = NULL) {
      ## clear previous results
      self$results <- list()
      if (!is.null(smoother)) { ## smoother update if necessary
        super$display("- set custom smoother (this will overwrite the old calibrator)")
        private$center_model <- smoother ## to keep center updated
        ## check if also the calibrator needs changes
        if (!is.null(calibrator)) {
          super$display("- Setting the new calibrator")
          private$center_model$calibrator <- parse_calibrator(calibrator) ## to keep center updated
        }
        private$init_centering(private$center_model)
      } else if (!is.null(calibrator)) { ## calibrator update if necessary
        super$display("- Setting the new calibrator")
        private$center_model$calibrator <- parse_calibrator(calibrator) ## to keep center updated
        super$cpp_model$set_calibrator(
          Calibration(private$center_model$calibrator$name),
          private$center_model$calibrator$parameters,
          private$center_model$calibrator$lambda
        )
      }
      # centering by using the loaded smoother and calibrator
      super$display("- Centering")
      super$cpp_model$init()
      super$cpp_model$solve()
      ## save fit results
      super$display("- Saving the results")
      self$results$mean <- as.matrix(super$cpp_model$mean())
      self$results$centered <- as.matrix(super$cpp_model$centered())
      self$results$calibration$lambda_opt <- hyperparameters(super$cpp_model$optimum()[1])
    }
  ),
  private = list(
    ## center_init_list
    center_model = NULL,
    ## initialize utils
    sanity_check_data = function(data) {
      ## domain informations
      super$sanity_check_domain(data)
      ## regression data
      private$sanity_check_centering_data(data)
    },
    sanity_check_centering_data = function(data) {
      self$X <- as.matrix(data$X)
      if (!is.null(data$w)) {
        self$w <- as.matrix(data$w)
      }
      ## TODO: check coherence between X and locations dimensions
      ## TODO: check coherence between X and w dimensions
    },
    init_centering = function(center_init_list) {
      ## init all the cpp objects
      super$cpp_model <- functional_models_factory(domain, self$model_traits, center_init_list$smoother)
      ## set data into the cpp model
      super$cpp_model$set_data(self$X)
      super$display("  Data have been set.")
      if (!is.null(self$w)) {
        super$cpp_model$set_weights(self$w)
        super$display("  Weights have been set.")
      }
      ## set locations into the smoother
      super$cpp_model$set_spatial_locations(self$locations)
      super$display("  Locations have been set.")
      ## calibrator init
      super$display("- Initializing the calibrator")
      super$cpp_model$set_calibrator(
        Calibration(center_init_list$smoother$calibrator$name),
        center_init_list$smoother$calibrator$parameters,
        center_init_list$smoother$calibrator$lambda
      )
    }
  )
)

## library friendly interface
fCentering_pro <- function(data, center = centering, VERBOSE = FALSE) {
  model <- fCentering_class$new(
    data = data,
    center = center,
    VERBOSE = VERBOSE
  )
  return(model)
}

## user friendly interfcace
#' @export
fCentering <- function(data,
                       penalty = simple_laplacian_penalty(),
                       smoother = smoothing(),
                       calibrator = hyperparameters(space = 1e-4),
                       VERBOSE = FALSE) {
  center_init_list <- fCentering_init(
    penalty = penalty,
    smoother = smoother,
    calibrator = parse_calibrator(calibrator)
  )
  return(fCentering_pro(data, center_init_list, VERBOSE))
}

## FunctionalModel models

fdaPDE_Functional_Model <- R6::R6Class(
  classname = "fdaPDE_Functional_Model",
  inherit = fdaPDE_Base_Model,
  public = list(
    ## regression data
    data = list(),
    ## options
    CENTER = NULL,
    init = function() {
      super$cpp_model$init()
    },
    Psi = function() {
      return(self$Psi_mat)
    },
    R0 = function() {
      return(self$R0_mat)
    },
    evaluate = function(f) {
      return(as.matrix(self$Psi() %*% f))
    },
    fitted_loadings = function() {
      return(self$evaluate(self$results$loadings))
    }
  ),
  private = list(
    ## functional_model_init_list
    functional_model = NULL,
    ## centering device instance
    R_center = NULL,
    ## initialize utils
    sanity_check_data = function(data) {
      ## domain informations
      super$sanity_check_domain(data)
      ## functional data
      private$sanity_check_functional_data(data)
    },
    sanity_check_functional_data = function(data) {
      ## X data
      self$data$X <- as.matrix(data$X)
      ## Y data (if any)
      if ("Y" %in% names(data)) {
        self$data$Y <- as.matrix(data$Y)
      }
    },
    init_model = function(functional_model_init_list) {
      ## set problem specific regularization_type & sampling_type
      functional_model_init_list$regularization_type <- self$model_traits$regularization_type
      functional_model_init_list$sampling_type <- self$model_traits$sampling_type
      ## save fPCA_model
      private$functional_model <- functional_model_init_list
      ## funcitona_model init
      super$cpp_model <- functional_models_factory(
        self$domain,
        self$model_traits,
        private$functional_model
      )
      ## set the selected solver
      private$set_solver()
      ## set locations
      super$cpp_model$set_spatial_locations(as.matrix(self$locations))
      super$display("  Locations have been set.")
      ## set data
      super$cpp_model$set_data(self$data)
      super$display("  Data have been set.")
    },
    ## centering utils
    init_centering = function(center_init_params) {
      ## init data to center
      data_center <- functional_data(
        domain = self$domain,
        locations = self$locations,
        X = self$data$X
      )
      ## init centering device
      private$R_center <- fCentering_pro(
        data = data_center,
        center = center_init_params,
        VERBOSE = self$VERBOSE
      )
    },
    center = function() {
      ## X center
      private$R_center$fit()
      self$results$calibration$lambda_centering_opt <- private$R_center$results$calibration$lambda_opt
      self$results$X_mean <- private$R_center$results$mean
      self$data$X <- private$R_center$results$centered
      ## Y center (if any)
      if ("Y" %in% names(self$data)) {
        scale_result <- scale(self$data$Y, center = TRUE, scale = FALSE)
        self$results$Y_mean <- as.matrix(colMeans(self$data$Y))
        self$data$Y <- as.matrix(scale_result)
      }
    },
    ## fit utils
    fit = function() {
      ## fitting the statistical model
      super$cpp_model$init()
      self$R0_mat <- super$cpp_model$R0()
      self$Psi_mat <- super$cpp_model$Psi()
      super$cpp_model$solve()
    },
    ## utils
    get_functional_model = function() {
      return(private$functional_model)
    },
    set_functional_model = function(functional_model) {
      private$functional_model <- functional_model
    },
    set_solver = function() {
      super$cpp_model$set_solver(
        SolutionPolicy(private$functional_model$solver$policy),
        private$functional_model$solver$parameters,
        Calibration(private$functional_model$solver$calibrator$name),
        private$functional_model$solver$calibrator$parameters,
        private$functional_model$solver$calibrator$lambda
      )
    },
    set_ncomp = function(n_comp) {
      super$cpp_model$set_ncomp(n_comp)
    },
    ## fPCA cpp_model interfaces
    set_lambda = function(lambda) {
      super$cpp_model$set_lambda(lambda)
    },
    get_scores = function() {
      return(as.matrix(super$cpp_model$scores()))
    },
    get_loadings = function() {
      return(as.matrix(super$cpp_model$loadings()))
    },
    ## fPLS cpp_model interfaces
    fitted = function(h = 0) {
      return(as.matrix(super$cpp_model$fitted(h)))
    },
    reconstructed = function(h = 0) {
      return(as.matrix(super$cpp_model$reconstructed(h)))
    },
    Beta = function(h = 0) {
      return(as.matrix(super$cpp_model$Beta(h)))
    },
    save_fPLS_results = function(MODE = "fPLS-R") {
      switch(MODE, 
        "fPLS-R" = {
          self$results$V <- as.matrix(super$cpp_model$Y_space_directions())
          self$results$W <- as.matrix(super$cpp_model$X_space_directions())
          self$results$C <- as.matrix(super$cpp_model$Y_loadings())
          self$results$R <- as.matrix(super$cpp_model$X_loadings())
          self$results$T <- as.matrix(super$cpp_model$X_latent_scores())
          self$results$U <- as.matrix(super$cpp_model$Y_latent_scores())
        },
        "fPLS-A" = {
          self$results$V <- as.matrix(super$cpp_model$Y_space_directions())
          self$results$W <- as.matrix(super$cpp_model$X_space_directions())
          self$results$S <- as.matrix(super$cpp_model$Y_loadings())
          self$results$R <- as.matrix(super$cpp_model$X_loadings())
          self$results$T <- as.matrix(super$cpp_model$X_latent_scores())
          self$results$U <- as.matrix(super$cpp_model$Y_latent_scores())
        },
        "fPLS-SB" = {
          self$results$V <- as.matrix(super$cpp_model$Y_space_directions())
          self$results$W <- as.matrix(super$cpp_model$X_space_directions())
          self$results$S <- as.matrix(super$cpp_model$Y_loadings())
          self$results$R <- as.matrix(super$cpp_model$X_loadings())
          self$results$T <- as.matrix(super$cpp_model$X_latent_scores())
          self$results$U <- as.matrix(super$cpp_model$Y_latent_scores())
        }
      )

    }
  )
)

### fPCA -----


fPCA_class <- R6::R6Class(
  classname = "fPCA",
  inherit = fdaPDE_Functional_Model,
  public = list(
    ## constructor
    initialize = function(data, fPCA_model, VERBOSE) {
      ## save options
      super$set_verbosity(VERBOSE)
      self$CENTER <- fPCA_model$CENTER
      super$display("\n\nfPCA model")
      ## save fPCA_model
      super$set_functional_model(fPCA_model)
      ## inputs sanity check
      super$display("- Performing inputs sanity check")
      smoother_params <- super$sanity_check_data(data)
      ## centering device initialization
      if (self$CENTER) {
        super$display("- Initializing the centering device")
        super$display("\n\n%% Start centering sub-routine %%")
        super$init_centering(super$get_functional_model()$center)
        super$center()
        super$display("\n\n%% End centering sub-routine %%\n\n")
      }
      ## model initialization
      super$display("- Initializing the model")
      super$init_model(fPCA_model)
    },
    ## utils
    fit = function(calibrator = NULL, solver = NULL, n_pc = NULL) {
      ## get functional model for updating it if necessary
      fPCA_model <- super$get_functional_model()
      ## calibrator update
      if (!is.null(solver)) {
        super$display("- Setting the new solver")
        fPCA_model$solver <- solver
        ## calibrator update
        if (!is.null(calibrator)) {
          super$display("- Setting the new calibrator")
          fPCA_model$solver$calibrator <- parse_calibrator(calibrator)
        }
        ## set changes
        super$set_functional_model(fPCA_model)
        ## set new solver
        super$set_solver()
      } else if (!is.null(calibrator)) { ## calibrator update
        ## update the calibrator of the current solver
        fPCA_model$solver$calibrator <- parse_calibrator(calibrator)
        ## set changes
        super$set_functional_model(fPCA_model)
        ## set new solver
        super$set_solver()
      }
      ## n_pc update
      if (!is.null(n_pc)) {
        fPCA_model$solver$parameters$n_pc <- n_pc
        super$set_functional_model(fPCA_model)
        super$set_ncomp(n_pc)
      }
      ## fit
      super$display("- Fitting the model")
      super$fit()
      ## save results
      super$display("- Saving the results")
      self$results$scores <- super$get_scores()
      self$results$loadings <- super$get_loadings()
      n_stat_units <- nrow(self$data$X)
      if (self$CENTER) {
        self$results$X_hat <- self$results$scores %*% t(self$results$loadings) + rep(1, n_stat_units) %*% t(self$results$X_mean)
      } else {
        self$results$X_hat <- self$results$scores %*% t(self$results$loadings)
      }
    }
  )
)

fPCA_pro <- function(data,
                     fPCA_model = fPCA_init(),
                     VERBOSE = FALSE) {
  model <- fPCA_class$new(
    data = data,
    fPCA_model = fPCA_model,
    VERBOSE = VERBOSE
  )
  return(model)
}

#' @export
fPCA <- function(data,
                 penalty = simple_laplacian_penalty(),
                 center = centering(),
                 solver = sequential(),
                 VERBOSE = FALSE) {
  ## overwriting smoother defaults with the required ones
  ## necessary to guarantee the consistency with the other statistical models
  ## while providing a clear interface to the final user
  fPCA_model <- fPCA_init()
  fPCA_model$penalty <- penalty
  if (is.null(center) || is.logical(center)) {
    fPCA_model$center <- centering
  } else {
    fPCA_model$center <- center
  }
  fPCA_model$CENTER <- ifelse(
    is.null(center),
    TRUE, ## fPCA defaults center to TRUE
    ifelse(is.logical(center),
      center,
      TRUE
    )
  )
  fPCA_model$solver <- solver
  ## return wrapped model
  return(fPCA_pro(
    data = data,
    fPCA_model = fPCA_model,
    VERBOSE = VERBOSE
  ))
}


### fPLS -----


fPLS_class <- R6::R6Class(
  classname = "fPLS",
  inherit = fdaPDE_Functional_Model,
  public = list(
    MODE = NULL,
    ## constructor
    initialize = function(data, fPLS_model, VERBOSE) {
      ## save options
      super$set_verbosity(VERBOSE)
      self$CENTER <- fPLS_model$CENTER
      self$MODE <- fPLS_model$MODE
      super$display("\n\nfPLS model")
      ## save fPLS_model
      super$set_functional_model(fPLS_model)
      ## inputs sanity check
      super$display("- Performing inputs sanity check")
      smoother_params <- super$sanity_check_data(data)
      ## centering device initialization
      if (self$CENTER) {
        super$display("- Initializing the centering device")
        super$display("\n\n%% Start centering sub-routine %%")
        super$init_centering(super$get_functional_model()$center)
        super$center()
        super$display("\n\n%% End centering sub-routine %%\n\n")
      }
      ## model initialization
      super$display("- Initializing the model")
      super$init_model(fPLS_model)
    },
    ## utils
    fit = function(calibrator = NULL, solver = NULL, n_comp = NULL) {
      ## get functional model for updating it ief necessary
      fPLS_model <- super$get_functional_model()
      ## calibrator update
      if (!is.null(solver)) {
        super$display("- Setting the new solver")
        fPLS_model$solver <- solver
        ## check if also the calibrator need to be updated
        if (!is.null(calibrator)) {
          super$display("- Setting the new calibrator")
          fPLS_model$solver$calibrator <- parse_calibrator(calibrator)
        }
        ## set changes
        super$set_functional_model(fPLS_model)
        ## set new solver
        super$set_solver()
      } else if (!is.null(calibrator)) { ## calibrator update
        fPLS_model$solver$calibrator <- parse_calibrator(calibrator)
        ## set changes
        super$set_functional_model(fPLS_model)
        ## set new solver
        super$set_solver()
      }
      ## n_comp update
      if (!is.null(n_comp)) {
        fPLS_model$parameters$n_comp <- n_comp
        super$set_functional_model(fPLS_model)
        super$set_ncomp(n_comp)
      }
      ## fit
      super$display("- Fitting the model")
      super$fit()
      ## save results
      super$display("- Saving the results")
      ## self$results$B_hat <- super$B()
      n_stat_units <- nrow(self$data$X)
      self$results$Y_hat <- list()
      self$results$X_hat <- list()
      for(h in 1:fPLS_model$parameters$n_comp) {
        if(self$CENTER) {
          self$results$Y_hat[[h]] <- super$fitted(h) + rep(1, n_stat_units) %*% t(self$results$Y_mean)
          self$results$X_hat[[h]] <- super$reconstructed(h) + rep(1, n_stat_units) %*% t(self$results$X_mean)
        } else {
          self$results$Y_hat[[h]] <- super$fitted(h)
          self$results$X_hat[[h]] <- super$reconstructed(h)
        }
      }
      if(self$MODE == "fPLS-R") {
        self$results$Beta_hat <- list()
        for(h in 1:fPLS_model$parameters$n_comp) {
          self$results$Beta_hat[[h]] <- super$Beta(h)
        }
      }
      super$save_fPLS_results(self$MODE)
    }
  )
)

fPLS_pro <- function(data,
                     fPLS_model = fPLS_init(),
                     VERBOSE = FALSE) {
  model <- fPLS_class$new(
    data = data,
    fPLS_model = fPLS_model,
    VERBOSE = VERBOSE
  )
  return(model)
}

#' @export
fPLS <- function(data,
                 penalty = simple_laplacian_penalty(),
                 center = centering(),
                 solver = sequential(),
                 smoother = smoothing(),
                 mode = "fPLS-R",
                 VERBOSE = FALSE) {
  ## overwriting smoother defaults with the required ones
  ## necessary to guarantee the consistency with the other statistical models
  ## while providing a clear interface to the final user
  fPLS_model <- fPLS_init()
  fPLS_model$penalty <- penalty
  if (is.null(center) || is.logical(center)) {
    fPLS_model$center <- centering
  } else {
    fPLS_model$center <- center
  }
  fPLS_model$CENTER <- ifelse(
    is.null(center),
    TRUE, ## fPCA defaults center to TRUE
    ifelse(is.logical(center),
      center,
      TRUE
    )
  )
  fPLS_model$MODE <- mode
  fPLS_model$solver <- solver
  fPLS_model$smoother <- smoother
  ## return wrapped model
  return(fPLS_pro(
    data = data,
    fPLS_model = fPLS_model,
    VERBOSE = VERBOSE
  ))
}
