#' @name joinme_r6_classes
#' @aliases JoinMeFit JoinMeDynPred SummaryJoinMeFit SummaryJoinMeDynPred
#' @title joinme R6 Classes
#'
#' @importFrom R6 R6Class
#'
#' @description
#' R6 classes used to hold joinme fit and prediction results with mutable state.
#' 
#' **JoinMeFit**: holds fitted model, Stan data, and formulas.
#' Public fields: fit, stan_data, formulaLong, formulaEvent, formulaVCov, config, call, tmax,
#' dataLong, dataEvent.
#' Methods: initialize, cache_get, cache_set.
#'
#' **JoinMeDynPred**: holds dynamic predictions and metadata.
#' Public fields: predictions, quantiles, draws, data, metadata, call, tmax, n_samples.
#' Methods: initialize, cache_get, cache_set.
#'
#' **SummaryJoinMeFit**: holds cached summary tables and diagnostics.
#' Public fields: tables, diagnostics, metadata.
#' Methods: initialize.
#'
#' **SummaryJoinMeDynPred**: holds summary tables for predictions.
#' Public fields: tables, metadata.
#' Methods: initialize.
#'
#' @keywords internal
#' @export JoinMeFit JoinMeDynPred SummaryJoinMeFit SummaryJoinMeDynPred
NULL

# File overview:
# - R6 containers for fit, prediction, and summary objects.
# - Lightweight cache support for expensive summaries.

#' @noRd
JoinMeFit <- R6::R6Class(
  classname = "JoinMeFit",
  public = list(
    fit = NULL,
    stan_data = NULL,
    formulaLong = NULL,
    formulaEvent = NULL,
    formulaVCov = NULL,
    config = NULL,
    call = NULL,
    tmax = NULL,
    dataLong = NULL,
    dataEvent = NULL,
    initialize = function(fit, stan_data, formulaLong, formulaEvent, formulaVCov, config, call, tmax, dataLong, dataEvent) {
      # Store core fit objects and metadata for downstream methods
      self$fit <- fit
      self$stan_data <- stan_data
      self$formulaLong <- formulaLong
      self$formulaEvent <- formulaEvent
      self$formulaVCov <- formulaVCov
      self$config <- config
      self$call <- call
      self$tmax <- tmax
      self$dataLong <- dataLong
      self$dataEvent <- dataEvent
    },
    cache_get = function(key) {
      private$cache[[key]]
    },
    cache_set = function(key, value) {
      private$cache[[key]] <- value
      invisible(self)
    }
  ),
  private = list(
    cache = list()
  ),
  active = list(
    family = function() {
      sapply(self$config$family_long, .family_name)
    }
  )
)

#' @noRd
JoinMeDynPred <- R6::R6Class(
  classname = "JoinMeDynPred",
  public = list(
    predictions = NULL,
    quantiles = NULL,
    draws = NULL,
    data = NULL,
    metadata = NULL,
    call = NULL,
    tmax = NULL,
    n_samples = NULL,
    initialize = function(predictions, quantiles, draws, data, metadata, call, tmax, n_samples) {
      # Store prediction outputs and metadata for plotting/summaries
      self$predictions <- predictions
      self$quantiles <- quantiles
      self$draws <- draws
      self$data <- data
      self$metadata <- metadata
      self$call <- call
      self$tmax <- tmax
      self$n_samples <- n_samples
    },
    cache_get = function(key) {
      private$cache[[key]]
    },
    cache_set = function(key, value) {
      private$cache[[key]] <- value
      invisible(self)
    }
  ),
  private = list(
    cache = list()
  )
)

#' @noRd
SummaryJoinMeFit <- R6::R6Class(
  classname = "summary_JoinMeFit",
  public = list(
    tables = NULL,
    diagnostics = NULL,
    metadata = NULL,
    initialize = function(tables, diagnostics, metadata) {
      # Store summary tables and diagnostics for printing
      self$tables <- tables
      self$diagnostics <- diagnostics
      self$metadata <- metadata
    }
  )
)

#' @noRd
SummaryJoinMeDynPred <- R6::R6Class(
  classname = "summary_JoinMeDynPred",
  public = list(
    tables = NULL,
    metadata = NULL,
    initialize = function(tables, metadata) {
      # Store prediction summary tables
      self$tables <- tables
      self$metadata <- metadata
    }
  )
)
