#' @name JoiNMe_r6_classes
#' @aliases JoiNMeFit JoiNMeFit JoiNMeDynPred PredJoiNMeFit SummaryJoiNMeFit SummaryJoiNMeDynPred
#' @title JoiNMe R6 Classes
#'
#' @importFrom R6 R6Class
#'
#' @description
#' R6 classes used to hold JoiNMe fit and prediction results with mutable state.
#'
#' `JoiNMeFit` is the canonical exported fit class. `JoiNMeFit` remains as a
#' compatibility alias.
#' 
#' **JoiNMeFit**: holds fitted model, Stan data, and formulas.
#' Public fields: fit, stan_data, formulaLong, formulaEvent, formulaVCov, config, call, tmax,
#' dataLong, dataEvent.
#' Methods: initialize, cache_get, cache_set.
#'
#' **JoiNMeDynPred**: holds dynamic predictions and metadata.
#' Public fields: predictions, quantiles, draws, data, metadata, call, tmax, n_samples.
#' Methods: initialize, cache_get, cache_set.
#'
#' **SummaryJoiNMeFit**: holds cached summary tables and diagnostics.
#' Public fields: tables, diagnostics, metadata.
#' Methods: initialize.
#'
#' **SummaryJoiNMeDynPred**: holds summary tables for predictions.
#' Public fields: tables, metadata.
#' Methods: initialize.
#'
#' @keywords internal
#' @export JoiNMeFit JoiNMeFit JoiNMeDynPred PredJoiNMeFit SummaryJoiNMeFit SummaryJoiNMeDynPred
NULL


#' @noRd
JoiNMeFit <- R6::R6Class(
  classname = "JoiNMeFit",
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
      class(self) <- unique(c("JoiNMeFit", "JoiNMeFit", class(self)))
    },
    cache_get = function(key) {
      private$cache[[key]]
    },
    cache_set = function(key, value) {
      private$cache[[key]] <- value
      invisible(self)
    },
    cache_clear = function() {
      private$cache <- list()
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
JoiNMeDynPred <- R6::R6Class(
  classname = "JoiNMeDynPred",
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
      class(self) <- unique(c("PredJoiNMeFit", "JoiNMeDynPred", class(self)))
    },
    cache_get = function(key) {
      private$cache[[key]]
    },
    cache_set = function(key, value) {
      private$cache[[key]] <- value
      invisible(self)
    },
    cache_clear = function() {
      private$cache <- list()
      invisible(self)
    }
  ),
  private = list(
    cache = list()
  )
)

#' @noRd
SummaryJoiNMeFit <- R6::R6Class(
  classname = "summary_JoiNMeFit",
  public = list(
    tables = NULL,
    diagnostics = NULL,
    metadata = NULL,
    initialize = function(tables, diagnostics, metadata) {
      # Store summary tables and diagnostics for printing
      self$tables <- tables
      self$diagnostics <- diagnostics
      self$metadata <- metadata
      class(self) <- unique(c("summary_JoiNMeFit", "summary_JoiNMeFit", class(self)))
    }
  )
)

#' @noRd
SummaryJoiNMeDynPred <- R6::R6Class(
  classname = "summary_JoiNMeDynPred",
  public = list(
    tables = NULL,
    metadata = NULL,
    initialize = function(tables, metadata) {
      # Store prediction summary tables
      self$tables <- tables
      self$metadata <- metadata
      class(self) <- unique(c("summary_PredJoiNMeFit", "summary_JoiNMeDynPred", class(self)))
    }
  )
)

#' @export
JoiNMeFit <- JoiNMeFit

#' @export
PredJoiNMeFit <- JoiNMeDynPred
