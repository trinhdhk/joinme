# File overview:
# - Lightweight generics that keep JoiNMe independent of rstantools at load time.
# - S3 methods remain compatible with external generics of the same names.

#' Log-likelihood values for JoiNMe objects
#'
#' @param object An object supporting pointwise log-likelihood extraction.
#' @param ... Additional arguments passed to class-specific methods.
#'
#' @return Class-specific log-likelihood output.
#' @export
log_lik <- function(object, ...) {
  UseMethod("log_lik")
}

#' Posterior linear predictor draws
#'
#' @param object An object supporting posterior linear-predictor extraction.
#' @param ... Additional arguments passed to class-specific methods.
#'
#' @return Class-specific posterior linear-predictor output.
#' @export
posterior_linpred <- function(object, ...) {
  UseMethod("posterior_linpred")
}

#' Posterior expected predictor draws
#'
#' @param object An object supporting posterior expected-predictor extraction.
#' @param ... Additional arguments passed to class-specific methods.
#'
#' @return Class-specific posterior expected-predictor output.
#' @export
posterior_epred <- function(object, ...) {
  UseMethod("posterior_epred")
}

#' Posterior predictive draws
#'
#' @param object An object supporting posterior predictive extraction.
#' @param ... Additional arguments passed to class-specific methods.
#'
#' @return Class-specific posterior predictive output.
#' @export
posterior_predict <- function(object, ...) {
  UseMethod("posterior_predict")
}
