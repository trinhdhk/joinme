#' The 'JoiNMe' package.
#'
#' @description Joint Mixed-Effects (JoiNMe) Model for Multivariate Longitudinal and Survival Data
#'
#' @details
#' Key entry points are `joinme()` for model fitting, `predict()` for dynamic
#' prediction, and `plot()` for visualisation of predicted trajectories and
#' survival curves.
#'
#' @name JoiNMe-package
#' @aliases JoiNMe
#' @import methods
#' @import Rcpp
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom rstan sampling
#' @importFrom rstantools rstan_config log_lik posterior_linpred posterior_epred posterior_predict
#' @importFrom brms conditional_effects posterior_summary
#' @useDynLib joinme, .registration = TRUE
#'
#' @references
#' Stan Development Team (NA). RStan: the R interface to Stan. R package version 2.36.0.9000. https://mc-stan.org
#'
"_PACKAGE"  # roxygen entry point for package-level docs
