#' The 'joinme' package.
#'
#' @description Joint Mixed-Effects (JoinME) Model for Multivariate Longitudinal and Survival Data
#'
#' @details
#' Key entry points are `joinme()` for model fitting, `predict()` for dynamic
#' prediction, and `plot()` for visualisation of predicted trajectories and
#' survival curves.
#'
#' @name joinme-package
#' @aliases joinme
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#' @importFrom rstantools rstan_config
#' @importFrom RcppParallel RcppParallelLibs
#' @useDynLib joinme, .registration = TRUE
#'
#' @references
#' Stan Development Team (NA). RStan: the R interface to Stan. R package version 2.36.0.9000. https://mc-stan.org
#'
"_PACKAGE"  # roxygen entry point for package-level docs
## @useDynLib joinme, .registration = TRUE
