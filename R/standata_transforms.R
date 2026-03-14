# Transform Specification Helper
#
# Provides utilities to build and manage transformation specifications
# for the joinme Stan model's association term.
# File overview:
# - Map user transform specs into Stan-friendly bytecode/spline data.
# - Validate functional bytecode and spline shapes before sampling.
.canonicalise_transform_type <- function(type) {
  type <- as.character(type %||% "identity")[1]
  if (identical(type, "ispline_penalized")) {
    return("ispline_penalised")
  }
  if (identical(type, "ispline_expit_penalized")) {
    return("ispline_expit_penalised")
  }
  if (identical(type, "ispline_exp_penalised")) {
    return("ispline_expit_penalised")
  }
  if (identical(type, "ispline_exp_penalized")) {
    return("ispline_expit_penalised")
  }
  type
}

#' @keywords internal
.transform_uses_expit_input <- function(spec_or_type) {
  type <- if (is.list(spec_or_type)) spec_or_type$type %||% "identity" else spec_or_type
  .canonicalise_transform_type(type) %in% c("ispline_expit", "ispline_expit_penalised")
}

#' @keywords internal
.is_ispline_transform_type <- function(spec_or_type) {
  type <- if (is.list(spec_or_type)) spec_or_type$type %||% "identity" else spec_or_type
  .canonicalise_transform_type(type) %in% c(
    "ispline",
    "ispline_penalised",
    "pmonospline",
    "pmono",
    "ispline_expit",
    "ispline_expit_penalised"
  )
}

#' @keywords internal
.transform_input_for_spec <- function(x, spec_or_type) {
  x <- as.numeric(x)
  if (.transform_uses_expit_input(spec_or_type)) {
    return(stats::plogis(x))
  }
  x
}

#' @keywords internal
.validate_expit_domain_values <- function(values, arg_name) {
  values <- as.numeric(values)
  if (!length(values)) {
    return(values)
  }
  if (any(!is.finite(values))) {
    cli::cli_abort(c(
      x = "{.arg {arg_name}} must contain only finite numeric values.",
      i = "For {.val ispline_expit} and {.val ispline_expit_penalised}, provide spline inputs on the expit scale in [0, 1]."
    ))
  }
  if (any(values < 0 | values > 1)) {
    cli::cli_abort(c(
      x = "{.arg {arg_name}} must lie on the expit scale for expit-based spline transforms.",
      i = "Supply values in [0, 1]; {.fn joinme} applies {.fn plogis} only to the raw association feature being transformed."
    ))
  }
  values
}

#' Transform Specification Helper
#'
#' Provides utilities to build and manage transformation specifications
#' for the joinme Stan model's association term.
#'
#' @description
#' Transform can be specified in four modes:
#'   - Mode 0: Identity (no transformation)
#'   - Mode 1: Functional (arbitrary nested functions via parser)
#'   - Mode 2: I-spline basis (monotonic, smooth)
#'   - Mode 2 (penalised): I-spline coefficients estimated with smoothness penalty
#'   - Mode 4: I-spline basis evaluated on `expit(x)` for a bounded, numerically
#'     stable spline input domain
#'   - Mode 3: Piecewise-linear (user-defined interpolation)
#'
#' Each transformation term (CV_total, CS_total, CV_mean, CS_mean, CV_marker,
#' CS_marker, corr) can have its own independent specification, enabling
#' flexibility in model building.
#'
#' @param transform_list List with elements cv_total, cs_total, cv_mean, cs_mean,
#'   cv_marker, cs_marker, corr, each specifying a transformation. See details.
#' @param default_mode Default transformation mode if not specified.
#'
#' @details
#' Each element of transform_list should be a list with:
#'   - type: "identity", "functional", "ispline", "ispline_penalised" (or alias
#'     "ispline_penalized"), "ispline_expit", "ispline_expit_penalised" (aliases
#'     `"ispline_exp_penalised"`, `"ispline_exp_penalized"`,
#'     `"ispline_expit_penalized"`), or "pwlin"
#'   - Additional fields depend on type:
#'     - functional: expr (quosure, formula, quoted expression, or string)
#'     - ispline: knots (vector), coeff (vector), degree (int)
#'     - ispline_penalised: knots or n_knots, degree, lambda, optional x, optional y, optional weights
#'       * if y is supplied, joinme fits the monotone spline to the training pairs (x, y) in R
#'         using the same anchored endpoint convention as the Stan-estimated path
#'         (first coefficient = 0, last coefficient = 1)
#'       * if y is omitted, Stan estimates the monotone spline coefficients directly and lambda controls smoothness
#'     - ispline_expit / ispline_expit_penalised: identical to the I-spline
#'       variants above, except the raw association feature is first mapped
#'       through `plogis(x)` and the spline basis is then evaluated on that
#'       bounded expit-scale input. User-supplied training `x` values and
#'       explicit `knots` for these transform types must therefore already be
#'       specified on the expit scale in `[0, 1]`.
#'     - pwlin: x (vector), y (vector)
#'
#' Defaults and minimal examples:
#'   - identity (default if term omitted): `list(type = "identity")`
#'   - functional: `list(type = "functional", expr = ~ log1p(x))`
#'   - ispline: `list(type = "ispline", knots = c(-1, 0, 1), coeff = c(0, 0.3, 0.8, 1.1, 1.3))`
#'     with `degree = 3` by default when omitted.
#'   - ispline_penalised (plug-in fit):
#'     `list(type = "ispline_penalised", x = seq(-2, 2, length.out = 50), y = exp(seq(-2, 2, length.out = 50)), n_knots = 6, degree = 3, lambda = 1)`
#'   - ispline_penalised (Stan-estimated):
#'     `list(type = "ispline_penalised", x = seq(-2, 2, length.out = 50), n_knots = 6, degree = 3, lambda = 1)`
#'   - ispline_expit:
#'     `list(type = "ispline_expit", knots = c(0.05, 0.5, 0.95), coeff = c(0, 0.25, 0.8, 1.0, 1.1))`
#'   - ispline_expit_penalised (plug-in fit):
#'     `list(type = "ispline_expit_penalised", x = seq(0.02, 0.98, length.out = 50), y = seq(0, 1, length.out = 50), n_knots = 6, degree = 3, lambda = 1)`
#'   - ispline_expit_penalised (Stan-estimated):
#'     `list(type = "ispline_expit_penalised", x = seq(0.02, 0.98, length.out = 50), n_knots = 6, degree = 3, lambda = 1)`
#'   - pwlin:
#'     `list(type = "pwlin", x = c(-2, -1, 0, 1, 2), y = c(0.2, 0.5, 1, 0.5, 0.2))`
#'
#' Default values used internally:
#'   - omitted term -> identity mode (`default_mode = 0`),
#'   - `ispline`: `degree = 3` if omitted,
#'   - `ispline_penalised`: `n_knots = 6`, `degree = 3`, `lambda = 1` if omitted,
#'   - `pwlin` and `functional`: no additional defaults beyond their required fields.
#'
#' For `ispline_expit*` declarations, the internal spline basis lives on the
#' bounded interval $(0, 1)$ after applying `plogis()` to the raw association
#' feature. User-supplied training `x` values and explicit `knots` should be on
#' that same bounded scale so the fitted spline domain matches the runtime basis.
#' This often yields better-conditioned knot placement and smoother optimisation
#' for steep nonlinear transforms.
#'
#' For user-facing declarations, prefer `joinme_tf(...)`, which validates and
#' normalises channel specifications before they are passed here.
#'
#' @return List of standata entries for transformation parameters:
#'   - tf_mode_*
#'   - functional_ops_* and const_data_* (if mode 1)
#'   - knots_* and coeff_* and spline_degree_* (if mode 2, 3, or 4)
#'
#' @section Usage:
#' Use `build_standata_transforms()` to construct the data list entries,
#' and `validate_transforms()` to verify consistency before sampling.
build_standata_transforms <- function(
  transform_list = NULL,
  default_mode = 0  # 0 = identity
) {
  # Output list to be merged into the main standata object
  standata <- list()
  
  term_defaults <- list(
    cv_tot = list(mode_suffix = "cv_tot", short_suffix = "cv"),
    cs_tot = list(mode_suffix = "cs_tot", short_suffix = "cs"),
    cv_mean = list(mode_suffix = "cv_mean", short_suffix = "cv_mean"),
    cs_mean = list(mode_suffix = "cs_mean", short_suffix = "cs_mean"),
    cv_marker = list(mode_suffix = "cv_marker", short_suffix = "cv_marker"),
    cs_marker = list(mode_suffix = "cs_marker", short_suffix = "cs_marker"),
    corr = list(mode_suffix = "corr", short_suffix = "corr")
  )
  
  # Default empty specs for all terms (identity transformation)
  for (term in names(term_defaults)) {
    mode_suffix <- term_defaults[[term]]$mode_suffix
    short_suffix <- term_defaults[[term]]$short_suffix
    standata[[paste0("tf_mode_", mode_suffix)]] <- default_mode
    standata[[paste0("n_functional_ops_", short_suffix)]] <- 0
    standata[[paste0("functional_ops_", short_suffix)]] <- integer()
    standata[[paste0("n_const_", short_suffix)]] <- 0
    standata[[paste0("const_data_", short_suffix)]] <- numeric()
    standata[[paste0("n_knots_", short_suffix)]] <- 0
    standata[[paste0("knots_", short_suffix)]] <- numeric()
    standata[[paste0("n_coeff_", short_suffix)]] <- 0
    standata[[paste0("coeff_", short_suffix)]] <- numeric()
    standata[[paste0("spline_degree_", short_suffix)]] <- 1L
    standata[[paste0("estimate_spline_", short_suffix)]] <- 0L
    standata[[paste0("lambda_spline_", short_suffix)]] <- 0.0
    standata[[paste0("n_free_spline_", short_suffix)]] <- 0L
  }
  
  # Override with user specifications
  # - each term may supply a different mode and parameter payload
  if (!is.null(transform_list)) {
    for (term_name in names(transform_list)) {
      spec <- transform_list[[term_name]]
      term_info <- .resolve_transform_term(term_name)
      mode_suffix <- term_info$mode_suffix
      short_suffix <- term_info$short_suffix
      spec$type <- .canonicalise_transform_type(spec$type)
      
      if (is.null(spec) || spec$type == "identity") {
        standata[[paste0("tf_mode_", mode_suffix)]] <- 0
      } else if (spec$type == "functional") {
        bc <- parse_transform_expr(spec$expr)
        ops <- bc$bytecode %||% bc$opcodes
        standata[[paste0("tf_mode_", mode_suffix)]] <- 1
        standata[[paste0("n_functional_ops_", short_suffix)]] <- length(ops)
        standata[[paste0("functional_ops_", short_suffix)]] <- ops
        standata[[paste0("n_const_", short_suffix)]] <- bc$n_const
        standata[[paste0("const_data_", short_suffix)]] <- bc$const_data
      } else if (.is_ispline_transform_type(spec$type)) {
        if (spec$type %in% c("ispline_penalised", "pmonospline", "pmono", "ispline_expit_penalised")) {
          # Two supported semantics:
          # 1) legacy plug-in mode: user supplies x/y and we fit coefficients in R,
          # 2) Stan-estimated mode: user omits y and Stan estimates spline shape
          #    with a smoothness penalty controlled by lambda.
          if (!is.null(spec$y)) {
            spec <- .make_penalised_ispline_transform(spec)
            spec$estimate_spline <- 0L
            spec$lambda_spline <- 0.0
            spec$n_free_spline <- 0L
          } else {
            spec <- .make_stan_penalised_ispline_transform(spec)
          }
        } else {
          spec$estimate_spline <- 0L
          spec$lambda_spline <- 0.0
          spec$n_free_spline <- 0L
          if (.transform_uses_expit_input(spec$type) && !is.null(spec$knots)) {
            spec$knots <- .validate_expit_domain_values(spec$knots, "knots")
          }
        }
        standata[[paste0("tf_mode_", mode_suffix)]] <- if (.transform_uses_expit_input(spec$type)) 4L else 2L
        standata[[paste0("n_knots_", short_suffix)]] <- length(spec$knots)
        standata[[paste0("knots_", short_suffix)]] <- spec$knots
        standata[[paste0("n_coeff_", short_suffix)]] <- length(spec$coeff)
        standata[[paste0("coeff_", short_suffix)]] <- spec$coeff
        standata[[paste0("spline_degree_", short_suffix)]] <- spec$degree %||% 3L
        standata[[paste0("estimate_spline_", short_suffix)]] <- as.integer(spec$estimate_spline %||% 0L)
        standata[[paste0("lambda_spline_", short_suffix)]] <- as.numeric(spec$lambda_spline %||% 0.0)
        standata[[paste0("n_free_spline_", short_suffix)]] <- as.integer(spec$n_free_spline %||% 0L)
      } else if (spec$type == "pwlin") {
        standata[[paste0("tf_mode_", mode_suffix)]] <- 3
        standata[[paste0("n_knots_", short_suffix)]] <- length(spec$x)
        standata[[paste0("knots_", short_suffix)]] <- spec$x
        standata[[paste0("n_coeff_", short_suffix)]] <- length(spec$y)
        standata[[paste0("coeff_", short_suffix)]] <- spec$y
      } else {
        cli::cli_abort(c(
          x = "Unknown transformation type: {spec$type}.",
          i = "Use one of: identity, functional, ispline, ispline_penalised, ispline_expit, ispline_expit_penalised, pwlin."
        ))
      }
    }
  }
  
  return(standata)
}

##' @keywords internal
.resolve_transform_term <- function(term_name) {
  switch(term_name,
    cv_total = list(mode_suffix = "cv_tot", short_suffix = "cv"),
    cv_tot = list(mode_suffix = "cv_tot", short_suffix = "cv"),
    cs_total = list(mode_suffix = "cs_tot", short_suffix = "cs"),
    cs_tot = list(mode_suffix = "cs_tot", short_suffix = "cs"),
    cv_mean = list(mode_suffix = "cv_mean", short_suffix = "cv_mean"),
    cs_mean = list(mode_suffix = "cs_mean", short_suffix = "cs_mean"),
    cv_marker = list(mode_suffix = "cv_marker", short_suffix = "cv_marker"),
    cs_marker = list(mode_suffix = "cs_marker", short_suffix = "cs_marker"),
    corr = list(mode_suffix = "corr", short_suffix = "corr"),
    cli::cli_abort(c(
      x = "Unknown transform term: {term_name}.",
      i = "Use cv_total, cs_total, cv_mean, cs_mean, cv_marker, cs_marker, or corr."
    ))
  )
}

##' Validate Transformation Specifications
##'
##' Check functional bytecode integrity and spline coefficient shapes before sampling.
##'
##' @param standata Prepared standata list with transformation entries
##' @return TRUE if valid; stops with error message if invalid
##'
##' @keywords internal
validate_transforms <- function(standata) {
  term_map <- list(
    cv = list(mode_suffix = "cv_tot", short_suffix = "cv"),
    cs = list(mode_suffix = "cs_tot", short_suffix = "cs"),
    cv_mean = list(mode_suffix = "cv_mean", short_suffix = "cv_mean"),
    cs_mean = list(mode_suffix = "cs_mean", short_suffix = "cs_mean"),
    cv_marker = list(mode_suffix = "cv_marker", short_suffix = "cv_marker"),
    cs_marker = list(mode_suffix = "cs_marker", short_suffix = "cs_marker"),
    corr = list(mode_suffix = "corr", short_suffix = "corr")
  )
  for (term_name in names(term_map)) {
    mode_suffix <- term_map[[term_name]]$mode_suffix
    short_suffix <- term_map[[term_name]]$short_suffix
    mode <- standata[[paste0("tf_mode_", mode_suffix)]]
    
    if (mode == 1) {
      # Functional bytecode validation
      functional_ops <- standata[[paste0("functional_ops_", short_suffix)]]
      const_data <- standata[[paste0("const_data_", short_suffix)]]
      verify_opcodes(functional_ops, const_data)
    } else if (mode %in% c(2, 3, 4)) {
      # Spline validation
      knots <- standata[[paste0("knots_", short_suffix)]]
      coeff <- standata[[paste0("coeff_", short_suffix)]]
      
      if (length(knots) == 0) {
        cli::cli_abort(c(x = "Spline knots empty for {short_suffix}.", i = "Provide at least two knots."))
      }
      if (length(coeff) == 0) {
        cli::cli_abort(c(x = "Spline coefficients empty for {short_suffix}.", i = "Provide coefficients matching the knots."))
      }
      
      if (any(diff(knots) <= 0)) {
        cli::cli_abort(c(x = "Knots not strictly increasing for {short_suffix}.", i = "Ensure knots are sorted and unique."))
      }
      if (mode == 4 && any(knots < 0 | knots > 1)) {
        cli::cli_abort(c(
          x = "Expit-spline knots must lie in [0, 1] for {short_suffix}.",
          i = "Provide explicit expit-scale knots, or provide expit-scale x values so knots can be derived from them."
        ))
      }
      
      if (mode == 2) {
        # I-spline: coeff length should match knot/degree conventions.
        degree <- standata[[paste0("spline_degree_", short_suffix)]]
        expected_len_a <- length(knots) + degree
        expected_len_b <- length(knots) + degree - 1
        if (!(length(coeff) %in% c(expected_len_a, expected_len_b))) {
          cli::cli_abort(c(
            x = "I-spline coeff length mismatch for {short_suffix}.",
            i = "Expected {expected_len_a} or {expected_len_b} coefficients, got {length(coeff)}."
          ))
        }
      }
    }
  }
  
  return(TRUE)
}

##' Example: Build Transform Specification from Formula String
##'
##' Helper function to construct user-friendly transform specifications
##' for the association term.
##'
##' @export
##'
##' @examples
##' # Example 1: Functional transformation (bytecode-backed)
##' spec_cv <- list(
##'   type = "functional",
##'   expr = ~ (log(sqrt(x + 1/inv_logit(3*x - 3))))^2
##' )
##'
##' # Example 2: Piecewise-linear
##' spec_cs <- list(
##'   type = "pwlin",
##'   x = c(-2, -1, 0, 1, 2),
##'   y = c(0.1, 0.3, 1.0, 0.3, 0.1)  # smooth bump
##' )
##'
##' # Example 3: I-spline (monotonic)
##' # - knots define basis locations (first/last are boundary knots)
##' # - coeff defines the monotone shape directly
##' spec_corr <- list(
##'   type = "ispline",
##'   knots = c(-1, 0, 1),
##'   coeff = c(0, 0.5, 1, 1.2, 1.5),
##'   degree = 3
##' )
##'
##' # Example 4: Penalised monotone I-spline (fit in R)
##' # - x is the input scale of the raw association feature
##' # - y is the desired transformed output at each x
##' # - lambda controls smoothness (higher = smoother)
##' spec_corr_pen <- list(
##'   type = "ispline_penalised",
##'   x = seq(-2, 2, length.out = 50),
##'   y = exp(seq(-2, 2, length.out = 50)),
##'   n_knots = 6,
##'   degree = 3,
##'   lambda = 1.0
##' )
##'
##' # Example 5: Penalised monotone I-spline (Stan-estimated coefficients)
##' # - omit y so Stan learns the monotone shape directly
##' # - if knots are omitted, n_knots = 6 and x quantiles define them
##' spec_corr_pen_stan <- list(
##'   type = "ispline_penalised",
##'   x = seq(-2, 2, length.out = 50),
##'   n_knots = 6,
##'   degree = 3,
##'   lambda = 1.0
##' )
##'
##' # Example 6: Penalised monotone I-spline on expit(x)
##' spec_corr_expit <- list(
##'   type = "ispline_expit_penalised",
##'   x = seq(0.02, 0.98, length.out = 50),
##'   y = seq(0.02, 0.98, length.out = 50)^0.75,
##'   n_knots = 6,
##'   degree = 3,
##'   lambda = 1.0
##' )
##'
##' # Combine and build standata
##' transforms <- joinme_tf(
##'   cv_total = spec_cv,
##'   cs_total = spec_cs,
##'   corr = spec_corr
##' )
##' standata_tf <- build_standata_transforms(transforms)
example_transform_spec <- function() {
  cat("See docstring for examples\n")
}

##' Penalised Monotone I-spline Transform Spec
##'
##' @description
##' Fits a monotone I-spline transformation with a smoothness penalty and
##' returns a transform specification compatible with `build_standata_transforms()`.
##'
##' This helper is the legacy plug-in constructor for `type = "ispline_penalised"`.
##' Its defaults are `n_knots = 6`, `degree = 3`, `lambda = 1.0`,
##' `weights = NULL` (equal weights), and `diff_order = 2`.
##' Returned coefficients follow the same anchored convention as the Stan-
##' estimated path: the first coefficient is fixed at `0`, the last coefficient
##' is fixed at `1`, and interior coefficients are monotone increasing.
##' Example:
##' `penalised_ispline_transform(x = seq(-2, 2, length.out = 50), y = exp(seq(-2, 2, length.out = 50)))`
##'
##' @param x Numeric vector of input values on the raw feature scale to transform.
##' @param y Numeric vector of target transformed values at `x` (same length as `x`).
##' @param knots Optional numeric vector of knots (including boundary knots).
##' @param n_knots Integer. If knots are not provided, number of knots to use
##'   (including boundary knots). Default 6.
##' @param degree Integer spline degree (default 3).
##' @param lambda Non-negative smoothness penalty weight (default 1.0).
##'   Larger values produce smoother fitted transforms; smaller values allow
##'   more local curvature.
##' @param weights Optional non-negative weights (same length as x).
##' @param diff_order Integer difference order for penalty (default 2).
##'
##' @return A list suitable for `build_standata_transforms()`.
##'
##' @export
penalised_ispline_transform <- function(
  x,
  y,
  knots = NULL,
  n_knots = 6,
  degree = 3,
  lambda = 1.0,
  weights = NULL,
  diff_order = 2
) {
  spec <- list(
    type = "ispline_penalised",
    x = x,
    y = y,
    knots = knots,
    n_knots = n_knots,
    degree = degree,
    lambda = lambda,
    weights = weights,
    diff_order = diff_order
  )
  .make_penalised_ispline_transform(spec)
}

##' @keywords internal
##' @rdname penalised_ispline_transform
##' @export
# American alias for `penalised_ispline_transform()`.
penalized_ispline_transform <- function(...) {
  penalised_ispline_transform(...)
}

##' @keywords internal
.make_stan_penalised_ispline_transform <- function(spec) {
  degree <- as.integer(spec$degree %||% 3L)
  if (degree < 1L) {
    cli::cli_abort(c(
      x = "{.arg degree} must be >= 1.",
      i = "Typical choice: degree = 3 (cubic)."
    ))
  }

  lambda <- spec$lambda %||% 1.0
  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda < 0) {
    cli::cli_abort(c(
      x = "{.arg lambda} must be a non-negative numeric scalar.",
      i = "Example: lambda = 1.0."
    ))
  }

  knots_raw <- spec$knots
  knots <- knots_raw
  if (is.null(knots)) {
    if (!is.null(spec$x)) {
      x <- as.numeric(spec$x)
      if (.transform_uses_expit_input(spec)) {
        x <- .validate_expit_domain_values(x, "x")
      }
      n_knots <- as.integer(spec$n_knots %||% 6L)
      if (length(x) < 2L || n_knots < 2L) {
        cli::cli_abort(c(
          x = "Stan-estimated penalised I-spline needs either explicit {.arg knots} or enough {.arg x} values to derive them.",
          i = "Supply knots directly, or pass x with n_knots >= 2."
        ))
      }
      probs <- seq(0, 1, length.out = n_knots)
      knots <- as.numeric(stats::quantile(x, probs = probs, names = FALSE))
    } else {
      cli::cli_abort(c(
        x = "Stan-estimated penalised I-spline requires {.arg knots} when {.arg y} is omitted.",
        i = "Either provide knots directly, or provide x so knots can be derived from quantiles."
      ))
    }
  } else {
    knots <- as.numeric(knots)
    if (.transform_uses_expit_input(spec)) {
      knots <- .validate_expit_domain_values(knots, "knots")
    }
  }

  if (length(knots) < 2L || any(diff(knots) <= 0)) {
    cli::cli_abort(c(
      x = "{.arg knots} must be a strictly increasing vector with at least 2 values.",
      i = "Include boundary knots at the ends."
    ))
  }

  n_coeff <- as.integer(spec$n_coeff %||% (length(knots) + degree - 1L))
  if (n_coeff < 2L) {
    cli::cli_abort(c(
      x = "Penalised I-spline needs at least 2 coefficients.",
      i = "Increase the number of knots or spline degree."
    ))
  }

  list(
    type = if (.transform_uses_expit_input(spec)) "ispline_expit" else "ispline",
    knots = knots,
    raw_knots = if (is.null(knots_raw)) NULL else as.numeric(knots_raw),
    coeff = rep(0, n_coeff),
    degree = degree,
    estimate_spline = 1L,
    lambda_spline = as.numeric(lambda),
    n_free_spline = as.integer(max(0L, n_coeff - 1L))
  )
}

##' @keywords internal
.make_penalised_ispline_transform <- function(spec) {
  # Validate and coerce input pairs for monotone spline fitting
  if (is.null(spec$x) || is.null(spec$y)) {
    cli::cli_abort(c(
      x = "Penalised I-spline requires {.arg x} and {.arg y}.",
      i = "Provide input-output pairs for the monotone transform."
    ))
  }
  x <- as.numeric(spec$x)
  if (.transform_uses_expit_input(spec)) {
    x <- .validate_expit_domain_values(x, "x")
  }
  y <- as.numeric(spec$y)
  if (length(x) != length(y) || length(x) < 2) {
    cli::cli_abort(c(
      x = "{.arg x} and {.arg y} must have the same length >= 2.",
      i = "Check the transform training data."
    ))
  }
  if (!is.numeric(spec$lambda) || length(spec$lambda) != 1 || spec$lambda < 0) {
    cli::cli_abort(c(
      x = "{.arg lambda} must be a non-negative numeric scalar.",
      i = "Example: lambda = 1.0."
    ))
  }
  lambda <- as.numeric(spec$lambda)
  degree <- as.integer(spec$degree %||% 3L)
  if (degree < 1) {
    cli::cli_abort(c(
      x = "{.arg degree} must be >= 1.",
      i = "Typical choice: degree = 3 (cubic)."
    ))
  }

  # Resolve knot locations (explicit or quantile-based)
  knots_raw <- spec$knots
  if (is.null(spec$knots)) {
    n_knots <- as.integer(spec$n_knots %||% 6L)
    if (n_knots < 2) {
      cli::cli_abort(c(
        x = "{.arg n_knots} must be >= 2.",
        i = "Include boundary knots at the ends."
      ))
    }
    probs <- seq(0, 1, length.out = n_knots)
    knots <- as.numeric(stats::quantile(x, probs = probs, names = FALSE))
  } else {
    knots <- as.numeric(spec$knots)
    if (.transform_uses_expit_input(spec)) {
      knots <- .validate_expit_domain_values(knots, "knots")
    }
  }
  if (length(knots) < 2 || any(diff(knots) <= 0)) {
    cli::cli_abort(c(
      x = "{.arg knots} must be a strictly increasing vector with at least 2 values.",
      i = "Include boundary knots at the ends."
    ))
  }

  internal_knots <- if (length(knots) > 2) knots[2:(length(knots) - 1)] else numeric(0)
  boundary_knots <- c(knots[1], knots[length(knots)])

  # splines2 provides the iSpline basis used for penalised fitting
  if (!requireNamespace("splines2", quietly = TRUE)) {
    cli::cli_abort(c(
      x = "Package {.pkg splines2} is required for penalised I-splines.",
      i = "Install splines2 or provide explicit coeff/knots for type = 'ispline'."
    ))
  }

  basis <- splines2::iSpline(x,
    knots = internal_knots,
    degree = degree,
    intercept = TRUE,
    Boundary.knots = boundary_knots
  )
  basis <- as.matrix(basis)

  weights <- spec$weights
  if (!is.null(weights)) {
    weights <- as.numeric(weights)
    if (length(weights) != length(x) || any(weights < 0)) {
      cli::cli_abort(c(
        x = "{.arg weights} must be non-negative and the same length as {.arg x}.",
        i = "Remove weights or provide a valid vector."
      ))
    }
  } else {
    weights <- rep(1, length(x))
  }

  diff_order <- as.integer(spec$diff_order %||% 2L)
  if (diff_order < 1) {
    cli::cli_abort(c(
      x = "{.arg diff_order} must be >= 1.",
      i = "Typical choice: diff_order = 2."
    ))
  }

  n_coef <- ncol(basis)
  if (n_coef < 2L) {
    cli::cli_abort(c(
      x = "Penalised I-spline basis must have at least 2 coefficients.",
      i = "Increase the number of knots or the spline degree."
    ))
  }
  Dmat <- diff(diag(n_coef), differences = diff_order)

  w_sqrt <- sqrt(weights)
  B_w <- basis * w_sqrt
  y_w <- y * w_sqrt

  init_raw <- tryCatch(
    as.numeric(qr.solve(crossprod(B_w), crossprod(B_w, y_w))),
    error = function(e) rep(0, n_coef)
  )

  .anchored_coeff_from_z <- function(z) {
    z <- as.numeric(z)
    z <- z - max(z)
    delta <- exp(z)
    delta <- delta / sum(delta)
    c(0, cumsum(delta))
  }

  init_coeff <- pmax(0, init_raw)
  init_coeff <- init_coeff - init_coeff[1]
  if (!is.finite(init_coeff[n_coef]) || init_coeff[n_coef] <= 0) {
    init_coeff <- seq(0, 1, length.out = n_coef)
  } else {
    init_coeff <- init_coeff / init_coeff[n_coef]
  }
  init_coeff <- cummax(pmin(pmax(init_coeff, 0), 1))
  init_coeff[1] <- 0
  init_coeff[n_coef] <- 1
  delta_init <- diff(init_coeff)
  delta_init <- pmax(delta_init, .Machine$double.eps)
  delta_init <- delta_init / sum(delta_init)
  z_init <- log(delta_init)

  fn <- function(z) {
    b <- .anchored_coeff_from_z(z)
    r <- y_w - B_w %*% b
    pen <- if (lambda > 0 && nrow(Dmat) > 0) Dmat %*% b else 0
    0.5 * sum(r ^ 2) + 0.5 * lambda * sum(pen ^ 2)
  }

  opt <- optim(
    z_init,
    fn,
    method = "BFGS",
    control = list(maxit = 1000)
  )

  if (opt$convergence != 0) {
    cli::cli_warn(c(
      x = "Penalised I-spline optimisation did not fully converge.",
      i = "Consider increasing lambda or adjusting knots."
    ))
  }

  list(
    type = if (.transform_uses_expit_input(spec)) "ispline_expit" else "ispline",
    knots = knots,
    raw_knots = if (is.null(knots_raw)) NULL else as.numeric(knots_raw),
    coeff = as.numeric(.anchored_coeff_from_z(opt$par)),
    degree = degree
  )
}

##' Default Backward-Compatible Transforms
##'
##' If no transforms specified, create identity specs so existing models
##' continue to work unchanged.
##'
##' @keywords internal
default_identity_transforms <- function() {
  list(
    cv_tot = list(type = "identity"),
    cs_tot = list(type = "identity"),
    corr = list(type = "identity")
  )
}
