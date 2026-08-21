#' @name JoiNMe_simulation
#' @title Simulation for JoiNMe
#'
#' @importFrom stats rnorm rexp runif plogis rbeta rbinom rnbinom rpois rt integrate uniroot
#'
#' @description
#' Simulator for the JoiNMe joint model aligned with the v11 Stan semantics.
#'
#' Key features:
#' - Multivariate longitudinal outcomes in long format with irregular times.
#' - Survival hazard depends on CV_total(t) = CV_mean(t) + CV_marker(t).
#' - Survival times drawn by inverse transform sampling using:
#'   * numeric integration (integrate) for cumulative hazard
#'   * numeric root solving (uniroot) for event time.
#'
#' @keywords internal
NULL

#' Softplus function
#' @description Numerically stable implementation of softplus transformation
#' 
#' @details
#' \deqn{\operatorname{softplus}(x) = \log(1+\exp(x))}
#' 
#' @param x Numeric input
#' @export
softplus <- function(x) {
  ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))
}

#' Separate tied event-interval endpoints after simulation
#'
#' @description
#' Numerical root finding can return the lower boundary of an event-time
#' interval when the underlying event occurs closer to zero than the solver's
#' tolerance.  A counting-process likelihood requires a strictly positive
#' interval length.  This helper moves an exactly tied stop time by a
#' scientifically negligible distance that remains representable on the time
#' scale.  Stops below their starts are deliberately left unchanged so that a
#' genuinely invalid interval is still rejected by [joinme_standata()].
#'
#' @param interval_start Numeric vector of event-interval start times.
#' @param interval_stop Numeric vector of event-interval stop times.
#' @param absolute_offset Small positive separation used near time zero.
#'
#' @return Numeric stop times, with exact ties separated from their starts.
#' @keywords internal
#' @noRd
.sim_separate_tied_event_endpoints <- function(
  interval_start,
  interval_stop,
  absolute_offset = 1e-9
) {
  interval_start <- as.numeric(interval_start) # event-interval lower endpoints in the simulation time unit
  interval_stop <- as.numeric(interval_stop) # event-interval upper endpoints before numerical tie correction
  absolute_offset <- as.numeric(absolute_offset) # minimum separation used when the event time is at or near zero
  if (length(interval_start) == 1L && length(interval_stop) != 1L) {
    interval_start <- rep(interval_start, length(interval_stop))
  }
  if (length(interval_start) != length(interval_stop)) {
    cli::cli_abort("Internal event-interval starts and stops must have the same length.")
  }
  if (
    length(absolute_offset) != 1L ||
      !is.finite(absolute_offset) ||
      absolute_offset <= 0
  ) {
    cli::cli_abort("The internal event-interval offset must be one positive finite number.")
  }
  if (any(!is.finite(interval_start)) || any(!is.finite(interval_stop))) {
    cli::cli_abort("Simulated event-interval endpoints must be finite.")
  }

  tied_endpoint <- interval_stop == interval_start # exact numerical ties requiring a strictly positive interval length
  if (!any(tied_endpoint)) return(interval_stop)

  time_magnitude <- pmax(1, abs(interval_start[tied_endpoint])) # local scale used to choose a representable floating-point increment
  representable_offset <- 8 * .Machine$double.eps * time_magnitude # separation exceeding several floating-point units at the current time scale
  endpoint_offset <- pmax(absolute_offset, representable_offset) # final negligible separation on both small and large time scales
  separated_stop <- interval_start[tied_endpoint] + endpoint_offset # corrected stop values for the tied intervals only
  if (any(!is.finite(separated_stop)) || any(separated_stop <= interval_start[tied_endpoint])) {
    cli::cli_abort("Could not represent a stop time strictly greater than its simulated start time.")
  }
  interval_stop[tied_endpoint] <- separated_stop
  interval_stop
}

#' Draw standardised shrinkage latents
#'
#' @param n Number of draws.
#' @param shrinkage Integer family code: 0 Student-t(6), 1 Laplace, 2 Normal.
#' @return Numeric vector of standardised draws.
#' @keywords internal
#' @noRd
.sim_draw_standard_shrinkage <- function(n, shrinkage) {
  n <- as.integer(n)
  if (n <= 0L) return(numeric(0))
  if (shrinkage == 1L) {
    # Inverse-CDF construction for double_exponential(0, 1).
    return(sample(c(-1, 1), n, replace = TRUE) * stats::rexp(n, rate = 1))
  }
  if (shrinkage == 2L) return(stats::rnorm(n))
  stats::rt(n, df = 6)
}

#' Draw marker-weight departures from a declared prior family
#'
#' @description
#' Generates the same centred, unit-scale marker-weight departure used by the
#' fitting programme.  Ordinary Student-t, Normal and Laplace declarations are
#' direct marginal draws.  A regularised horseshoe first draws its local,
#' global and finite-slab scales, then applies the same transformation as the
#' Stan prior module.  Returning the hierarchy alongside the departures makes
#' a simulation study fully auditable without changing the effective weights.
#'
#' @param n Number of marker-specific departures.
#' @param prior A checked `joinme_prior_spec` whose location and ordinary scale
#'   have already been fixed at zero and one.
#' @param n_sets Number of marker-weight sets. The family name `"student_t"`
#'   draws one `2 + Gamma(2, 0.1)` degrees-of-freedom value shared by every
#'   set, exactly as in Stan.
#'
#' @return A named list containing `departure`, `raw`, and any horseshoe scales.
#' @keywords internal
#' @noRd
.sim_draw_marker_weight_prior <- function(n, prior, n_sets = 1L) {
  number_departures <- as.integer(n) # number of marker-by-weight-set deviations required by the simulation
  number_sets <- as.integer(n_sets) # shared or term-specific weight sets governing contiguous marker groups in Stan
  if (length(number_sets) != 1L || is.na(number_sets) || number_sets < 1L ||
      number_departures %% number_sets != 0L) {
    cli::cli_abort("{.arg n_sets} must be positive and divide the number of marker-weight departures exactly.")
  }
  if (number_departures <= 0L) {
    return(list(
      departure = numeric(0),
      raw = numeric(0),
      df = numeric(0),
      df_was_fitted = logical(0),
      local_scale = numeric(0),
      global_scale = numeric(0),
      slab_multiplier = numeric(0)
    ))
  }

  if (identical(prior$family, "student_t")) {
    fitted_df <- 2 + stats::rgamma(1L, shape = 2, rate = 0.1) # shifted-Gamma law for the single tail parameter shared by all marker-weight sets in Stan
    departure <- stats::rt(number_departures, df = fitted_df) # set-major departure vector whose entries all use the same realised degrees of freedom
    return(list(
      departure = departure,
      raw = departure,
      df = fitted_df,
      df_was_fitted = TRUE
    ))
  }
  if (identical(prior$family, "normal")) {
    departure <- stats::rnorm(number_departures) # centred unit-scale Gaussian deviations
    return(list(departure = departure, raw = departure))
  }
  if (identical(prior$family, "laplace")) {
    departure <- sample(c(-1, 1), number_departures, replace = TRUE) *
      stats::rexp(number_departures, rate = 1) # centred unit-scale-parameter Laplace deviations
    return(list(departure = departure, raw = departure))
  }

  raw <- stats::rnorm(number_departures) # standard-Normal coefficient seeds in the regularised horseshoe
  local_scale <- abs(stats::rt(number_departures, df = prior$df)) # coefficient-specific half-Student-t local scales
  global_scale <- abs(stats::rt(1L, df = prior$global_df)) * prior$global_scale # shared half-Student-t global scale for this marker-weight block
  slab_multiplier <- 1 / stats::rgamma(
    1L,
    shape = 0.5 * prior$slab_df,
    rate = 0.5 * prior$slab_df
  ) # inverse-gamma finite-slab variance multiplier used by Stan
  squared_slab <- prior$slab_scale^2 * slab_multiplier # realised finite-slab variance
  regularised_local_scale <- sqrt(
    squared_slab * local_scale^2 /
      (squared_slab + global_scale^2 * local_scale^2)
  ) # local scales after the finite slab has regularised their upper tails
  departure <- raw * global_scale * regularised_local_scale # effective centred marker-weight deviations

  list(
    departure = departure,
    raw = raw,
    local_scale = local_scale,
    global_scale = global_scale,
    slab_multiplier = slab_multiplier
  )
}

#' Draw or fix a population coefficient vector from a simulation declaration
#'
#' @description
#' A `prior_*()` object is interpreted as a population-generating distribution.
#' A numeric declaration is interpreted as the population coefficient itself.
#' This distinction is confined to simulation. A separate fitting-prior object
#' supplies the fitting declaration retained in `truth$recovery`.
#'
#' @param declaration Numeric fixed values or a `joinme_prior_spec`.
#' @param coefficient_names Names and order of the required coefficients.
#' @param context Statistical component named in diagnostic messages.
#'
#' @return A named numeric vector in `coefficient_names` order.
#' @keywords internal
#' @noRd
.sim_resolve_prior_or_fixed <- function(declaration, coefficient_names, context) {
  coefficient_names <- as.character(coefficient_names) # fitted model-matrix order defining the population coefficient vector
  number_coefficients <- length(coefficient_names) # dimension of the requested population regression component
  if (number_coefficients == 0L) return(stats::setNames(numeric(0), coefficient_names))

  if (is.numeric(declaration) && !inherits(declaration, "joinme_prior_spec")) {
    if (!is.null(dim(declaration))) {
      cli::cli_abort("Fixed {.arg {context}} coefficients must be supplied as a numeric vector.")
    }
    supplied_names <- names(declaration) # optional model-matrix labels supplied for exact alignment
    fixed_values <- as.numeric(declaration) # population values held fixed throughout generation
    if (!is.null(supplied_names)) {
      unscoped_coefficient_names <- sub("^(all::|family=[^:]+::)", "", coefficient_names) # readable formula labels accepted for scoped distributional designs
      if (setequal(supplied_names, unscoped_coefficient_names) && !anyDuplicated(unscoped_coefficient_names)) {
        supplied_names <- coefficient_names[match(supplied_names, unscoped_coefficient_names)]
      }
      unknown_names <- setdiff(supplied_names, coefficient_names)
      if (length(unknown_names) > 0L || anyDuplicated(supplied_names)) {
        cli::cli_abort("Fixed {.arg {context}} coefficients have unknown or duplicated names.")
      }
      if (!setequal(supplied_names, coefficient_names)) {
        cli::cli_abort("Named fixed {.arg {context}} coefficients must name every required coefficient exactly once.")
      }
      fixed_values <- fixed_values[match(coefficient_names, supplied_names)]
    } else if (length(fixed_values) == 1L) {
      fixed_values <- rep(fixed_values, number_coefficients)
    } else if (length(fixed_values) != number_coefficients) {
      cli::cli_abort("Fixed {.arg {context}} coefficients must have length one or {number_coefficients}.")
    }
    if (any(!is.finite(fixed_values))) cli::cli_abort("Fixed {.arg {context}} coefficients must be finite.")
    return(stats::setNames(fixed_values, coefficient_names))
  }

  if (!inherits(declaration, "joinme_prior_spec")) {
    cli::cli_abort("{.arg {context}} must be a numeric fixed value or a prior_*() declaration.")
  }
  expand_parameter <- function(value, name) {
    value <- as.numeric(unlist(value, use.names = FALSE))
    if (length(value) == 1L) value <- rep(value, number_coefficients)
    if (length(value) != number_coefficients || any(!is.finite(value))) {
      cli::cli_abort("The {.arg {context}} {name} must have length one or {number_coefficients}.")
    }
    value
  }
  location <- expand_parameter(declaration$mu %||% 0, "location") # population distribution location aligned with coefficient order
  scale <- expand_parameter(declaration$scale %||% 1, "scale") # population distribution scale aligned with coefficient order
  if (any(scale <= 0)) cli::cli_abort("The {.arg {context}} scale must be positive.")
  family <- declaration$family # declared population-generating family
  draws <- switch(
    family,
    normal = stats::rnorm(number_coefficients, location, scale),
    student_t = location + scale * stats::rt(number_coefficients, df = declaration$df),
    laplace = location + scale * sample(c(-1, 1), number_coefficients, replace = TRUE) * stats::rexp(number_coefficients),
    horseshoe = {
      local_scale <- abs(stats::rt(number_coefficients, df = declaration$df)) # coefficient-specific half-Student-t scales
      global_scale <- abs(stats::rt(1L, df = declaration$global_df)) * declaration$global_scale # shared global population scale
      slab_multiplier <- 1 / stats::rgamma(1L, 0.5 * declaration$slab_df, 0.5 * declaration$slab_df) # finite-slab variance multiplier
      squared_slab <- declaration$slab_scale^2 * slab_multiplier # realised slab variance
      regularised_local <- sqrt(squared_slab * local_scale^2 / (squared_slab + global_scale^2 * local_scale^2)) # regularised local scales
      location + scale * stats::rnorm(number_coefficients) * global_scale * regularised_local
    },
    cli::cli_abort("Unsupported population-generating family {.val {family}} for {.arg {context}}.")
  )
  stats::setNames(as.numeric(draws), coefficient_names)
}

#' Resolve intercept and slope population declarations
#'
#' @param raw_component Unnormalised component retained by [jm_truth()].
#' @param checked_component Normalised fitting-prior component.
#' @param coefficient_names Model-matrix coefficient labels.
#' @param context Statistical component used in messages.
#'
#' @return Named population coefficient vector.
#' @keywords internal
#' @noRd
.sim_resolve_regression_component <- function(raw_component, checked_component, coefficient_names, context) {
  raw_component <- raw_component %||% list() # analyst's prior-or-fixed simulation declarations
  if (is.numeric(raw_component) && !inherits(raw_component, "joinme_prior_spec")) {
    return(.sim_resolve_prior_or_fixed(raw_component, coefficient_names, context))
  } # a bare numeric component fixes the complete model-matrix vector
  if (inherits(raw_component, "joinme_prior_spec")) {
    return(.sim_resolve_prior_or_fixed(raw_component, coefficient_names, context))
  }
  intercept_positions <- which(grepl("(^|::)\\(Intercept\\)$", coefficient_names)) # population intercept columns under the fitted design convention
  slope_positions <- setdiff(seq_along(coefficient_names), intercept_positions) # all non-intercept population columns
  output <- stats::setNames(numeric(length(coefficient_names)), coefficient_names) # complete population vector assembled role by role
  if (length(intercept_positions) > 0L) {
    output[intercept_positions] <- .sim_resolve_prior_or_fixed(
      raw_component$intercept %||% checked_component$intercept,
      coefficient_names[intercept_positions], paste0(context, "$intercept")
    )
  }
  if (length(slope_positions) > 0L) {
    output[slope_positions] <- .sim_resolve_prior_or_fixed(
      raw_component$slope %||% checked_component$slope,
      coefficient_names[slope_positions], paste0(context, "$slope")
    )
  }
  output
}

#' Combine global and component population declarations for simulation
#'
#' @param simulation_request Unnormalised declaration retained by [jm_truth()].
#' @param component_name Scientific component name.
#' @param roles Regression roles required by the component.
#'
#' @return A named list containing the most specific declaration for each role.
#' @keywords internal
#' @noRd
.sim_population_component <- function(simulation_request, component_name, roles = c("intercept", "slope")) {
  component <- simulation_request[[component_name]] %||% list() # component-specific population declaration
  if ((is.numeric(component) || inherits(component, "joinme_prior_spec")) && length(roles) > 1L) {
    return(component)
  } # a bare multi-role declaration governs the complete regression vector in model-matrix order
  if (inherits(component, "joinme_prior_spec") || is.numeric(component)) {
    component <- stats::setNames(rep(list(component), length(roles)), roles)
  }
  if (!is.list(component)) component <- list()
  stats::setNames(lapply(roles, function(role_name) {
    component[[role_name]] %||% simulation_request[[role_name]]
  }), roles) # component declaration with the global population fallback filled role by role
}

#' Resolve formulaDist population coefficients, including scoped declarations
#'
#' @param simulation_request Unnormalised declaration retained by [jm_truth()].
#' @param checked_parameter Normalised priors for one distributional parameter.
#' @param parameter_name Canonical formulaDist left-hand-side name.
#' @param coefficient_names Scoped model-matrix column names.
#' @param marker_levels Longitudinal marker labels.
#' @param family_names Family label for each marker.
#'
#' @return A named population coefficient vector in fitted design order.
#' @keywords internal
#' @noRd
.sim_distributional_population <- function(
  simulation_request,
  checked_parameter,
  parameter_name,
  coefficient_names,
  marker_levels,
  family_names
) {
  raw_collection <- simulation_request$distributional %||% list() # all unnormalised formulaDist declarations
  raw_default <- raw_collection[[parameter_name]] # parameter-wide fixed values or generating distributions
  base_component <- .sim_population_component(
    c(simulation_request[c("intercept", "slope")], stats::setNames(list(raw_default), parameter_name)),
    parameter_name
  )
  output <- .sim_resolve_regression_component(
    base_component, checked_parameter, coefficient_names, parameter_name
  ) # parameter-wide coefficients before family or marker refinements

  selector_names <- setdiff(names(raw_collection) %||% character(0), parameter_name)
  for (selector_name in selector_names) {
    selector <- tryCatch(.parse_dist_selector_text(selector_name, allow_marker = TRUE), error = function(error) NULL)
    if (is.null(selector) || !identical(selector$param, parameter_name) || is.null(selector$scope_type)) next
    scope_family <- if (identical(selector$scope_type, "family")) {
      .canonical_family_name(selector$scope_value)
    } else {
      marker_position <- match(selector$scope_value, marker_levels)
      if (is.na(marker_position)) next
      family_names[[marker_position]]
    } # fitted family block selected directly or through its marker label
    positions <- which(grepl(paste0("^family=", scope_family, "::"), coefficient_names))
    if (length(positions) == 0L) next
    scoped_request <- raw_collection[[selector_name]] # population override for this fitted family block
    stripped_names <- sub("^family=[^:]+::", "", coefficient_names[positions])
    if (inherits(scoped_request, "joinme_prior_spec") || is.numeric(scoped_request)) {
      output[positions] <- .sim_resolve_prior_or_fixed(
        scoped_request, stripped_names, selector_name
      ) # bare scoped declaration follows the selected family block's complete model-matrix order
    } else {
      scoped_request <- scoped_request %||% list()
      intercept_positions <- which(grepl("(^|::)\\(Intercept\\)$", stripped_names)) # intercept columns within this selected family block
      slope_positions <- setdiff(seq_along(stripped_names), intercept_positions) # non-intercept columns within this selected family block
      if (!is.null(scoped_request$intercept) && length(intercept_positions) > 0L) {
        output[positions[intercept_positions]] <- .sim_resolve_prior_or_fixed(
          scoped_request$intercept,
          stripped_names[intercept_positions],
          paste0(selector_name, "$intercept")
        )
      }
      if (!is.null(scoped_request$slope) && length(slope_positions) > 0L) {
        output[positions[slope_positions]] <- .sim_resolve_prior_or_fixed(
          scoped_request$slope,
          stripped_names[slope_positions],
          paste0(selector_name, "$slope")
        )
      }
    }
  }
  output
}

#' Declare a time-varying simulation covariate generator
#'
#' @description
#' Helper used inside `covariate_formulas` in `simulate_joinme()` to create
#' subject-specific stepwise covariates over time.
#' Elsewhere, this function does not evaluate anything.
#' Only a specification list is returned.
#'
#' @param fun Random generator function (for example `rnorm`, `rbinom`).
#' @param n_step Integer vector of length 1 or 2.
#'   - length 1: each subject receives exactly `n_step` step periods.
#'   - length 2: each subject draws its number of periods uniformly from
#'     `min(n_step):max(n_step)`.
#'   Ignored when `steps` is supplied.
#' @param steps Optional vector of positive integer time points where the
#'   step changes occur exactly. Defaults to `NULL`.
#'   Values greater than subject-specific stop times are ignored.
#' @param ... Additional arguments forwarded to `fun`.
#'
#' @return Specification fed into `simulate_joinme()`.
#' @export
time_varyring <- function(fun, n_step = NULL, steps = NULL, ...) {
  structure(
    list(fun = fun, n_step = n_step, steps = steps, args = list(...)),
    class = "JoiNMe_time_varyring_spec"
  )
}

#' Weibull baseline hazard factory
#'
#' @param shape Weibull shape (>0).
#' @param scale Weibull scale (>0).
#' @return function h0(t)
#' @export
weibull_h0 <- function(shape = 1.4, scale = 6.0) {
  force(shape)
  force(scale)
  function(t) {
    t <- pmax(t, 0)
    (shape / scale) * (t / scale)^(shape - 1)
  }
}

#' Find a bracketing interval for root solving
#' @keywords internal
#' @noRd
.find_bracket <- function(f, lower = 0, upper = 1, upper_max = 50, expand = 1.7, max_expand = 60L) {
  f_lower <- f(lower)
  if (!is.finite(f_lower) || f_lower > 0) {
    return(NULL)
  }
  f_upper <- f(upper)
  if (!is.finite(f_upper)) {
    return(NULL)
  }
  it <- 0L
  while (f_upper < 0 && upper < upper_max && it < max_expand) {
    upper <- min(upper_max, upper * expand)
    f_upper <- f(upper)
    if (!is.finite(f_upper)) {
      return(NULL)
    }
    it <- it + 1L
  }
  if (f_upper >= 0) list(lower = lower, upper = upper) else NULL
}

#' Simulate joint data (Student-t longitudinal; CV_total survival association)
#'
#' @param n_id Number of subjects.
#' @param D Number of markers.
#' @param n_t Measurements per (id, marker).
#' @param seed RNG seed.
#' @param include_marker_only Include marker-only random effects block.
#' @param Q_idm Dimension of marker-by-id basis.
#' @param R_id Dimension of id-level basis.
#' @param R_mk Dimension of marker-only basis.
#' @param beta Fixed effects (Intercept, time, x1) on original scale.
#' @param gamma_w Hazard covariates (e.g. x1, x2).
#' @param alpha_cv_total Association coefficient for CV_total.
#' @param h0 Baseline hazard function h0(t).
#' @param time_cens Administrative censoring time.
#' @param nu Student-t df.
#' @param sigma_y Student-t scale.
#' @param sigma_fun Optional function to generate observation-level scales.
#'   Signature: function(id, times, dataEvent) -> numeric vector.
#' @param formulaDist Optional distributional regression formulas.
#' @param beta_sigma Optional coefficients for sigma regression.
#' @param beta_nu Optional coefficients for nu regression.
#' @param integration_control list passed to integrate().
#' @param root_control list controlling bracketing.
#'
#' @return list(dataLong, dataEvent, truth, helpers, tmax)
#' @export
simulate_joinme_joint_student_t_cvtotal <- function(
  n_id = 50,
  D = 10,
  n_t = 10,
  seed = .Random.seed[[1]],
  include_marker_only = TRUE,
  Q_idm = 2,
  R_id = 2,
  R_mk = 1,
  beta = c("(Intercept)" = 1.0, "time" = 0.5, "x1" = 0.4),
  gamma_w = c("x1" = 0.3, "x2" = -0.2),
  alpha_cv_total = 0.6,
  h0 = weibull_h0(shape = 1.4, scale = 6.0),
  time_cens = 8.0,
  nu = 5,
  sigma_y = 0.4,
  sigma_fun = NULL,
  formulaDist = NULL,
  beta_sigma = NULL,
  beta_nu = NULL,
  integration_control = list(rel.tol = 1e-6, subdivisions = 2000L, stop.on.error = TRUE),
  root_control = list(t_init = 1.0, t_max = 50.0, expand = 1.7, max_expand = 60L)
) {
  # Workflow: simulate covariates -> random effects -> longitudinal -> survival
  # - matches Stan semantics for CV_total association
  set.seed(seed)
  # Step: configure optional parallel execution for simulation.

  dist_formulas <- .normalize_formula_dist(formulaDist)
  .validate_dist_formula_scopes(dist_formulas, "student_t")

  # Event-level data (exogenous)
  dataEvent <- data.frame(
    id = seq_len(n_id),
    x1 = rnorm(n_id),
    x2 = rnorm(n_id),
    stringsAsFactors = FALSE
  )
  dataEvent$Xcov1 <- 0.0

  # ---- id RE u_id
  tau_u <- rep(0.6, R_id)
  if (R_id >= 2) tau_u[2] <- 0.4
  Corr_u <- diag(R_id)
  if (R_id == 2) Corr_u <- matrix(c(1, 0.2, 0.2, 1), 2, 2)
  Lcorr_u <- t(chol(Corr_u))
  z_u <- matrix(rnorm(n_id * R_id), n_id, R_id)
  L_u <- diag(tau_u, R_id, R_id) %*% Lcorr_u
  u_id <- z_u %*% t(L_u)

  # ---- marker-only v_marker
  if (include_marker_only && R_mk > 0) {
    tau_v <- rep(0.4, R_mk)
    Corr_v <- diag(R_mk)
    Lcorr_v <- t(chol(Corr_v))
    z_v <- matrix(rnorm(D * R_mk), D, R_mk)
    L_v <- diag(tau_v, R_mk, R_mk) %*% Lcorr_v
    v_marker <- z_v %*% t(L_v)
  } else {
    tau_v <- numeric(0)
    Lcorr_v <- matrix(0, 0, 0)
    z_v <- matrix(0, D, 0)
    v_marker <- matrix(0, D, 0)
  }

  # ---- marker-by-id latent z_w (id-specific, iid standard normal)
  z_w_lat <- array(rnorm(n_id * D * Q_idm), dim = c(n_id, D, Q_idm))
  z_w <- z_w_lat

  # ---- id-specific covariance L_i (diagonal)
  alpha_L <- rep(-0.2, Q_idm)
  lambda_L <- rep(0.3, Q_idm)
  z_L <- matrix(rnorm(n_id * Q_idm), nrow = n_id, ncol = Q_idm)

  L_i <- array(0.0, dim = c(n_id, Q_idm, Q_idm))
  for (i in seq_len(n_id)) {
    for (q in seq_len(Q_idm)) {
      lp <- alpha_L[q] + lambda_L[q] * z_L[i, q]
      L_i[i, q, q] <- softplus(lp)
    }
  }

  # w_idm[i,d] = L_i[i] %*% z_w[i,d]
  w_idm <- array(0.0, dim = c(n_id, D, Q_idm))
  for (i in seq_len(n_id)) {
    Li <- matrix(L_i[i, , ], Q_idm, Q_idm)
    for (d in seq_len(D)) w_idm[i, d, ] <- as.numeric(Li %*% z_w[i, d, ])
  }

  # ---- design rows (original time)
  fixed_row <- function(i, t) c("(Intercept)" = 1, time = t, x1 = dataEvent$x1[i])
  id_row <- function(t) if (R_id == 1) c(1) else c(1, t)
  mk_row <- function(i, t) if (!include_marker_only || R_mk <= 0) numeric(0) else c(dataEvent$x1[i])
  idm_row <- function(t) if (Q_idm == 1) c(1) else c(1, t)
  haz_row <- function(i) c(x1 = dataEvent$x1[i], x2 = dataEvent$x2[i])

  # CV components (vectorized over t)
  cv_mean_i <- function(i, t) {
    # fixed_part: (1, t, x1) * beta
    # id_part: (1, t) * u_id[i,]
    x1_val <- dataEvent$x1[i]
    fixed_part <- beta[1] + beta[2] * t + beta[3] * x1_val
    id_part <- if (R_id == 1) u_id[i, 1] else u_id[i, 1] + u_id[i, 2] * t
    as.numeric(fixed_part + id_part)
  }

  cv_marker_i <- function(i, t) {
    # z_idm: (1, t) if Q_idm=2 else (1)
    # z_mk: (x1) if R_mk=1 else (empty)
    # Total CV_marker = mean over markers of [ mk_part + idm_part ]
    x1_val <- dataEvent$x1[i]
    lt <- length(t)
    per_marker_sum <- numeric(lt)

    for (d in seq_len(D)) {
      mk_part <- if (include_marker_only && R_mk > 0) x1_val * v_marker[d, 1] else 0
      idm_part <- if (Q_idm == 1) w_idm[i, d, 1] else w_idm[i, d, 1] + w_idm[i, d, 2] * t
      per_marker_sum <- per_marker_sum + mk_part + idm_part
    }
    per_marker_sum / D
  }
  cv_total_i <- function(i, t) cv_mean_i(i, t) + cv_marker_i(i, t)

  eta_w_i <- function(i) sum(haz_row(i) * gamma_w)

  hazard_i <- function(i, t) {
    h0(t) * exp(eta_w_i(i) + alpha_cv_total * cv_total_i(i, t))
  }

  cumhaz_i <- function(i, t) {
    if (t <= 0) {
      return(0)
    }
    out <- do.call(integrate, c(list(f = function(u) hazard_i(i, u), lower = 0, upper = t), integration_control))
    as.numeric(out$value)
  }

  survival_prob_i <- function(i, t) exp(-cumhaz_i(i, t))

  draw_event_time <- function(i) {
    U <- runif(1)
    target_H <- -log(U)
    f <- function(t) cumhaz_i(i, t) - target_H
    br <- .find_bracket(f,
      lower = 0, upper = root_control$t_init, upper_max = root_control$t_max,
      expand = root_control$expand, max_expand = root_control$max_expand
    )
    if (is.null(br)) {
      return(list(time = time_cens, event = 0L, bracketing_failed = TRUE, U = U))
    }
    T <- uniroot(f, lower = br$lower, upper = br$upper)$root
    if (T > time_cens) {
      list(time = time_cens, event = 0L, bracketing_failed = FALSE, U = U)
    } else {
      list(time = T, event = 1L, bracketing_failed = FALSE, U = U)
    }
  }

  # Draw event times by inverse transform sampling
  death <- lapply(seq_len(n_id), draw_event_time)
  dataEvent$time <- vapply(death, `[[`, numeric(1), "time")
  dataEvent$event <- vapply(death, `[[`, integer(1), "event")

  # ---- longitudinal long format
  marker_levels <- paste0("m", seq_len(D))
  dataLong <- do.call(rbind, lapply(seq_len(n_id), function(i) {
    do.call(rbind, lapply(seq_len(D), function(d) {
      t_max_i <- min(dataEvent$time[i], time_cens)
      t_obs <- sort(runif(n_t, 0, max(1e-8, t_max_i)))
      df <- data.frame(
        id = i,
        marker = factor(marker_levels[d], levels = marker_levels),
        time = t_obs,
        x1 = dataEvent$x1[i],
        x2 = dataEvent$x2[i],
        stringsAsFactors = FALSE
      )
      mu <- numeric(length(t_obs))
      for (k in seq_along(t_obs)) {
        tt <- t_obs[k]
        mu_fixed <- sum(fixed_row(i, tt) * beta)
        mu_id <- sum(id_row(tt) * u_id[i, ])
        mu_mk <- 0
        if (include_marker_only && R_mk > 0) mu_mk <- sum(mk_row(i, tt) * v_marker[d, ])
        mu_idm <- sum(idm_row(tt) * w_idm[i, d, ])
        mu[k] <- mu_fixed + mu_id + mu_mk + mu_idm
      }
      sig <- if (is.null(sigma_fun)) {
        rep(sigma_y, length(mu))
      } else {
        sigma_fun(i, t_obs, dataEvent)
      }
      if (is.null(sigma_fun) && !is.null(dist_formulas$sigma) && !is.null(beta_sigma)) {
        X_sigma <- .build_dist_matrix(dist_formulas$sigma, df, family_by_row = rep("student_t", nrow(df)))$X
        sig <- as.numeric(exp(X_sigma %*% beta_sigma))
      }
      if (length(sig) != length(mu)) {
        stop("sigma_fun must return a vector with length equal to times.")
      }
      nu_vec <- rep(nu, length(mu))
      if (!is.null(dist_formulas$nu) && !is.null(beta_nu)) {
        X_nu <- .build_dist_matrix(dist_formulas$nu, df, family_by_row = rep("student_t", nrow(df)))$X
        nu_vec <- 2 + exp(as.numeric(X_nu %*% beta_nu))
      }
      df$y <- vapply(seq_along(mu), function(k) mu[k] + rt(1, df = nu_vec[k]) * sig[k], numeric(1))
      df
    }))
  }))

  tmax <- max(dataEvent$time)

  # Store latent truth for debugging or benchmark comparisons
  truth <- list(
    beta = beta,
    gamma_w = gamma_w,
    alpha_cv_total = alpha_cv_total,
    nu = nu,
    sigma_y = sigma_y,
    beta_sigma = beta_sigma,
    beta_nu = beta_nu,
    formulaDist = dist_formulas,
    include_marker_only = include_marker_only,
    tau_u = tau_u,
    Lcorr_u = Lcorr_u,
    z_u = z_u,
    u_id = u_id,
    tau_v = tau_v,
    Lcorr_v = Lcorr_v,
    z_v = z_v,
    v_marker = v_marker,
    z_w_lat = z_w_lat,
    z_w = z_w,
    alpha_L = alpha_L,
    lambda_L = lambda_L,
    z_L = z_L,
    L_i = L_i,
    w_idm = w_idm
  )

  # Return helper closures for CV and hazard functions
  helpers <- list(
    cv_mean = cv_mean_i,
    cv_marker = cv_marker_i,
    cv_total = cv_total_i,
    hazard = hazard_i,
    cumhaz = cumhaz_i,
    survival_prob = survival_prob_i
  )

  list(dataLong = dataLong, dataEvent = dataEvent, truth = truth, helpers = helpers, tmax = tmax)
}

#' Resolve marker-specific family declarations for simulation
#'
#' @description
#' Converts character and structured family declarations into the
#' marker-aligned codes, inverse-link instructions, and fixed skew-Laplace
#' quantiles used by `simulate_joinme()`. The same family parser and fixed-`tau`
#' normaliser are used by fitting, which prevents simulation and estimation
#' from assigning different meanings to a family declaration.
#'
#' @param families Character vector or list containing one family declaration
#'   per marker.
#' @param D Positive integer number of longitudinal markers.
#'
#' @return A list containing marker-aligned `family_codes`, `link_names`,
#'   `inv_link_specs`, `use_tau_fixed`, and `tau_fixed`.
#' @keywords internal
#' @noRd
.sim_resolve_family_specs <- function(families, D) {
  if (.is_single_family_spec(families)) {
    families <- rep(list(families), D)
  } else if (length(families) == 1L && !is.list(families)) {
    families <- rep(list(families), D)
  } else if (!is.list(families)) {
    families <- as.list(families)
  }

  if (length(families) != D) {
    cli::cli_abort(
      "families must have length {D} (number of markers), got {length(families)}"
    )
  }

  specifications <- lapply(families, .extract_family_and_link)
  family_codes <- vapply(
    specifications,
    function(specification) specification$family_code,
    integer(1)
  )
  fixed_tau <- .normalise_fixed_tau_by_marker(
    use_tau_fixed = as.integer(vapply(
      specifications,
      function(specification) !is.na(specification$tau_fixed),
      logical(1)
    )),
    tau_fixed = vapply(
      specifications,
      function(specification) {
        if (is.na(specification$tau_fixed)) 0.5 else specification$tau_fixed
      },
      numeric(1)
    ),
    family_codes = family_codes
  )
  link_names <- vapply(
    specifications,
    function(specification) specification$link_name,
    character(1)
  )
  link_names[is.na(link_names)] <- "custom"

  list(
    family_codes = as.integer(family_codes),
    link_names = link_names,
    inv_link_specs = lapply(
      specifications,
      function(specification) specification$inv_link_bc
    ),
    use_tau_fixed = fixed_tau$use_tau_fixed,
    tau_fixed = fixed_tau$tau_fixed
  )
}

#' Simulate joint model data
#' 
#' @description
#' Simulator for `JoiNMe` that mirrors the fitting syntax as closely as
#' possible. The simulator supports:
#' - multivariate outcomes via marker-level families,
#' - multi-structure random effects from `formulaLong` (id, marker, marker-by-id),
#' - formula-based covariate generation through user-provided random generators,
#' - association-driven survival via `assoc`/`formulaAssoc` terms consistent with
#'   the model fit interface.
#'
#' The implementation is model-matrix based end-to-end, so every simulated component
#' is generated from the same formula machinery used during fitting.
#'
#' Population coefficients are declared through `truth`. A `prior_*()` object
#' generates one population coefficient vector and a finite numeric vector
#' fixes that vector exactly. Every realised coefficient is stored in `truth`.
#' Random effects, marker-weight departures, observations and event times remain
#' conditional random realisations.
#'
#' @param formulaLong Longitudinal formula (same role as in `joinme()`).
#'   Grouping terms may use `weighted(group, weights = <column>)` to mirror
#'   fitting syntax. The referenced weight column
#'   should be available in generated covariates (for example via
#'   `covariate_formulas`). In nested marker terms, outer `( ... || marker )`
#'   keeps marker-only and marker-by-id blocks independent, while inner
#'   `( ... || id )` makes the marker-by-id covariance diagonal.
#' @param formulaEvent Event/survival formula (same role as in `joinme()`). Set
#'   this to `NULL` to simulate only the nested longitudinal process. In that
#'   case `dataEvent` is `NULL`, the complete scheduled longitudinal history is
#'   retained, and `assoc` and `formulaAssoc` must be absent.
#' @param formulaVCov Optional covariance-regression specification for the
#'   id-specific marker-by-id covariance factor. A formula is shared by both
#'   covariance components; `list(sd = ~ ..., corr = ~ ...)` supplies
#'   independent observed-covariate regressions with exactly the same syntax as
#'   [joinme()] and [joinme_mix()]. Each formula is evaluated on event-level covariates (one row per subject),
#'   must not include random-effect bars `( ... | ... )`, and must not include
#'   the longitudinal time variable.
#'
#'   The `sd` formula drives subject-specific standard deviations and the
#'   `corr` formula drives Cholesky-correlation-factor rows used to build
#'   `L_i = SD_i * K_i`; their population intercepts, slopes and latent
#'   loadings are declared under `truth$vcov`.
#'   The default `~ 1` gives intercept-only regressions for both components.
#' @param formulaDist Optional distributional regression formulas (same role as in `joinme()`).
#'   Supported LHS parameters are `sigma`, `nu`, `phi`, `alpha` (aliases:
#'   `alpha_skew`, `skew`), `kappa`, and `tau`.
#'
#'   Family-scoped syntax is supported with square brackets:
#'   `param[family=<name>] ~ ...`, for example
#'   `sigma[family=student_t] ~ 1 + time`.
#'
#'   If no bracket is used (e.g. `sigma ~ 1 + time`), the formula applies to all
#'   rows where that parameter exists. If a bracket is used, it only applies to
#'   rows of that family and is estimated separately from other scopes.
#' @param formulaAssoc Optional formula for association terms in hazard, e.g. `~ cv_total + corr + vcov`.
#'   If provided, it overrides `assoc`.
#' @param transforms Optional transform specifications for association terms.
#'   Supports the same structure as `joinme(..., transforms=...)` for
#'   `cv_total`, `cv_mean`, `cv_marker`, `cs_total`, `cs_mean`, `cs_marker`, `corr`, `vcov`.
#'   Simulation applies `cv_total` and `cv_marker` transforms at marker level before
#'   weighted averaging (aligned with fit/predict Stan semantics).
#'   Functional transforms support arithmetic and common nonlinear functions,
#'   including `inv_logit`/`expit`/`sigmoid`, `exp`, `log`, `sqrt`, `power`,
#'   `cbrt`, `softplus`/`log1p_exp`, trigonometric and hyperbolic functions,
#'   the standard normal CDF (`Phi`, `pnorm`), and the standard normal
#'   quantile (`inv_Phi`, `qnorm`, `probit`).
#'
#'   Quick parameterisation reference:
#'   - `list(type = "identity")` or omitted term: identity transform (default).
#'   - `list(type = "functional", expr = ~ log1p(x))`: functional transform;
#'     `expr` is required.
#'   - `list(type = "ispline", knots = c(-1, 0, 1), coeff = c(0, 0.3, 0.8, 1.1, 1.3), degree = 3)`:
#'     direct monotone I-spline; `degree` defaults to `3` if omitted.
#'   - `list(type = "ispline_penalised", x = seq(-2, 2, length.out = 50), y = exp(seq(-2, 2, length.out = 50)), n_knots = 6, degree = 3, lambda = 1)`:
#'     penalised monotone I-spline; defaults are `n_knots = 6`, `degree = 3`,
#'     and `lambda = 1` when omitted.
#'   - `list(type = "ispline_expit", knots = c(0.05, 0.5, 0.95), coeff = c(0, 0.25, 0.8, 1.0, 1.1), degree = 3)`:
#'     monotone I-spline evaluated on `plogis(x)`; explicit knots are supplied
#'     on the expit scale.
#'   - `list(type = "ispline_expit_penalised", x = seq(0.02, 0.98, length.out = 50), y = seq(0.02, 0.98, length.out = 50)^0.8, n_knots = 6, degree = 3, lambda = 1)`:
#'     penalised monotone I-spline on `plogis(x)` in plug-in mode.
#'   - `list(type = "pwlin", x = c(-2, -1, 0, 1, 2), y = c(0.2, 0.5, 1, 0.5, 0.2))`:
#'     piecewise-linear transform; `x` and `y` are required.
#'     Simulation deliberately treats these as fixed interpolation pairs. This
#'     differs from model fitting, where knot ordinates are estimated as an
#'     ordered simplex construction and `y` no longer fixes the curve.
#'
#'   For monotone spline transforms during simulation:
#'   - `type = "ispline"`: provide `knots` and `coeff` directly (plus optional
#'     `degree`).
#'   - `type = "ispline_penalised"` (alias: `"ispline_penalized"`): provide
#'     training pairs `x` and `y`, plus
#'     smoothness penalty `lambda`. The simulator fits the monotone spline in R
#'     before evaluating the transformed association. The fitted plug-in spline
#'     uses anchored endpoint coefficients matching the Stan-estimated path:
#'     increasing splines run from `0` to `1`, while decreasing splines run from
#'     `1` to `0`. If `direction` is omitted, the simulator infers it from the
#'     supplied `(x, y)` pairs.
#'   - `type = "ispline_expit"` / `"ispline_expit_penalised"`: same semantics
#'     as the ordinary I-spline variants above, except the spline basis is built
#'     on `plogis(x)` for a bounded-domain representation. Training `x` values
#'     and explicit `knots` for these transform types are specified on that
#'     bounded expit scale. This is useful
#'     when the raw association feature has long tails or steep nonlinear effects.

#' @param truth A [jm_truth()] declaration for population generation. A numeric
#'   intercept or slope fixes its
#'   population coefficient vector; a `prior_*()` declaration draws that
#'   vector once. This applies to `longitudinal`, `survival`, `baseline`,
#'   `assoc_coef`, `vcov`, `functional`, and named `formulaDist` components. Its
#'   `marker_weights$offset` component supplies the marker-specific
#'   constant contribution and its `marker_weights$family` component generates centred
#'   marker-specific weight departures, because those departures are otherwise
#'   random simulation parameters. That field accepts a family name such as
#'   `"student_t"`, `"normal"`, `"laplace"`, or `"horseshoe"`. The names
#'   `"constant"` and `"none"` use the offset without a random departure. The
#'   family name `"student_t"` draws one value shared by all active sets as
#'   `2 + Gamma(2, 0.1)`, matching fitting. Departure location and ordinary scale remain zero and
#'   one. All declarations and their realised population coefficients are
#'   retained in the realised simulation truth.
#'   Set `marker_weights$family` to `"constant"` or `"none"` when the
#'   declared offset is the complete marker weight. No common location or
#'   marker-specific departure is then drawn. Under a stochastic family, the
#'   effective weight is `offset + marker_weight_mean + departure`, where the
#'   departure is drawn directly from the declared centred unit-scale family.
#'   No additional marker-weight scale is used: the survival association slope
#'   already scales the weighted marker feature.
#' @param shrinkage Integer selecting the random-effect component distribution:
#'   `0` denotes Student-t with 6 degrees of freedom, `1` Laplace, and `2`
#'   Normal. Marker-weight departures no longer use this switch; their family
#'   is declared by `truth = jm_truth(marker_weights = list(family = "normal"))`
#'   or another supported family name.
#' @param n_id Number of subjects.
#' @param families Marker-specific family names or `jm_family()` declarations.
#'   Use `jm_family()` entries to supply custom `link`/`inv_link` expressions
#'   or to fix the skew-Laplace quantile for a marker, for example
#'   `jm_family("skew_laplace", tau = 0.8)`.
#'   A named probit link applies `Phi` as its inverse link. In formula
#'   expressions, `Phi`/`pnorm` are the standard normal CDF and
#'   `inv_Phi`/`qnorm`/`probit` are the standard normal quantile. Simulation
#'   evaluates the same bytecode instructions and fixed-quantile selection as
#'   fitting and prediction.
#' @param marker_levels Optional marker names; defaults to `m1`, `m2`, ...
#' @param times_obs Scheduled observation time grid used for every `(id, marker)`
#'   before optional visit-time jitter and post-event censoring are applied.
#' @param obs_time_noise_sd Optional Gaussian noise SD added independently to each
#'   simulated observation time after the visit schedule is chosen; noisy visit
#'   times are clipped to `[0, time_cens]`.
#' @param censor_longitudinal_after_event Logical; if TRUE (default), simulated
#'   longitudinal observations are truncated at the subject event time. If FALSE,
#'   the longitudinal schedule may continue after the event time up to
#'   `time_cens`.
#' @param seed RNG seed.
#' @param covariate_formulas Named or LHS formulas used to generate event-level covariates,
#'   e.g. `list(x1 ~ rnorm(n_id), x2 ~ rt(n_id, df = 5))`.
#'   Time-varying step covariates are supported via
#'   `time_varyring(fun, n_step, steps, ...)`, for example
#'   `list(x ~ time_varyring(rnorm, c(3, 6), mean = 0, sd = 1))` or
#'   `list(x ~ time_varyring(rnorm, steps = c(1, 3, 5), mean = 0, sd = 1))`.
#'   In this mode, each id receives stepwise periods over `[0, time_cens]`, and
#'   each period value is sampled from `fun`.
#' @param assoc Association components (same names as fit): `cv_total`, `cv_mean`,
#'   `cv_marker`, `cs_total`, `cs_mean`, `cs_marker`, `corr`, `vcov`.
#'
#'   Meaning of each channel:
#'   - `cv_*`: current-value channels from longitudinal trajectories,
#'   - `cs_*`: current-slope channels from finite differences (`eps_cs`),
#'   - `corr`: off-diagonal entries of the subject-specific Cholesky-correlation
#'     factor `K` from marker-by-id random effects.
#'   - `vcov`: the same off-diagonal `K` entries together with the
#'     subject-specific standard deviations.
#'
#'   You can also provide `formulaAssoc = ~ ...` to select channels; when present,
#'   it overrides `assoc`.
#'   Covariance-regression population values belong to `truth$vcov`. Its `sd`
#'   and `corr` components each accept `intercept`, `slope`, and `latent`
#'   declarations. Their element-wise form is
#'   `eta_{im} = alpha_m + x_i^T beta_m + lambda_m z_{im}`, with
#'   `lambda_m >= 0` and `z_{im} ~ Normal(0, 1)`. Diagonal entries apply
#'   `diag_link` to give positive subject-specific standard deviations. Off-diagonal
#'   entries are mapped through `tanh` and then assembled row by row into a valid
#'   subject-specific Cholesky-correlation factor `K_i`. The covariance factor is
#'   reconstructed as `L_i = SD_i * K_i`. Marker-by-id latent seeds are sampled as
#'   iid standard normal values and the final marker-by-id effects are obtained as
#'   `b_id = L_i z_id`.
#'
#'   In the ordinary non-mixture model, \eqn{\lambda_m} is the
#'   conditional standard deviation of the unexplained subject heterogeneity
#'   on covariance-predictor coordinate \eqn{m}, before applying `diag_link` or
#'   `tanh`. It is not itself an entry of `L_i`, a covariance, or a correlation.
#'   With `class_type = "corr"` or `"vcov"`, the selected `z_{im}` has a
#'   class-specific location and scale. Conditional on class \eqn{g}, its
#'   contribution to `eta_{im}` consequently has location
#'   `lambda_m * mix_location[g, m]` and distributional scale
#'   `lambda_m * mix_scale[g, m]`.
#'
#'   The covariance-regression dimension follows the marker-by-id random-effect
#'   dimension `Q_idm`:
#'   - if covariance is full: `M = Q_idm * (Q_idm + 1) / 2` lower-tri entries,
#'   - if covariance is forced diagonal (`||` in nested id-marker term): `M = Q_idm`.
#'
#'   `dist` block syntax:
#'   - `re_params$dist[[param]]` applies to all random-effect terms for that
#'     distributional parameter,
#'   - `re_params$dist[[param]]$terms[[j]]` optionally sets term-specific
#'     controls (same `sd`/`corr` fields as above).
#' @param vcov_diag_link Link applied to covariance-regression scale predictors;
#'   either `"softplus"` or `"exp"`.
#' @param family_params Family-specific constants used when the corresponding
#'   parameter has no distributional regression. In particular, specify the
#'   Beta mean/sample-size model with \code{beta = list(kappa = ...)} and the
#'   skew double exponential model with
#'   \code{skew_double_exponential = list(sigma = ..., tau = ...)}.
#'   \code{kappa} must be positive. \code{tau} must lie in \eqn{(0,1)}, with
#'   \eqn{0.5} giving the symmetric double exponential distribution.
#'   A marker-specific `tau` in `jm_family()` takes precedence over this shared
#'   family constant and over a `tau` distributional regression for that
#'   marker, matching the fitted likelihood.
#' @param time_cens Administrative censoring horizon.
#' @param eps_cs Finite-difference step for slope-type associations (`cs_*`).
#' @param integration_control Control list passed to `integrate()`.
#' @param root_control Root finding control for inverse-CDF sampling.
#' @param quadrature_nodes Optional Gauss-Kronrod node count for simulation-side
#'   survival integration. Allowed values: 7, 15, 31, 41, 51, 61.
#' @param n_workers Number of workers to use when `use_mirai = TRUE`.
#' @param use_mirai Logical; when TRUE and `n_workers > 1`, use `mirai` for
#'   parallel simulation if available. Parallel jobs are seeded deterministically
#'   from the main `seed` for reproducible simulations.
#' @param id_var,marker_var,time_var,y_var,event_time_var,event_var Column names aligned
#'   with `joinme_standata()` defaults.
#' @param .mixture_specification Private latent-progress simulation
#'   specification assembled by [simulate_joinme_mix()]. Users should call
#'   [simulate_joinme_mix()] rather than supplying this argument directly.
#' @param left_truncation_max Optional non-negative delayed-entry bound. When > 0,
#'   each subject receives a sampled entry in
#'   `[0, min(left_truncation_max, stop_time))`.
#' @param truncate_longitudinal_before_entry Logical; when TRUE and delayed entry
#'   is active (`left_truncation_max > 0`), longitudinal rows observed before the
#'   sampled entry time are removed.
#'
#' @return A list containing `dataLong`, `dataEvent`, `truth`, `marker_info`,
#'   `helpers`, and `tmax`. Marker-weight truth
#'   distinguishes base, common-location, departure, and effective values. Baseline-hazard truth is
#'   stored in `truth$baseline_hazard`; when a log-linear representation exists,
#'   its resolved coefficients are also in `truth$stan_fit$bs_gamma_c`.
#'   Marker-specific fixed skew-Laplace quantiles are recorded in
#'   `truth$stan_fit$use_tau_fixed` and `truth$stan_fit$tau_fixed`; realised
#'   row-specific quantiles are recorded in
#'   `truth$distributional$rowwise$tau`.
#'   Hazard-scale association coefficients are stored as `alpha_cv_total`,
#'   `alpha_cs_total`, `alpha_cv_mean`, `alpha_cs_mean`, `alpha_corr`, and
#'   `alpha_vcov`, matching the fitted posterior output names.
#'   `truth$recovery` contains the paired fitting entry-point name and a
#'   directly reusable argument list. Thus
#'   `do.call(joinme, c(sim$truth$recovery$arguments, list(seed = 1)))`
#'   recreates the full fitted specification without retyping formulae, priors,
#'   association controls, or covariance settings.
#'
#' @examples
#' \dontrun{
#' sim <- simulate_joinme(
#'   n_id = 30,
#'   families = list(
#'     jm_family("gaussian"),
#'     jm_family("student_t", inv_link = ~ inv_logit(x / 2)),
#'     jm_family("skew_normal")
#'   ),
#'   quadrature_nodes = 31,
#'   formulaDist = list(
#'     sigma[family=gaussian] ~ 1 + x1,
#'     sigma[family=student_t] ~ 1 + time,
#'     nu[family=student_t] ~ 1,
#'     alpha_skew[family=skew_normal] ~ 1 + x1
#'   )
#' )
#' 
#' sim_split <- simulate_joinme(
#'   n_id = 20,
#'   formulaEvent = survival::Surv(time_start, time_stop, event) ~ x1 + x2,
#'   left_truncation_max = 1.0,
#'   seed = 99
#' )
#'
#' sim_tv <- simulate_joinme(
#'   n_id = 20,
#'   formulaEvent = survival::Surv(time_start, time_stop, event) ~ x + x2,
#'   covariate_formulas = list(
#'     x ~ time_varyring(rnorm, steps = c(1, 3, 5), mean = 0, sd = 1),
#'     x2 ~ rnorm(n_id)
#'   ),
#'   left_truncation_max = 0.5
#' )
#' }
#' @export
simulate_joinme <- function(
  formulaLong = y ~ 1 + time + x1 +
    (1 + time | id) +
    (0 + x1 + (1 + time | id) | marker),
  formulaEvent = survival::Surv(time, event) ~ x1 + x2,
  formulaVCov = ~ 1,
  formulaDist = NULL,
  formulaAssoc = NULL,
  transforms = NULL,
  truth = joinme_truth(),
  n_id = 50,
  families = c("gaussian", "student_t", "binomial"),
  marker_levels = NULL,
  times_obs = seq(0, 5, length.out = 8),
  obs_time_noise_sd = 0,
  censor_longitudinal_after_event = TRUE,
  left_truncation_max = 0,
  truncate_longitudinal_before_entry = TRUE,
  seed = .Random.seed[[1]],
  covariate_formulas = list(
    x1 ~ rnorm(n_id),
    x2 ~ rnorm(n_id)
  ),
  shrinkage = 0L,
  assoc = c("cv_total"),
  vcov_diag_link = "softplus",
  family_params = list(
    gaussian = list(sigma = 1.0),
    student_t = list(sigma = 1.5, nu = 4),
    bernoulli = list(),
    binomial = list(trials = 10),
    poisson = list(),
    negbin2 = list(phi = 2),
    skew_normal = list(sigma = 1.0, alpha = 0),
    double_exponential = list(sigma = 1.0),
    skew_double_exponential = list(sigma = 1.0, tau = 0.5),
    beta = list(kappa = 10),
    cumulative_logit = list(cutpoints = c(-1, 1))
  ),
  time_cens = 8.0,
  eps_cs = 1e-3,
  integration_control = list(rel.tol = 1e-6, subdivisions = 2000L, stop.on.error = TRUE),
  root_control = list(t_init = 1.0, t_max = 50.0, expand = 1.7, max_expand = 60L),
  quadrature_nodes = 15L,
  n_workers = 1L,
  use_mirai = TRUE,
  id_var = "id",
  marker_var = "marker",
  time_var = "time",
  y_var = "y",
  event_time_var = "time",
  event_var = "event",
  .mixture_specification = NULL
) {
  simulation_call <- match.call(expand.dots = TRUE) # supplied arguments used to distinguish the default association from an explicit survival request
  if (!inherits(truth, "joinme_truth")) {
    cli::cli_abort(c(
      x = "{.arg truth} must be created with {.fn jm_truth}.",
      i = "Fixed generating values and their between-simulation distributions belong in that declaration."
    ))
  }
  simulation_truth_request <- unclass(truth) # exact fixed-or-random population declarations resolved after formula dimensions are known
  prior_specification <- attr(truth, "fitting_priors", exact = TRUE) # probability distributions retained solely for the optional recovery fit
  if (!inherits(prior_specification, "joinme_priors")) {
    cli::cli_abort("The {.arg truth} declaration does not contain valid recovery priors.")
  }
  re_params <- simulation_truth_request$re_params # ordinary random-effect covariance truths drawn or fixed once per simulated data set
  baseline_request <- simulation_truth_request$basehaz # baseline-hazard truth represented by one mutually exclusive public form
  h0 <- if (is.function(baseline_request)) baseline_request else NULL # custom hazard function, when explicitly declared
  formulaBasehaz <- if (inherits(baseline_request, "formula")) baseline_request else NULL # log-linear formula hazard, when declared
  baseline_hazard <- if (is.null(h0) && is.null(formulaBasehaz)) baseline_request else NULL # named parametric or spline hazard declaration
  longitudinal_only <- is.null(formulaEvent) # whether generation contains only the nested longitudinal process
  if (longitudinal_only) {
    if ("formulaAssoc" %in% names(simulation_call) && !is.null(formulaAssoc)) {
      cli::cli_abort("{.arg formulaAssoc} requires a survival process in {.fn simulate_joinme}.")
    }
    if ("assoc" %in% names(simulation_call) && length(assoc) > 0L) {
      cli::cli_abort("{.arg assoc} requires a survival process in {.fn simulate_joinme}.")
    }
    formulaEvent <- survival::Surv(time, event) ~ 1 # internal subject scaffold used only while common design matrices are constructed
    formulaAssoc <- NULL # no association channel exists without an event process
    assoc <- character(0) # empty fitted association specification returned in the recovery call
    h0 <- function(t) rep(0, length(t)) # exactly zero internal hazard prevents event generation
    censor_longitudinal_after_event <- FALSE # the complete scheduled longitudinal history is retained
    root_control <- list(
      t_init = time_cens,
      t_max = time_cens,
      expand = 1.7,
      max_expand = 1L
    ) # finite neutral bounds for the common inverse-event-time routine
  }

  # Algorithm overview (statistical simulation workflow):
  # Parse formulas and design structures to mirror the fitted-model
  #         likelihood parameterisation exactly.
  # Simulate subject-level exogenous covariates, including optional
  #         piecewise-constant time-varying processes.
  # Draw fixed/random effects and covariance-regression parameters,
  #         then build latent longitudinal trajectories.
  # Construct the event hazard from baseline + covariate + association
  #         channels and sample event times by inverse-transform sampling.
  # Generate delayed-entry times (left truncation), construct event
  #         interval rows automatically when required by the data-generating
  #         process, and sample longitudinal observations.
  # Draw marker responses from family-specific observation models,
  #         collect truth objects, and return simulation outputs.
  set.seed(seed)
  marker_weight_offsets <- prior_specification$marker_weights$offset # declared constant contribution shared with the fitting interface
  marker_weight_sets_shared <- prior_specification$marker_weights$shared # whether simulated weighted association terms share one marker-weight set
  constant_marker_weights <- identical(prior_specification$marker_weights$family$family, "constant") # constant family suppresses all simulated and fitted marker-weight variation
  if (!is.list(re_params) || is.null(names(re_params))) {
    cli::cli_abort("{.arg re_params} must be a named list.")
  }
  unknown_random_effect_blocks <- setdiff(names(re_params), c("id", "marker", "dist")) # unsupported random-effect simulation controls
  if (length(unknown_random_effect_blocks) > 0L) {
    cli::cli_abort(c(
      x = "Unknown {.arg re_params} component{?s}: {.field {unknown_random_effect_blocks}}.",
      i = "Use {.field id}, {.field marker}, or {.field dist}; covariance-regression coefficients belong to {.arg truth$vcov}."
    ))
  }
  shrinkage <- as.integer(shrinkage)
  if (length(shrinkage) != 1L || is.na(shrinkage) || !(shrinkage %in% 0:2)) {
    cli::cli_abort(c(
      x = "{.arg shrinkage} must be one of 0, 1, or 2.",
      i = "The mappings are 0 = Student-t(6), 1 = Laplace, and 2 = Normal."
    ))
  }
  obs_time_noise_sd <- as.numeric(obs_time_noise_sd %||% 0)
  left_truncation_max <- as.numeric(left_truncation_max %||% 0)
  if (length(left_truncation_max) != 1L || !is.finite(left_truncation_max) || left_truncation_max < 0) {
    cli::cli_abort(c(
      x = "{.arg left_truncation_max} must be one non-negative finite number.",
      i = "Use 0 to disable delayed entry in simulated event data."
    ))
  }
  if (!is.logical(truncate_longitudinal_before_entry) || length(truncate_longitudinal_before_entry) != 1L || is.na(truncate_longitudinal_before_entry)) {
    cli::cli_abort(c(
      x = "{.arg truncate_longitudinal_before_entry} must be TRUE/FALSE.",
      i = "Set TRUE to drop longitudinal rows before delayed-entry times."
    ))
  }
  if (length(obs_time_noise_sd) != 1L || !is.finite(obs_time_noise_sd) || obs_time_noise_sd < 0) {
    cli::cli_abort(c(
      x = "{.arg obs_time_noise_sd} must be one non-negative finite number.",
      i = "Use 0 to disable observation-time jitter."
    ))
  }
  if (!is.logical(censor_longitudinal_after_event) || length(censor_longitudinal_after_event) != 1L || is.na(censor_longitudinal_after_event)) {
    cli::cli_abort(c(
      x = "{.arg censor_longitudinal_after_event} must be TRUE/FALSE.",
      i = "Set TRUE to truncate longitudinal measurements at the event time, or FALSE to keep them through {.arg time_cens}."
    ))
  }
  times_obs <- as.numeric(times_obs %||% numeric(0))
  times_obs <- times_obs[is.finite(times_obs)]
  if (length(times_obs) == 0L) {
    cli::cli_abort(c(
      x = "{.arg times_obs} must contain at least one finite scheduled observation time.",
      i = "The number of longitudinal observations is inferred directly from {.arg times_obs}."
    ))
  }
  times_obs <- sort(times_obs)

  .sim_surv_lhs_arity <- function(formula_event) {
    lhs <- formula_event[[2]]
    if (!is.call(lhs)) return(NA_integer_)
    head_expr <- lhs[[1]]
    is_surv <- identical(as.character(head_expr), "Surv") ||
      (is.call(head_expr) && identical(as.character(head_expr[[1]]), "::") && identical(as.character(head_expr[[3]]), "Surv"))
    if (!is_surv) return(NA_integer_)
    max(0L, length(lhs) - 1L)
  }
  # Base seed for deterministic parallel jobs when mirai is enabled.
  seed_base <- as.integer(seed %||% 1L)
  n_workers <- as.integer(n_workers)
  if (!is.finite(n_workers) || n_workers < 1L) n_workers <- 1L
  use_mirai <- isTRUE(use_mirai) && n_workers > 1L
  if (use_mirai && !requireNamespace("mirai", quietly = TRUE)) {
    cli::cli_warn("Package {.pkg mirai} not installed; falling back to serial simulation.")
    use_mirai <- FALSE
  }
  if (use_mirai) {
    mirai::daemons(n_workers)
    on.exit(mirai::daemons(0L), add = TRUE)
  }

  # Reuse the package-neutral interpreter module.  Binding the functions into
  # this simulation closure also makes the complete evaluator available when
  # the closure is serialised to a mirai worker.
  .sim_eval_bytecode_scalar <- eval_bytecode_scalar
  .sim_eval_bytecode_vector <- eval_bytecode_vector

  # A mirai process begins with a clean R session.  The simulation functions
  # below are closures and therefore carry the quantities defined inside
  # simulate_joinme(), but they do not carry functions found by searching the
  # package namespace.  Supply that small, explicit set of functions to every
  # worker.  This keeps parallel simulation equivalent to serial simulation
  # without asking the worker to load an installed copy of joinme, which may
  # differ from the source tree being examined during package development.
  #
  # The bytecode functions form one self-contained statistical transformation
  # module.  All members of that module are supplied because a functional
  # association may use either the fast vector evaluation or the general stack
  # evaluator.  The model-matrix functions preserve fitted spline bases and
  # factor contrasts while the cumulative hazard is evaluated.  Finally,
  # .find_bracket locates the event-time root itself.
  .sim_mirai_helpers <- list(
    "%||%" = `%||%`,
    ".find_bracket" = .find_bracket,
    ".is_model_matrix_template" = .is_model_matrix_template,
    ".mm" = .mm,
    ".mm_event" = .mm_event,
    ".bytecode_opcodes" = .bytecode_opcodes,
    ".bytecode_normal_ops" = .bytecode_normal_ops,
    ".bytecode_unary_ops" = .bytecode_unary_ops,
    "verify_bytecode" = verify_bytecode,
    ".normalize_bytecode_program" = .normalize_bytecode_program,
    ".eval_canonical_bytecode_vector" = .eval_canonical_bytecode_vector,
    "eval_bytecode_scalar" = eval_bytecode_scalar,
    "eval_bytecode_vector" = eval_bytecode_vector
  )

  #' @param x Vector of inputs to process.
  #' @param fun Function to apply.
  #' @param ... Additional arguments passed to fun.
  #' @return List of results.
  .sim_parallel_lapply <- function(x, fun, ...) {
    # Parallel map with mirai; falls back to serial when mirai is disabled.
    # Each job sets a deterministic seed derived from the main seed.
    if (!use_mirai) return(lapply(x, fun, ...))
    args <- list(...)
    jobs <- lapply(seq_along(x), function(idx) {
      xi <- x[[idx]]
      seed_i <- seed_base + as.integer(idx)
      # Named values supplied through mirai's dots become bindings in the
      # worker's global environment.  This detail is important: a serialised
      # closure searches that environment after its own simulation bindings,
      # whereas values supplied through mirai's local `.args` environment are
      # not visible through the closure's lexical parent.  The fitted formulae
      # and simulation state remain inside fun's closure; only package-level
      # functions need to be added here.
      worker_arguments <- c(
        list(fun = fun, xi = xi, args = args, seed_i = seed_i),
        .sim_mirai_helpers
      )
      do.call(
        mirai::mirai,
        c(
          list(.expr = quote({
            set.seed(seed_i)
            do.call(fun, c(list(xi), args))
          })),
          worker_arguments
        )
      )
    })
    # lapply(jobs, mirai::collect_mirai)
    mirai::collect_mirai(jobs, options = c('.stop', '.progress'))
  }

  .sim_random_signed <- function(n, min_abs = 0.15, max_abs = 0.7) {
    signs <- sample(c(-1, 1), size = n, replace = TRUE)
    magnitudes <- stats::runif(n, min = min_abs, max = max_abs)
    signs * magnitudes
  }

  .sim_resolve_re_block_cfg <- function(cfg, K, label) {
    if (K <= 0) {
      return(list(sd = numeric(0), corr = matrix(0.0, 0, 0), cov = matrix(0.0, 0, 0), Lcorr = matrix(0.0, 0, 0)))
    }

    sd_vec <- cfg$sd
    if (is.null(sd_vec)) {
      sd_vec <- stats::rexp(K, rate = 1) # one population SD per coefficient under the Stan model's exponential(1) prior
    } else {
      sd_vec <- as.numeric(sd_vec)
      if (length(sd_vec) == 1L) sd_vec <- rep(sd_vec, K)
      if (length(sd_vec) != K || any(!is.finite(sd_vec)) || any(sd_vec <= 0)) {
        cli::cli_abort(c(
          x = "Random-effect SD specification is invalid for {label}.",
          i = "Expected {K} finite positive value(s)."
        ))
      }
    }

    corr <- cfg$corr
    if (is.null(corr)) {
      lkj_draw <- .sim_draw_lkj_correlation(
        dimension = K,
        eta = simulation_truth_request$lkj$eta
      ) # one population correlation matrix under the same LKJ law as Stan
      corr <- lkj_draw$corr
      Lcorr <- lkj_draw$Lcorr
    } else {
      corr <- as.matrix(corr)
      if (!all(dim(corr) == c(K, K)) || any(!is.finite(corr))) {
        cli::cli_abort(c(
          x = "Random-effect correlation matrix is invalid for {label}.",
          i = "Expected a finite {K}x{K} matrix."
        ))
      }
      if (!isTRUE(all.equal(corr, t(corr), tolerance = 1e-8))) {
        cli::cli_abort(c(
          x = "Random-effect correlation matrix must be symmetric for {label}.",
          i = "Ensure corr equals its transpose."
        ))
      }
      if (any(abs(diag(corr) - 1) > 1e-8)) {
        cli::cli_abort(c(
          x = "Random-effect correlation matrix must have unit diagonal for {label}.",
          i = "All diagonal entries must equal 1."
        ))
      }
      if (any(abs(corr) > 1 + 1e-8)) {
        cli::cli_abort(c(
          x = "Random-effect correlation entries must lie in [-1, 1] for {label}.",
          i = "Check the supplied correlation matrix."
        ))
      }
      eig <- tryCatch(eigen(corr, symmetric = TRUE, only.values = TRUE)$values, error = function(e) NA_real_)
      if (any(!is.finite(eig)) || min(eig) <= 1e-8) {
        cli::cli_abort(c(
          x = "Random-effect correlation matrix must be positive definite for {label}.",
          i = "Supply a valid correlation matrix."
        ))
      }
      Lcorr <- t(chol(corr)) # lower Cholesky factor of the fixed correlation truth
    }

    cov_mat <- diag(as.numeric(sd_vec), K, K) %*% corr %*% diag(as.numeric(sd_vec), K, K)

    list(sd = sd_vec, corr = corr, cov = cov_mat, Lcorr = Lcorr)
  }

  .sim_inverse_diag_link <- function(x, diag_link) {
    x <- as.numeric(x)
    if (identical(diag_link, "exp")) {
      return(log(x))
    }
    log(expm1(x))
  }

  .sim_eval_covariate_formulas <- function(n_subjects, formulas, base_df) {
    if (is.null(formulas) || length(formulas) == 0) return(base_df)
    if (!is.list(formulas)) formulas <- list(formulas)
    if (is.null(names(formulas))) names(formulas) <- rep("", length(formulas))

    out <- base_df
    tv_specs <- list()

    .sim_draw_time_varyring <- function(fun, n_step = NULL, steps = NULL, fun_args, n_subjects, time_horizon) {
      steps <- steps %||% NULL
      steps_clean <- NULL
      if (!is.null(steps)) {
        steps_num <- as.numeric(steps)
        if (length(steps_num) == 0L || any(!is.finite(steps_num)) || any(steps_num <= 0) || any(abs(steps_num - round(steps_num)) > 1e-8)) {
          cli::cli_abort(c(
            x = "{.fn time_varyring}: {.arg steps} must be a vector of positive integers.",
            i = "Example: {.code steps = c(1, 3, 5)}."
          ))
        }
        steps_clean <- sort(unique(as.integer(round(steps_num))))
      } else {
        n_step <- as.integer(n_step)
        if (!(length(n_step) %in% c(1L, 2L)) || any(!is.finite(n_step)) || any(n_step < 1L)) {
          cli::cli_abort(c(
            x = "{.fn time_varyring} requires {.arg n_step} with length 1 or 2 and positive integers when {.arg steps} is NULL.",
            i = "Use {.code n_step = 4} or {.code n_step = c(3, 7)}."
          ))
        }
      }

      if (is.null(steps_clean)) {
        if (length(n_step) == 1L) {
          n_steps_by_id <- rep.int(n_step, n_subjects)
        } else {
          n_min <- min(n_step)
          n_max <- max(n_step)
          n_steps_by_id <- sample.int(n_max - n_min + 1L, size = n_subjects, replace = TRUE) + n_min - 1L
        }
      } else {
        n_steps_by_id <- rep.int(length(steps_clean) + 1L, n_subjects)
      }

      out_specs <- vector("list", n_subjects)
      for (ii in seq_len(n_subjects)) {
        n_segments <- max(1L, as.integer(n_steps_by_id[ii]))
        if (!is.null(steps_clean)) {
          steps_i <- steps_clean[steps_clean < time_horizon]
          breaks <- sort(unique(c(0, steps_i, time_horizon)))
          if (length(breaks) < 2L) {
            breaks <- c(0, time_horizon)
          }
          n_segments <- length(breaks) - 1L
        } else {
          if (n_segments <= 1L) {
            breaks <- c(0, time_horizon)
          } else {
            inner <- sort(stats::runif(n_segments - 1L, min = 0, max = time_horizon))
            breaks <- c(0, inner, time_horizon)
          }
        }
        n_draw <- length(breaks) - 1L
        vals <- do.call(fun, c(list(n = n_draw), fun_args))
        if (length(vals) == 1L) vals <- rep(vals, n_draw)
        if (length(vals) != n_draw) {
          cli::cli_abort(c(
            x = "{.fn time_varyring}: generator returned length {length(vals)} for {.val {n_draw}} segments.",
            i = "Ensure {.arg fun} returns one value per segment (or a scalar)."
          ))
        }
        out_specs[[ii]] <- list(
          breaks = as.numeric(breaks),
          values = as.numeric(vals),
          steps = if (is.null(steps_clean)) integer(0) else as.integer(steps_clean)
        )
      }

      out_specs
    }

    for (j in seq_along(formulas)) {
      fj <- formulas[[j]]
      if (!inherits(fj, "formula")) fj <- stats::as.formula(fj)

      lhs_name <- names(formulas)[j]
      if (!nzchar(lhs_name)) {
        if (length(fj) < 3) {
          cli::cli_abort(c(
            x = "Each entry in {.arg covariate_formulas} must be named or have an LHS.",
            i = "Example: list(x1 ~ rnorm(n_id), x2 ~ runif(n_id))."
          ))
        }
        lhs_vars <- all.vars(fj[[2]])
        if (length(lhs_vars) != 1) {
          cli::cli_abort(c(
            x = "Covariate generation LHS must contain exactly one variable.",
            i = "Use syntax like x1 ~ rnorm(n_id)."
          ))
        }
        lhs_name <- lhs_vars[1]
      }

      rhs_expr <- if (length(fj) >= 3) fj[[3]] else fj[[2]]
      eval_env <- list2env(c(as.list(out), list(n_id = n_subjects, id = seq_len(n_subjects))), parent = parent.frame())

      is_tv_call <- is.call(rhs_expr) && identical(as.character(rhs_expr[[1]]), "time_varyring")
      if (is_tv_call) {
        rhs_args <- as.list(rhs_expr)[-1]
        rhs_arg_names <- names(rhs_args) %||% rep("", length(rhs_args))

        idx_fun <- which(rhs_arg_names == "fun")[1]
        if (is.na(idx_fun)) idx_fun <- 1L
        fn <- eval(rhs_args[[idx_fun]], envir = eval_env)
        if (!is.function(fn)) {
          cli::cli_abort(c(
            x = "{.fn time_varyring}: {.arg fun} must evaluate to a function.",
            i = "Use generators such as {.fn rnorm}, {.fn rbinom}, or a custom function."
          ))
        }

        idx_steps <- which(rhs_arg_names == "steps")[1]
        idx_n_step <- which(rhs_arg_names == "n_step")[1]
        if (is.na(idx_n_step)) {
          pos_candidates <- setdiff(seq_along(rhs_args), c(idx_fun, idx_steps))
          if (length(pos_candidates) > 0L) idx_n_step <- pos_candidates[1]
        }

        steps <- if (!is.na(idx_steps)) eval(rhs_args[[idx_steps]], envir = eval_env) else NULL
        n_step <- if (!is.na(idx_n_step)) eval(rhs_args[[idx_n_step]], envir = eval_env) else NULL
        if (is.null(steps) && is.null(n_step)) {
          cli::cli_abort(c(
            x = "{.fn time_varyring} requires either {.arg n_step} or {.arg steps}.",
            i = "Example: {.code x ~ time_varyring(rnorm, c(3, 6), mean = 0, sd = 1)} or {.code x ~ time_varyring(rnorm, steps = c(1, 3, 5), mean = 0, sd = 1)}."
          ))
        }

        consumed <- unique(c(idx_fun, idx_n_step, idx_steps))
        consumed <- consumed[is.finite(consumed) & consumed >= 1L]
        extra_idx <- setdiff(seq_along(rhs_args), consumed)
        extra_args <- if (length(extra_idx) > 0L) lapply(rhs_args[extra_idx], eval, envir = eval_env) else list()

        tv_specs[[lhs_name]] <- .sim_draw_time_varyring(
          fun = fn,
          n_step = n_step,
          steps = steps,
          fun_args = extra_args,
          n_subjects = n_subjects,
          time_horizon = time_cens
        )
        sampled <- vapply(tv_specs[[lhs_name]], function(spec) spec$values[1], numeric(1))
      } else {
        sampled <- eval(rhs_expr, envir = eval_env)
      }

      if (length(sampled) == 1) sampled <- rep(sampled, n_subjects)
      if (length(sampled) != n_subjects) {
        cli::cli_abort(c(
          x = "Covariate formula for '{lhs_name}' returned length {length(sampled)}.",
          i = "Expected length {.val {n_subjects}} or scalar."
        ))
      }
      out[[lhs_name]] <- sampled
    }
    attr(out, "time_varyring_specs") <- tv_specs
    out
  }

  .sim_draw_re_block <- function(n_group, K, cfg, label, resolved = NULL) {
    if (K <= 0) return(matrix(0.0, n_group, 0))
    resolved <- resolved %||% .sim_resolve_re_block_cfg(cfg, K, label)
    sd_vec <- resolved$sd
    corr <- resolved$corr
    Sigma <- diag(as.numeric(sd_vec), K, K) %*% corr %*% diag(as.numeric(sd_vec), K, K)
    .sim_draw_mvn(n_group, Sigma)
  }

  .sim_draw_mvn <- function(n_group, Sigma) {
    k_dim <- ncol(Sigma)
    if (k_dim <= 0L || n_group <= 0L) {
      return(matrix(0.0, n_group, k_dim))
    }
    R <- chol(Sigma)
    Z <- matrix(stats::rnorm(n_group * k_dim), n_group, k_dim)
    Z %*% R
  }

  .sim_rhs_matrix <- function(rhs_list, data) {
    if (length(rhs_list) == 0) {
      return(matrix(0.0, nrow(data), 0))
    }
    mats <- lapply(rhs_list, function(rhs) .mm(rhs, data))
    do.call(cbind, mats)
  }

  .sim_get_family_param <- function(fam_name, param_name, fallback) {
    val <- family_params[[fam_name]][[param_name]]
    if (is.null(val)) fallback else val
  }

  .sim_random_assoc_parts <- function(assoc, M_corr, M_vcov) {
    scalar <- assoc_coef_scalar
    scalar_terms <- intersect(names(scalar), assoc)
    if (length(scalar_terms) > 0) {
      for (term in scalar_terms) {
        scalar[[term]] <- if (term %in% c("cv_total", "cv_marker", "cs_total", "cs_marker")) {
          stats::runif(1, min = 0.15, max = 0.75)
        } else {
          .sim_random_signed(1, min_abs = 0.1, max_abs = 0.75)
        }
      }
    }
    list(
      scalar = scalar,
      corr = if (M_corr > 0 && "corr" %in% assoc) .sim_random_signed(M_corr, min_abs = 0.1, max_abs = 0.5) else rep(0.0, M_corr),
      vcov = if (M_vcov > 0 && "vcov" %in% assoc) .sim_random_signed(M_vcov, min_abs = 0.1, max_abs = 0.8) else rep(0.0, M_vcov)
    )
  }

  #' @keywords internal
  #' @param x Numeric vector on link scale.
  #' @param inv_link_bc Parsed inverse-link bytecode.
  #' @return Numeric vector on response scale.
  .sim_apply_inv_link_bc <- function(x, inv_link_bc) {
    # Apply a parsed inverse-link bytecode using the shared interpreter.
    .sim_eval_bytecode_vector(
      x = x,
      bytecode = inv_link_bc$bytecode,
      const_data = inv_link_bc$const_data %||% numeric(0)
    )
  }

  .sim_sample_ordinal <- function(eta, cutpoints) {
    cp <- sort(as.numeric(cutpoints))
    cdf_vals <- stats::plogis(cp - eta)
    probs <- c(cdf_vals[1], diff(cdf_vals), 1 - cdf_vals[length(cdf_vals)])
    probs <- pmax(probs, 1e-12)
    probs <- probs / sum(probs)
    sample.int(length(probs), size = 1, prob = probs)
  }

  .sim_assoc_from_formula <- function(f) {
    if (is.null(f)) return(NULL)
    vars <- all.vars(f)
    unique(vars[vars %in% c("cv_total", "cv_mean", "cv_marker", "cs_total", "cs_mean", "cs_marker", "corr", "vcov")])
  }

  .sim_normalize_transform_list <- function(transform_list) {
    unclass(.normalise_joinme_tf_input(transform_list, validate = FALSE))
  }

  .sim_validate_fixed_transform_shifts <- function(transforms) {
    if (is.null(transforms) || !length(transforms)) {
      return(transforms)
    }
    has_fit_only_iota <- vapply(transforms, function(spec) {
      isTRUE((spec$n_iota_intercept %||% 0L) > 0L) ||
        isTRUE((spec$n_iota_slope %||% 0L) > 0L) ||
        isTRUE(spec$estimate_iota_intercept) ||
        isTRUE(spec$estimate_iota_slope)
    }, logical(1))
    if (any(has_fit_only_iota)) {
      bad_terms <- names(transforms)[has_fit_only_iota]
      cli::cli_abort(c(
        x = "Fit-only functional transform affine shifts are not available in {.fn simulate_joinme}.",
        i = "Remove {.arg intercept}/{.arg slope} from {.fn joinme_tf} for simulation, or encode fixed shifts directly in the transform expression for: {.val {paste(bad_terms, collapse = ', ')}}."
      ))
    }
    transforms
  }

  .sim_canonicalise_transform_type <- function(type) {
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

  .sim_is_ispline_transform_type <- function(spec_or_type) {
    type <- if (is.list(spec_or_type)) spec_or_type$type %||% "identity" else spec_or_type
    .sim_canonicalise_transform_type(type) %in% c(
      "ispline",
      "ispline_penalised",
      "pmonospline",
      "pmono",
      "ispline_expit",
      "ispline_expit_penalised"
    )
  }

  .sim_transform_uses_expit_input <- function(spec_or_type) {
    type <- if (is.list(spec_or_type)) spec_or_type$type %||% "identity" else spec_or_type
    .sim_canonicalise_transform_type(type) %in% c("ispline_expit", "ispline_expit_penalised")
  }

  .sim_validate_expit_domain_values <- function(values, arg_name) {
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
        i = "Supply values in [0, 1]; {.fn simulate_joinme} applies {.fn plogis} only to the raw association feature being transformed."
      ))
    }
    values
  }

  .sim_transform_input_for_spec <- function(x, spec_or_type) {
    x <- as.numeric(x)
    if (.sim_transform_uses_expit_input(spec_or_type)) {
      return(stats::plogis(x))
    }
    x
  }

  .sim_resolve_monotone_direction <- function(direction = NULL, default = "increasing") {
    direction <- direction %||% default
    if (is.numeric(direction) && length(direction) == 1L && is.finite(direction)) {
      if (direction > 0) return(1L)
      if (direction < 0) return(-1L)
    }
    direction_chr <- tolower(trimws(as.character(direction)[1]))
    if (direction_chr %in% c("increasing", "increase", "inc", "+1", "1")) return(1L)
    if (direction_chr %in% c("decreasing", "decrease", "dec", "-1")) return(-1L)
    cli::cli_abort(c(
      x = "{.arg direction} must be either {.val increasing} or {.val decreasing}.",
      i = "Numeric aliases {.val 1} and {.val -1} are also accepted."
    ))
  }

  .sim_infer_monotone_direction <- function(x, y) {
    ord <- order(x)
    x <- as.numeric(x)[ord]
    y <- as.numeric(y)[ord]
    keep <- is.finite(x) & is.finite(y)
    x <- x[keep]
    y <- y[keep]
    if (length(y) < 2L) return(1L)
    if (y[length(y)] < y[1]) return(-1L)
    1L
  }

  .sim_make_penalised_ispline_transform <- function(spec) {
    if (is.null(spec$x) || is.null(spec$y)) {
      cli::cli_abort(c(
        x = "Penalised I-spline requires {.arg x} and {.arg y}.",
        i = "Provide input-output pairs for the monotone transform."
      ))
    }

    x <- as.numeric(spec$x)
    if (.sim_transform_uses_expit_input(spec)) {
      x <- .sim_validate_expit_domain_values(x, "x")
    }
    y <- as.numeric(spec$y)
    if (length(x) != length(y) || length(x) < 2L) {
      cli::cli_abort(c(
        x = "{.arg x} and {.arg y} must have the same length >= 2.",
        i = "Check the transform training data."
      ))
    }

    spline_direction <- .sim_resolve_monotone_direction(spec$direction, default = .sim_infer_monotone_direction(x, y))

    lambda <- as.numeric(spec$lambda %||% 1.0)
    if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda < 0) {
      cli::cli_abort(c(
        x = "{.arg lambda} must be a non-negative numeric scalar.",
        i = "Example: lambda = 1.0."
      ))
    }

    degree <- as.integer(spec$degree %||% 3L)
    if (degree < 1L) {
      cli::cli_abort(c(
        x = "{.arg degree} must be >= 1.",
        i = "Typical choice: degree = 3 (cubic)."
      ))
    }

    knots_raw <- spec$knots
    if (is.null(knots_raw)) {
      n_knots <- as.integer(spec$n_knots %||% 6L)
      if (n_knots < 2L) {
        cli::cli_abort(c(
          x = "{.arg n_knots} must be >= 2.",
          i = "Include boundary knots at the ends."
        ))
      }
      probs <- seq(0, 1, length.out = n_knots)
      knots <- as.numeric(stats::quantile(x, probs = probs, names = FALSE))
    } else {
      knots <- as.numeric(knots_raw)
      if (.sim_transform_uses_expit_input(spec)) {
        knots <- .sim_validate_expit_domain_values(knots, "knots")
      }
    }

    if (length(knots) < 2L || any(diff(knots) <= 0)) {
      cli::cli_abort(c(
        x = "{.arg knots} must be a strictly increasing vector with at least 2 values.",
        i = "Include boundary knots at the ends."
      ))
    }

    if (!requireNamespace("splines2", quietly = TRUE)) {
      cli::cli_abort(c(
        x = "Package {.pkg splines2} is required for penalised I-splines.",
        i = "Install splines2 or provide explicit coeff/knots for type = 'ispline'."
      ))
    }

    internal_knots <- if (length(knots) > 2L) knots[2:(length(knots) - 1L)] else numeric(0)
    boundary_knots <- c(knots[1], knots[length(knots)])
    basis <- splines2::iSpline(
      x,
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
    if (diff_order < 1L) {
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

    .sim_anchored_coeff_from_z <- function(z) {
      z <- as.numeric(z)
      z <- z - max(z)
      delta <- exp(z)
      delta <- delta / sum(delta)
      c(0, cumsum(delta))
    }

    init_coeff <- pmax(0, init_raw)
    init_coeff <- init_coeff - init_coeff[1]
    span <- max(init_coeff, na.rm = TRUE)
    if (!is.finite(span) || span <= 0) {
      init_coeff <- seq(0, 1, length.out = n_coef)
    } else {
      init_coeff <- init_coeff / span
      init_coeff <- pmin(pmax(init_coeff, 0), 1)
      init_coeff <- cummax(init_coeff)
      init_coeff[1] <- 0
      init_coeff[n_coef] <- 1
    }
    delta_init <- diff(init_coeff)
    delta_init <- pmax(delta_init, .Machine$double.eps)
    delta_init <- delta_init / sum(delta_init)
    z_init <- log(delta_init)

    fn <- function(z) {
      b <- .sim_anchored_coeff_from_z(z)
      fit_vals <- as.numeric(B_w %*% b)
      if (spline_direction < 0L) {
        fit_vals <- w_sqrt * (sum(b) - as.numeric(basis %*% b))
      }
      r <- y_w - fit_vals
      pen <- if (lambda > 0 && nrow(Dmat) > 0) Dmat %*% b else 0
      0.5 * sum(r^2) + 0.5 * lambda * sum(pen^2)
    }

    opt <- optim(z_init, fn, method = "BFGS", control = list(maxit = 1000))
    if (opt$convergence != 0) {
      cli::cli_warn(c(
        x = "Penalised I-spline optimisation did not fully converge.",
        i = "Consider increasing lambda or adjusting knots."
      ))
    }

    list(
      type = if (.sim_transform_uses_expit_input(spec)) "ispline_expit" else "ispline",
      knots = knots,
      raw_knots = if (is.null(knots_raw)) NULL else as.numeric(knots_raw),
      coeff = as.numeric(.sim_anchored_coeff_from_z(opt$par)),
      degree = degree,
      direction = if (spline_direction < 0L) "decreasing" else "increasing",
      spline_direction = spline_direction
    )
  }

  .sim_make_assoc_transform <- function(spec, term_name) {
    if (is.null(spec) || is.null(spec$type) || identical(spec$type, "identity")) {
      return(function(x) x)
    }

    tf_type <- .sim_canonicalise_transform_type(spec$type)

    if (identical(tf_type, "functional")) {
      bc <- parse_transform_expr(spec$expr)
      return(function(x) {
        .sim_eval_bytecode_vector(
          x = x,
          bytecode = bc$bytecode,
          const_data = bc$const_data %||% numeric(0)
        )
      })
    }

    if (.sim_is_ispline_transform_type(tf_type)) {
      if (tf_type %in% c("ispline_penalised", "pmonospline", "pmono", "ispline_expit_penalised")) {
        if (is.null(spec$y)) {
          cli::cli_abort(c(
            x = "Simulation requires both {.arg x} and {.arg y} for penalised spline transforms.",
            i = "Use plug-in mode in {.fn simulate_joinme} by supplying x/y pairs, or provide an explicit {.val ispline} transform instead.",
            i = "Stan-estimated penalised splines without y are supported in {.fn JoiNMe}, not in simulation."
          ))
        }
        spec <- .sim_make_penalised_ispline_transform(spec)
      }
      if (!requireNamespace("splines2", quietly = TRUE)) {
        cli::cli_abort(c(
          x = "Transform type {.val {tf_type}} for {.val {term_name}} requires {.pkg splines2}.",
          i = "Install {.pkg splines2} or use identity/functional/pwlin in simulation."
        ))
      }
      knots <- as.numeric(spec$knots)
      if (.sim_transform_uses_expit_input(spec)) {
        knots <- .sim_validate_expit_domain_values(knots, "knots")
      }
      coeff <- as.numeric(spec$coeff)
      degree <- as.integer(spec$degree %||% 3L)
      return(function(x) {
        x <- .sim_transform_input_for_spec(x, spec)
        if (length(knots) < 2) {
          cli::cli_abort(c(
            x = "I-spline transform for {.val {term_name}} requires at least 2 knots.",
            i = "Provide boundary knots as the first and last entries."
          ))
        }
        boundary <- c(knots[1], knots[length(knots)])
        x <- pmin(pmax(x, boundary[1]), boundary[2])
        internal_knots <- if (length(knots) > 2) knots[2:(length(knots) - 1)] else numeric(0)
        basis <- splines2::iSpline(x,
          knots = internal_knots,
          degree = degree,
          intercept = TRUE,
          Boundary.knots = boundary
        )
        coef_len <- ncol(basis)
        if (length(coeff) < coef_len) {
          coeff_use <- c(coeff, rep(0, coef_len - length(coeff)))
        } else {
          coeff_use <- coeff[seq_len(coef_len)]
        }
        vals <- as.numeric(basis %*% coeff_use)
        if (.sim_resolve_monotone_direction(spec$direction %||% spec$spline_direction %||% 1L, default = 1L) < 0L) {
          vals <- sum(coeff_use) - vals
        }
        vals
      })
    }

    if (identical(tf_type, "pwlin")) {
      xk <- as.numeric(spec$x)
      yk <- as.numeric(spec$y)
      if (length(xk) < 2 || length(yk) != length(xk)) {
        cli::cli_abort(c(
          x = "Piecewise transform for {.val {term_name}} needs matching x/y vectors (length >= 2).",
          i = "Use list(type='pwlin', x=..., y=...)."
        ))
      }
      return(function(x) {
        stats::approx(x = xk, y = yk, xout = as.numeric(x), method = "linear", rule = 2)$y
      })
    }

    cli::cli_abort(c(
      x = "Unsupported transform type for {.val {term_name}}: {.val {tf_type}}.",
      i = "Use one of identity, functional, ispline, ispline_penalised, ispline_expit, ispline_expit_penalised, pwlin."
    ))
  }

  .sim_make_baseline_hazard <- function(spec, formula_basehaz = NULL) {
    tag_truth <- function(fun, type, parameters = list(), formula = NULL, coefficients = NULL) {
      attr(fun, "basehaz_truth") <- list(
        type = type,
        parameters = parameters,
        formula = formula,
        coefficients = coefficients
      )
      fun
    }
    if (!is.null(formula_basehaz)) {
      f_bh <- stats::as.formula(formula_basehaz)
      t_ref <- unique(c(0, time_cens / 2, time_cens))
      baseline_clock_variables <- unique(c(time_var, event_time_var)) # accepted names for the common original-time baseline clock
      baseline_reference_data <- data.frame(.joinme_reference_time = t_ref)
      for (baseline_clock_variable in baseline_clock_variables) {
        baseline_reference_data[[baseline_clock_variable]] <- t_ref
      }
      baseline_template <- .make_model_matrix_template(
        f_bh,
        baseline_reference_data
      ) # fitted formula transformation on original simulated time
      X_ref <- .mm(baseline_template, baseline_reference_data)
      coef_bh <- .sim_resolve_regression_component(
        .sim_population_component(simulation_truth_request, "baseline"),
        prior_specification$baseline,
        colnames(X_ref),
        "baseline"
      ) # formula-based baseline-hazard population coefficients

      fun <- function(t) {
        df_t <- data.frame(.joinme_reference_time = as.numeric(t))
        for (baseline_clock_variable in baseline_clock_variables) {
          df_t[[baseline_clock_variable]] <- as.numeric(t)
        }
        X_bh <- .mm(baseline_template, df_t)
        lp <- as.numeric(X_bh %*% coef_bh)
        lp <- pmin(lp, log(.Machine$double.xmax) - 2)
        as.numeric(exp(lp))
      }
      return(tag_truth(fun, "formula", formula = f_bh, coefficients = coef_bh))
    }

    if (is.null(spec)) {
      spec <- list(type = "weibull")
    }
    if (is.character(spec)) {
      spec <- list(type = spec)
    }
    if (!is.list(spec) || is.null(spec$type)) {
      cli::cli_abort(c(
        x = "{.arg truth$basehaz} must be a character mode or named list with {.arg type}.",
        i = "Supported types: constant, linear, piecewise, weibull, spline."
      ))
    }

    mode <- tolower(as.character(spec$type)[1])
    if (mode == "constant") {
      rate <- as.numeric(spec$rate %||% spec$lambda %||% stats::runif(1, min = 0.01, max = 0.5))
      if (isTRUE(spec$log)) {
        rate <- exp(rate)
      }
      fun <- function(t) rep(rate, length(t))
      return(tag_truth(
        fun,
        "constant",
        parameters = list(rate = rate),
        formula = ~1,
        coefficients = stats::setNames(log(rate), "(Intercept)")
      ))
    }

    if (mode == "linear") {
      intercept <- as.numeric(spec$intercept %||% stats::rnorm(1, mean = -2.2, sd = 0.3))
      slope <- as.numeric(spec$slope %||% stats::rnorm(1, mean = 0.25, sd = 0.1))
      fun <- function(t) {
        t <- as.numeric(t)
        # softplus(intercept + slope * pmax(t, 0))
        exp(intercept + slope * pmax(t, 0))
      }
      return(tag_truth(
        fun,
        "linear",
        parameters = list(intercept = intercept, slope = slope),
        formula = ~1 + time,
        coefficients = c("(Intercept)" = intercept, time = slope)
      ))
    }

    if (mode %in% c("piecewise", "pwlin", "piecewise_linear")) {
      breaks <- sort(as.numeric(spec$breaks %||% c(2, 5)))
      rates <- as.numeric(spec$rates %||% sort(stats::runif(length(breaks) + 1L, min = 0.03, max = 0.3)))
      if (length(rates) != length(breaks) + 1L) {
        cli::cli_abort(c(
          x = "Piecewise baseline requires length(rates) = length(breaks) + 1.",
          i = "Got {.val {length(rates)}} rates and {.val {length(breaks)}} breaks."
        ))
      }
      fun <- function(t) {
        t <- pmax(as.numeric(t), 0)
        idx <- findInterval(t, vec = breaks, rightmost.closed = TRUE) + 1L
        rates[idx]
      }
      return(tag_truth(fun, "piecewise", parameters = list(breaks = breaks, rates = rates)))
    }

    if (mode == "weibull") {
      shape <- as.numeric(spec$shape %||% stats::runif(1, min = 0.8, max = 1.8))
      scale <- as.numeric(spec$scale %||% stats::runif(1, min = 4, max = 8))
      return(tag_truth(
        weibull_h0(shape = shape, scale = scale),
        "weibull",
        parameters = list(shape = shape, scale = scale)
      ))
    }

    if (mode %in% c("spline", "bs", "ns")) {
      basis_type <- tolower(as.character(spec$basis %||% if (mode == "ns") "ns" else "bs")[1])
      knots <- as.numeric(spec$knots %||% stats::quantile(seq(0, time_cens, length.out = 100), probs = c(0.25, 0.5, 0.75)))
      degree <- as.integer(spec$degree %||% 3L)
      intercept <- as.logical(spec$intercept %||% TRUE)

      build_basis <- function(t) {
        if (basis_type == "ns") {
          splines::ns(t, knots = knots, Boundary.knots = c(0, time_cens), intercept = intercept)
        } else {
          splines::bs(t, knots = knots, Boundary.knots = c(0, time_cens), degree = degree, intercept = intercept)
        }
      }

      # Resolve spline coefficients once so baseline hazard is deterministic over
      # repeated evaluations during integration/root-finding.
      t_ref <- unique(c(0, time_cens / 2, time_cens))
      B_ref <- as.matrix(build_basis(t_ref))
      basis_coefficient_names <- paste0("basis_", seq_len(ncol(B_ref))) # stable scientific labels independent of spline library column labels
      baseline_population_request <- if (
        is.null(simulation_truth_request$baseline) &&
          is.null(simulation_truth_request$intercept) &&
          is.null(simulation_truth_request$slope)
      ) {
        prior_normal(mu = c(-2.2, rep(0, max(0L, ncol(B_ref) - 1L))), scale = 0.25)
      } else {
        .sim_population_component(simulation_truth_request, "baseline")
      } # conservative default spline log-hazard centred on a low baseline rate
      coef_vec <- .sim_resolve_regression_component(
        baseline_population_request,
        prior_specification$baseline,
        basis_coefficient_names,
        "baseline"
      ) # spline log-hazard coefficients generated from the common baseline declaration

      fun <- function(t) {
        t <- pmax(as.numeric(t), 0)
        B <- as.matrix(build_basis(t))
        lp <- as.numeric(B %*% coef_vec)
        lp <- pmin(lp, log(.Machine$double.xmax) - 2)
        as.numeric(exp(lp))
      }
      return(tag_truth(
        fun,
        basis_type,
        parameters = list(knots = knots, degree = degree, intercept = intercept),
        coefficients = coef_vec
      ))
    }

    cli::cli_abort(c(
      x = "Unsupported baseline hazard mode: {.val {mode}}.",
      i = "Use one of: constant, linear, piecewise, weibull, spline."
    ))
  }

  # ---- Step 1: Define marker/family dimensions and parse model structure
  if (.is_single_family_spec(families)) {
    n_markers_requested <- length(marker_levels %||% "m1")
    families <- rep(list(families), n_markers_requested)
  } else if (
    length(families) == 1L &&
      !is.null(marker_levels) &&
      length(marker_levels) > 1L
  ) {
    families <- rep(families, length(marker_levels))
  }
  D <- length(families)
  if (D < 1) {
    cli::cli_abort(c(
      x = "{.arg families} must include at least one marker family.",
      i = "Example: families = c('gaussian', 'student_t')."
    ))
  }
  if (is.null(marker_levels)) marker_levels <- paste0("m", seq_len(D))
  if (length(marker_levels) != D) {
    cli::cli_abort(c(
      x = "{.arg marker_levels} length must match {.arg families}.",
      i = "Expected {D}, got {length(marker_levels)}."
    ))
  }
  # ---- Step 1a: Resolve marker-family mappings and inverse-link bytecode
  family_spec <- .sim_resolve_family_specs(families, D)
  family_codes <- family_spec$family_codes
  link_names <- family_spec$link_names
  inv_link_specs <- family_spec$inv_link_specs
  use_tau_fixed <- family_spec$use_tau_fixed
  tau_fixed <- family_spec$tau_fixed
  family_names <- vapply(family_codes, .family_code_to_name, character(1))
  skew_laplace_marker <- family_codes == 9L

  # Resolve distributional formulas before generating any random quantities.
  # This permits immediate rejection of a tau regression that no simulated
  # marker can use, avoiding unnecessary longitudinal and event-time sampling.
  dist_formulas <- .normalize_formula_dist(formulaDist)
  family_names_present <- vapply(
    sort(unique(family_codes)),
    .family_code_to_name,
    character(1)
  )
  .validate_dist_formula_scopes(dist_formulas, family_names_present)
  if (
    !is.null(dist_formulas$tau) &&
      any(skew_laplace_marker) &&
      all(use_tau_fixed[skew_laplace_marker] == 1L)
  ) {
    cli::cli_abort(c(
      x = "The {.code tau ~ ...} distributional regression has no simulated skew-Laplace marker.",
      i = "Remove the tau regression or omit {.arg tau} from at least one skew-Laplace family specification."
    ))
  }

  # ---- Step 1b: Parse longitudinal random-effect structure from formulaLong
  f_exp <- reformulas::expandDoubleVerts(formulaLong)
  bars <- reformulas::findbars(f_exp)
  if (length(bars) == 0) {
    cli::cli_abort(c(
      x = "{.arg formulaLong} must include random-effects terms.",
      i = "Include at least ( ... | id ) and typically a marker block."
    ))
  }
  fixed_formula <- reformulas::nobars(f_exp)
  fixed_rhs <- stats::update(fixed_formula, . ~ .)
  fixed_rhs[[2]] <- NULL

  grp_names <- vapply(bars, function(b) .group_name_from_expr(b[[3]]), character(1))
  indep_flags <- .get_re_independence(formulaLong, marker_var = marker_var, id_var = id_var)
  id_idx <- which(grp_names == id_var)
  id_rhs_list <- if (length(id_idx) > 0) .bar_terms_to_rhs_list(bars[id_idx]) else list()

  nested_terms <- .extract_nested_marker_terms(formulaLong, marker_var = marker_var, id_var = id_var)
  mk_rhs_list <- nested_terms$mk_rhs_list
  idm_rhs_list <- nested_terms$idm_rhs_list

  # ---- Step 2: Generate exogenous event-level covariates and event design
  dataEvent <- data.frame(id = seq_len(n_id), stringsAsFactors = FALSE)
  dataEvent <- .sim_eval_covariate_formulas(n_subjects = n_id, formulas = covariate_formulas, base_df = dataEvent)
  tv_cov_specs <- attr(dataEvent, "time_varyring_specs") %||% list()
  attr(dataEvent, "time_varyring_specs") <- NULL

  surv_lhs_arity <- .sim_surv_lhs_arity(formulaEvent)
  split_required <- (isTRUE(surv_lhs_arity >= 3L) || length(tv_cov_specs) > 0L || left_truncation_max > 0)
  event_layout <- if (split_required) "split" else "subject"

  if (identical(event_layout, "split") && !isTRUE(surv_lhs_arity >= 3L)) {
    cli::cli_abort(c(
      x = "Automatic split event rows require counting-process survival syntax in {.arg formulaEvent}.",
      i = "Use {.code survival::Surv(time_start, time_stop, event) ~ ...}."
    ))
  }
  names(dataEvent)[names(dataEvent) == "id"] <- id_var

  # Covariance-regression design for the subject-specific marker-by-id covariance.
  formulaVCov <- .get_vcov_formula(
    formulaVCov = formulaVCov,
    default = ~ 1,
    context = "simulate_joinme()"
  )
  vcov_design <- .build_vcov_design(
    formulaVCov = formulaVCov,
    dataEvent = dataEvent,
    time_var = time_var,
    context = "simulate_joinme()"
  )
  K_cov_sd <- vcov_design$K_cov_sd # standard-deviation predictor count shared with fitting
  Xcov_sd <- vcov_design$Xcov_sd # standard-deviation design in fitted subject order
  K_cov_corr <- vcov_design$K_cov_corr # off-diagonal correlation predictor count shared with fitting
  Xcov_corr <- vcov_design$Xcov_corr # correlation design in fitted subject order

  # ---- Build prototype matrices for deterministic coefficient alignment
  # Longitudinal formulae are evaluated on the original observation-time grid.
  # Explicit spline knots and boundaries therefore have their declared study-
  # time meaning, while implicit spline attributes are learnt from the same
  # scheduled time support that supplies the simulated longitudinal records.
  time_scale_internal <- max(c(time_cens, times_obs), na.rm = TRUE)
  if (!is.finite(time_scale_internal) || time_scale_internal <= 0) {
    time_scale_internal <- 1.0
  }

  longitudinal_reference_times <- sort(unique(as.numeric(times_obs)))
  prototype <- dataEvent[
    rep(1L, D * length(longitudinal_reference_times)),
    ,
    drop = FALSE
  ]
  prototype[[marker_var]] <- factor(
    rep(marker_levels, each = length(longitudinal_reference_times)),
    levels = marker_levels
  )
  prototype[[time_var]] <- rep(longitudinal_reference_times, times = D)

  fixed_template <- .make_model_matrix_template(
    fixed_rhs,
    prototype
  )
  id_templates <- lapply(id_rhs_list, function(rhs) {
    .make_model_matrix_template(rhs, prototype)
  })
  mk_templates <- lapply(mk_rhs_list, function(rhs) {
    .make_model_matrix_template(rhs, prototype)
  })
  idm_templates <- lapply(idm_rhs_list, function(rhs) {
    .make_model_matrix_template(rhs, prototype)
  })

  X_proto <- .mm(fixed_template, prototype)
  Z_id_proto <- .sim_rhs_matrix(id_templates, prototype)
  Z_mk_proto <- .sim_rhs_matrix(mk_templates, prototype)
  Z_idm_proto <- .sim_rhs_matrix(idm_templates, prototype)

  beta_long <- .sim_resolve_regression_component(
    .sim_population_component(simulation_truth_request, "longitudinal"),
    prior_specification$longitudinal,
    colnames(X_proto),
    "longitudinal"
  ) # population longitudinal coefficients drawn from prior_*() declarations or fixed by numeric role declarations
  beta_long_internal <- beta_long

  assoc_from_formula <- .sim_assoc_from_formula(formulaAssoc)
  assoc_effective <- if (!is.null(assoc_from_formula)) assoc_from_formula else assoc
  assoc_effective <- unique(assoc_effective)
  if (length(assoc_effective) == 0) assoc_effective <- "cv_total"

  marker_weight_spec <- .get_marker_weight_structure(
    marker_weight_offsets = marker_weight_offsets,
    marker_levels = marker_levels,
    active_terms = .active_weighted_assoc_terms(assoc_effective),
    marker_weight_sets_shared = marker_weight_sets_shared,
    estimate_marker_weights = !isTRUE(constant_marker_weights),
    context = "simulate_joinme()"
  )
  marker_weights_offset_by_term <- marker_weight_spec$offset_by_term
  marker_mean_names <- if (isTRUE(marker_weight_spec$marker_weight_sets_shared)) "shared" else marker_weight_spec$active_terms # population location labels in compact fitted-set order
  marker_weight_means <- if (isTRUE(constant_marker_weights) || length(marker_weight_spec$active_terms) == 0L) {
    rep(0, marker_weight_spec$n_sets)
  } else {
    raw_marker_weight_request <- simulation_truth_request$marker_weights %||% list() # unnormalised marker-weight declarations
    .sim_resolve_prior_or_fixed(
      raw_marker_weight_request$intercept %||% simulation_truth_request$intercept %||% prior_specification$marker_weights$intercept,
      marker_mean_names,
      "marker_weights$intercept"
    )
  } # population marker-weight locations generated by the intercept declaration; family remains reserved for departures
  names(marker_weight_means) <- if (isTRUE(marker_weight_spec$marker_weight_sets_shared)) {
    "shared"
  } else {
    marker_weight_spec$active_terms
  }
  z_marker_weight_sets <- matrix(
    0,
    nrow = marker_weight_spec$n_sets,
    ncol = D,
    dimnames = dimnames(marker_weight_spec$offset_matrix)
  ) # realised direct departures z_sd used in effective weights
  standardised_departure_sets <- z_marker_weight_sets # family-transformed unit-scale departures retained for direct comparison with Stan's prior effect
  marker_weight_prior_draw <- list(
    departure = numeric(0),
    raw = numeric(0),
    df = numeric(0),
    df_was_fitted = logical(0),
    local_scale = numeric(0),
    global_scale = numeric(0),
    slab_multiplier = numeric(0)
  ) # realised departure hierarchy, empty when marker weights are fixed or inactive
  if (!isTRUE(constant_marker_weights) && length(marker_weight_spec$active_terms) > 0L) {
    marker_weight_prior_draw <- .sim_draw_marker_weight_prior(
      marker_weight_spec$n_sets * D,
      prior = prior_specification$marker_weights$family,
      n_sets = marker_weight_spec$n_sets
    )
    standardised_departure_sets[] <- matrix(
      marker_weight_prior_draw$departure,
      nrow = marker_weight_spec$n_sets,
      ncol = D,
      byrow = TRUE
    ) # compact set-by-marker direct departures matching contiguous Stan set segments
    z_marker_weight_sets[] <- standardised_departure_sets
  }
  marker_weights_eff <- marker_weights_offset_by_term
  marker_weights_latent_by_term <- marker_weights_offset_by_term
  marker_weight_mean_by_term <- stats::setNames(vector("list", length(marker_weights_offset_by_term)), names(marker_weights_offset_by_term)) # term-indexed view of common locations retained in truth
  for (term_key in names(marker_weights_latent_by_term)) {
    marker_weights_latent_by_term[[term_key]][] <- 0
    marker_weight_mean_by_term[[term_key]] <- 0
    set_pos <- marker_weight_spec$set_index[[term_key]]
    if (set_pos > 0L) {
      marker_weight_mean_by_term[[term_key]] <- marker_weight_means[[set_pos]]
      marker_weights_latent_by_term[[term_key]] <- as.numeric(z_marker_weight_sets[set_pos, ])
      marker_weights_eff[[term_key]] <-
        as.numeric(marker_weights_offset_by_term[[term_key]]) +
        marker_weight_means[[set_pos]] +
        marker_weights_latent_by_term[[term_key]]
    }
    names(marker_weights_latent_by_term[[term_key]]) <- marker_levels
    names(marker_weights_eff[[term_key]]) <- marker_levels
  }
  .sim_public_marker_weights <- function(by_term) {
    if (isTRUE(marker_weight_spec$marker_weight_sets_shared)) {
      term_key <- if (length(marker_weight_spec$active_terms) > 0L) {
        marker_weight_spec$active_terms[[1L]]
      } else {
        .weighted_assoc_term_keys()[[1L]]
      }
      return(by_term[[term_key]])
    }
    by_term[marker_weight_spec$active_terms]
  }
  marker_weights_offset <- .sim_public_marker_weights(marker_weights_offset_by_term)
  marker_weights_latent <- .sim_public_marker_weights(marker_weights_latent_by_term)
  marker_weights_effective <- .sim_public_marker_weights(marker_weights_eff)

  .sim_time_varyring_value <- function(spec_i, t_val) {
    breaks <- as.numeric(spec_i$breaks)
    values <- as.numeric(spec_i$values)
    t_use <- min(max(as.numeric(t_val), breaks[1]), breaks[length(breaks)])
    seg <- findInterval(t_use, breaks, rightmost.closed = TRUE, all.inside = TRUE)
    seg <- min(max(seg, 1L), length(values))
    values[seg]
  }

  .sim_event_row_at <- function(i, t_val) {
    row <- dataEvent[i, , drop = FALSE]
    if (length(tv_cov_specs) > 0L) {
      for (nm in names(tv_cov_specs)) {
        row[[nm]] <- .sim_time_varyring_value(tv_cov_specs[[nm]][[i]], t_val)
      }
    }
    row
  }

  event_reference_times <- sort(unique(c(0, as.numeric(times_obs), as.numeric(time_cens)))) # original-time support defining spline knots and boundaries
  event_clock_variables <- unique(c(time_var, event_time_var)) # formula variables representing the common original-time event clock
  event_formula_training_data <- do.call(rbind, lapply(seq_len(n_id), function(i) {
    do.call(rbind, lapply(event_reference_times, function(reference_time) {
      row <- .sim_event_row_at(i, reference_time)
      for (event_clock_variable in event_clock_variables) {
        row[[event_clock_variable]] <- reference_time
      }
      row # spline attributes are learnt once from the complete original-time support
    }))
  }))
  event_design_template <- .make_event_model_matrix_template(
    formulaEvent,
    event_formula_training_data
  )
  event_template <- do.call(rbind, lapply(seq_len(n_id), function(i) {
    row <- .sim_event_row_at(i, 0)
    for (event_clock_variable in event_clock_variables) {
      row[[event_clock_variable]] <- 0
    }
    row
  }))
  W_event <- .mm_event(event_design_template, event_template)

  raw_survival_request <- .sim_population_component(simulation_truth_request, "survival", "slope") # population event-regression declaration with global fallback
  beta_event <- .sim_resolve_prior_or_fixed(
    raw_survival_request$slope %||% prior_specification$survival$slope,
    colnames(W_event),
    "survival$slope"
  ) # event-regression population coefficients in formulaEvent model-matrix order

  # ---- Step 3: Draw hierarchical random effects and covariance-regression factors
  K_id <- ncol(Z_id_proto)
  K_mk <- ncol(Z_mk_proto)
  K_idm <- ncol(Z_idm_proto)

  .validate_assoc_cov_structure(
    assoc = assoc_effective,
    q_idm = K_idm,
    diagonal_only = indep_flags$indep_idmarker_cov,
    context = "simulate_joinme()"
  )

  re_id_effective <- .sim_resolve_re_block_cfg(
    re_params[["id", exact = TRUE]] %||% list(),
    K_id,
    "id"
  )
  re_marker_effective <- .sim_resolve_re_block_cfg(
    re_params[["marker", exact = TRUE]] %||% list(),
    K_mk,
    "marker"
  )
  if (K_id > 1L && as.integer(indep_flags$indep_id_re %||% 0L) == 1L) {
    id_corr_offdiag <- re_id_effective$corr
    diag(id_corr_offdiag) <- 0
    if (any(abs(id_corr_offdiag) > 1e-8)) {
      cli::cli_warn(c(
        x = "Top-level {.code || id} requests independent subject-level random effects in {.fn simulate_joinme}.",
        i = "Ignoring nonzero off-diagonal entries in {.arg re_params$id$corr}."
      ))
    }
    re_id_effective$corr <- diag(K_id)
    re_id_effective$cov <- diag(as.numeric(re_id_effective$sd), K_id, K_id)^2
    re_id_effective$Lcorr <- diag(K_id)
  }
  re_id_effective_internal <- re_id_effective
  re_marker_effective_internal <- re_marker_effective
  # ---- Step 3a: prepare the optional latent-progress distribution
  # The mixture builder is invoked only after the random-effect formulae have
  # determined every available dimension. This is the earliest point at which
  # `class_dimensions` can be checked faithfully. The builder itself reuses
  # the fitting coordinate-layout machinery, so simulation and estimation
  # cannot silently disagree about which standardised coordinate an index
  # denotes.
  mixture_configuration <- .sim_prepare_mixture(
    specification = .mixture_specification,
    data_event = dataEvent,
    marker_prototype = prototype,
    dimensions = c(
      subject = K_id,
      marker = K_mk,
      covariance_basis = K_idm
    ),
    labels = list(
      subject = colnames(Z_id_proto) %||% character(0),
      marker = colnames(Z_mk_proto) %||% character(0),
      covariance = colnames(Z_idm_proto) %||% character(0)
    ),
    covariance_is_diagonal =
      as.integer(indep_flags$indep_idmarker_cov %||% 0L) == 1L,
    shrinkage = shrinkage,
    id_variable = id_var,
    marker_variable = marker_var
  ) # checked component layout, probabilities, allocations, locations and scales

  # Draw the standardised individual effects first, replace only the selected
  # coordinates conditional on the subject's common class, and finally apply
  # the ordinary random-effect covariance factor. This is precisely the
  # `u_i = L_u z_{u,i}` ordering used in Stan.
  z_id <- if (K_id > 0L) {
    matrix(stats::rnorm(n_id * K_id), nrow = n_id, ncol = K_id)
  } else {
    matrix(0, nrow = n_id, ncol = 0L)
  } # standardised individual random effects before covariance scaling
  z_id <- .sim_apply_mixture_to_latent(
    latent_matrix = z_id,
    level = "subject",
    allocation = mixture_configuration$allocation$subject %||% integer(0),
    mixture = mixture_configuration,
    shrinkage = shrinkage
  ) # selected component-conditional individual coordinates
  re_id <- if (K_id > 0L) {
    z_id %*% chol(re_id_effective_internal$cov)
  } else {
    matrix(0, nrow = n_id, ncol = 0L)
  } # realised subject effects on the original study-time basis

  # Repeat the same construction for marker-level standardised effects. Marker
  # weights remain population association quantities and do not introduce a
  # second allocation label.
  z_marker <- if (K_mk > 0L) {
    matrix(stats::rnorm(D * K_mk), nrow = D, ncol = K_mk)
  } else {
    matrix(0, nrow = D, ncol = 0L)
  } # standardised marker random effects before covariance scaling
  z_marker <- .sim_apply_mixture_to_latent(
    latent_matrix = z_marker,
    level = "marker",
    allocation = mixture_configuration$allocation$marker %||% integer(0),
    mixture = mixture_configuration,
    shrinkage = shrinkage
  ) # selected component-conditional marker coordinates
  re_marker <- if (K_mk > 0L) {
    z_marker %*% chol(re_marker_effective_internal$cov)
  } else {
    matrix(0, nrow = D, ncol = 0L)
  } # realised marker effects on the original study-time basis

  raw_vcov_request <- simulation_truth_request$vcov %||% list() # fixed values or generating distributions for covariance-regression parameters
  if (inherits(raw_vcov_request, "joinme_prior_spec") || is.numeric(raw_vcov_request) ||
      any(names(raw_vcov_request) %||% character(0) %in% c("intercept", "slope", "latent"))) {
    raw_vcov_request <- list(sd = raw_vcov_request, corr = raw_vcov_request)
  } # a shared declaration applies to both marginal-scale and correlation regressions
  raw_vcov_component <- function(part_name) {
    component <- raw_vcov_request[[part_name]] %||% list() # declarations specific to the selected covariance part
    if (inherits(component, "joinme_prior_spec") || is.numeric(component)) {
      component <- list(intercept = component, slope = component, latent = component)
    }
    if (!is.list(component)) component <- list()
    list(
      intercept = component$intercept %||% simulation_truth_request$intercept %||% prior_specification$vcov[[part_name]]$intercept,
      slope = component$slope %||% simulation_truth_request$slope %||% prior_specification$vcov[[part_name]]$slope,
      latent = component$latent %||% simulation_truth_request$slope %||% prior_specification$vcov[[part_name]]$latent
    ) # complete population declarations with global and checked defaults
  }
  re_idm_effective <- list(
    mode = "iid_standard_normal",
    sd = if (K_idm > 0) rep(1, K_idm) else numeric(0),
    corr = if (K_idm > 0) diag(K_idm) else matrix(0.0, 0, 0)
  )
  re_idm_flat <- if (K_idm > 0) {
    .sim_draw_mvn(n_id * D, diag(K_idm))
  } else {
    matrix(0.0, nrow = n_id * D, ncol = 0)
  }

  re_idm <- if (K_idm > 0) {
    aperm(array(re_idm_flat, dim = c(D, n_id, K_idm)), c(2, 1, 3))
  } else {
    array(0.0, dim = c(n_id, D, 0))
  }

  M_cov <- if (K_idm > 0) {
    if (as.integer(indep_flags$indep_idmarker_cov %||% 0L) == 1L) K_idm else K_idm * (K_idm + 1L) / 2L
  } else {
    0L
  }

  idx_row_cov <- integer(M_cov)
  idx_col_cov <- integer(M_cov)
  if (M_cov > 0) {
    if (as.integer(indep_flags$indep_idmarker_cov %||% 0L) == 1L) {
      idx_row_cov <- seq_len(M_cov)
      idx_col_cov <- seq_len(M_cov)
    } else {
      pos <- 1L
      for (r in seq_len(K_idm)) {
        for (c in seq_len(r)) {
          idx_row_cov[pos] <- r
          idx_col_cov[pos] <- c
          pos <- pos + 1L
        }
      }
    }
  }

  .sim_cov_lp_to_matrix <- function(lp_vec, q_idm, idx_row, idx_col, diag_link) {
    .cov_lp_to_chol(
      lp_vec = lp_vec,
      q_idm = q_idm,
      idx_row = idx_row,
      idx_col = idx_col,
      diag_link = diag_link
    )
  }

  .sim_canonicalize_cov_latent <- function(lambda, z) {
    lambda <- as.numeric(lambda)
    if (length(lambda) == 0L) {
      return(list(
        lambda = numeric(0),
        z = matrix(0.0, nrow = nrow(z), ncol = 0L),
        sign = numeric(0)
      ))
    }

    sign_vec <- ifelse(lambda < 0, -1, 1)
    lambda_out <- abs(lambda)
    z_out <- z[, seq_along(lambda), drop = FALSE]
    z_out <- sweep(z_out, 2L, sign_vec, `*`)

    list(lambda = lambda_out, z = z_out, sign = sign_vec)
  }

  diag_link_cov <- tolower(as.character(vcov_diag_link)[1])
  if (!diag_link_cov %in% c("softplus", "exp")) {
    cli::cli_abort(c(
      x = "{.arg vcov_diag_link} must be 'softplus' or 'exp'.",
      i = "Use 'softplus' (default) or 'exp'."
    ))
  }

  # STEP 3b-i: separate the scientific SD and correlation parameter blocks.
  #
  # A full lower triangle is still retained internally because class selection
  # and association code use that stable coordinate order. The public
  # generative syntax, however, mirrors formulaVCov directly: `sd` has Q rows
  # and `corr` has Q(Q-1)/2 rows. Superseded packed inputs are rejected because
  # they cannot express distinct standard-deviation and correlation formulae.
  M_corr_cov <- if (K_idm >= 2L && as.integer(indep_flags$indep_idmarker_cov %||% 0L) == 0L) {
    as.integer(K_idm * (K_idm - 1L) / 2L)
  } else {
    0L
  } # number of off-diagonal partial-correlation coordinates actually fitted
  diagonal_positions <- which(idx_row_cov == idx_col_cov) # packed locations governed by formulaVCov$sd
  correlation_positions <- which(idx_row_cov != idx_col_cov) # packed locations governed by formulaVCov$corr
  re_cov_sd_cfg <- raw_vcov_component("sd") # marginal-scale generating declarations
  re_cov_corr_cfg <- raw_vcov_component("corr") # off-diagonal correlation generating declarations
  sd_coordinate_names <- if (K_idm > 0L) {
    paste0("sd[", seq_len(K_idm), "]")
  } else {
    character(0)
  } # names of the fitted marginal-scale coordinates; explicitly empty when no marker-by-subject effect is present
  correlation_coordinate_names <- if (M_corr_cov > 0L) {
    paste0("corr[", seq_len(M_corr_cov), "]")
  } else {
    character(0)
  } # names of the fitted off-diagonal coordinates; avoids paste0() turning a zero-length index into the spurious label "corr[]"

  alpha_cov_sd <- .sim_resolve_prior_or_fixed(
    re_cov_sd_cfg$intercept, sd_coordinate_names, "vcov$sd$intercept"
  ) # SD-regression intercept for every marker-by-subject basis coordinate
  alpha_cov_corr <- .sim_resolve_prior_or_fixed(
    re_cov_corr_cfg$intercept, correlation_coordinate_names, "vcov$corr$intercept"
  ) # correlation-regression intercept for every row-major off-diagonal coordinate
  lambda_cov_sd <- .sim_resolve_prior_or_fixed(
    re_cov_sd_cfg$latent, sd_coordinate_names, "vcov$sd$latent"
  ) # unexplained subject heterogeneity loading on each unlinked SD predictor
  lambda_cov_corr <- .sim_resolve_prior_or_fixed(
    re_cov_corr_cfg$latent, correlation_coordinate_names, "vcov$corr$latent"
  ) # unexplained subject heterogeneity loading on each Fisher-like tanh predictor

  for (value_name in c("alpha_cov_sd", "alpha_cov_corr", "lambda_cov_sd", "lambda_cov_corr")) {
    if (any(!is.finite(get(value_name)))) {
      cli::cli_abort("Every covariance-regression intercept and loading must be finite.")
    }
  }
  beta_cov_sd <- matrix(.sim_resolve_prior_or_fixed(
    re_cov_sd_cfg$slope,
    as.vector(outer(sd_coordinate_names, colnames(Xcov_sd), paste, sep = ":")),
    "vcov$sd$slope"
  ), nrow = K_idm, ncol = K_cov_sd) # Q by K_cov_sd systematic SD effects
  beta_cov_corr <- matrix(.sim_resolve_prior_or_fixed(
    re_cov_corr_cfg$slope,
    as.vector(outer(correlation_coordinate_names, colnames(Xcov_corr), paste, sep = ":")),
    "vcov$corr$slope"
  ), nrow = M_corr_cov, ncol = K_cov_corr) # M_corr by K_cov_corr systematic correlation effects

  alpha_cov <- numeric(M_cov) # compatibility/reporting vector in packed lower-triangular order
  lambda_cov <- numeric(M_cov) # compatibility/reporting loading vector in the same packed order
  alpha_cov[diagonal_positions] <- alpha_cov_sd
  alpha_cov[correlation_positions] <- alpha_cov_corr
  lambda_cov[diagonal_positions] <- lambda_cov_sd
  lambda_cov[correlation_positions] <- lambda_cov_corr
  L_i <- array(0.0, dim = c(n_id, K_idm, K_idm))
  z_cov <- matrix(0.0, nrow = n_id, ncol = max(1L, M_cov))
  lambda_cov_sign <- rep(1, M_cov)
  if (M_cov > 0 && K_idm > 0) {
    z_cov_raw <- matrix(stats::rnorm(n_id * M_cov), nrow = n_id, ncol = M_cov)
    cov_latent <- .sim_canonicalize_cov_latent(lambda_cov, z_cov_raw)
    lambda_cov <- cov_latent$lambda
    z_cov <- cov_latent$z
    lambda_cov_sign <- cov_latent$sign
    lambda_cov_sd <- lambda_cov[diagonal_positions]
    lambda_cov_corr <- lambda_cov[correlation_positions]

    # The covariance-regression mixture acts on Stan's canonical `z_L`
    # coordinates, after any sign absorbed from a supplied negative loading.
    # Individual and covariance-regression classes share the same subject
    # allocation held in `mixture_configuration$allocation$subject`.
    covariance_class_type <- intersect(
      c("corr", "vcov"),
      mixture_configuration$class_type %||% character(0)
    ) # public covariance representation selected for the shared latent block
    z_cov <- .sim_apply_mixture_to_latent(
      latent_matrix = z_cov,
      level = if (length(covariance_class_type) == 1L) {
        covariance_class_type[[1L]]
      } else {
        ""
      },
      allocation =
        mixture_configuration$allocation$subject %||% integer(0),
      mixture = mixture_configuration,
      shrinkage = shrinkage
    ) # selected component-conditional covariance-regression latents

    # STEP 3b-ii: evaluate the two independent observed-covariate regressions
    # and place their predictors into the unchanged packed latent order.
    lp_cov <- matrix(alpha_cov, nrow = n_id, ncol = M_cov, byrow = TRUE)
    if (K_cov_sd > 0L) {
      lp_cov[, diagonal_positions] <- lp_cov[, diagonal_positions, drop = FALSE] +
        Xcov_sd %*% t(beta_cov_sd)
    }
    if (K_cov_corr > 0L && M_corr_cov > 0L) {
      lp_cov[, correlation_positions] <- lp_cov[, correlation_positions, drop = FALSE] +
        Xcov_corr %*% t(beta_cov_corr)
    }
    lp_cov <- lp_cov + sweep(z_cov, 2L, lambda_cov, `*`)

    for (i in seq_len(n_id)) {
      L_i[i, , ] <- .sim_cov_lp_to_matrix(
        lp_vec = lp_cov[i, ],
        q_idm = K_idm,
        idx_row = idx_row_cov,
        idx_col = idx_col_cov,
        diag_link = diag_link_cov
      )
    }
  }

  L_i_eff <- L_i

  re_idm_scaled <- array(0.0, dim = c(n_id, D, K_idm))
  if (K_idm > 0) {
    for (i in seq_len(n_id)) {
      Li <- matrix(L_i_eff[i, , ], K_idm, K_idm)
      re_idm_scaled[i, , ] <- re_idm[i, , ] %*% t(Li)
    }
  }

  label_re_matrix <- function(mat, row_labels, col_labels) {
    if (!is.matrix(mat)) {
      return(mat)
    }
    row_labels <- as.character(row_labels %||% character(0))
    col_labels <- as.character(col_labels %||% character(0))
    if (nrow(mat) == length(row_labels)) {
      rownames(mat) <- row_labels
    }
    if (ncol(mat) == length(col_labels)) {
      colnames(mat) <- col_labels
    }
    mat
  }

  label_re_array <- function(arr, id_labels, marker_labels, term_labels) {
    if (is.null(arr) || length(dim(arr)) != 3L) {
      return(arr)
    }

    dim_names <- dimnames(arr)
    if (is.null(dim_names)) {
      dim_names <- vector("list", 3L)
    }

    id_labels <- as.character(id_labels %||% character(0))
    marker_labels <- as.character(marker_labels %||% character(0))
    term_labels <- as.character(term_labels %||% character(0))

    if (dim(arr)[1] == length(id_labels)) {
      dim_names[[1]] <- id_labels
    }
    if (dim(arr)[2] == length(marker_labels)) {
      dim_names[[2]] <- marker_labels
    }
    if (dim(arr)[3] == length(term_labels)) {
      dim_names[[3]] <- term_labels
    }

    dimnames(arr) <- dim_names
    arr
  }

  re_id_public <- re_id
  re_marker_public <- re_marker
  re_idm_public <- re_idm_scaled

  id_labels_public <- as.character(dataEvent[[id_var]] %||% seq_len(n_id))
  id_term_labels <- colnames(Z_id_proto) %||% paste0("id_re_", seq_len(ncol(re_id_public)))
  marker_term_labels <- colnames(Z_mk_proto) %||% paste0("marker_re_", seq_len(ncol(re_marker_public)))
  idm_term_labels <- colnames(Z_idm_proto) %||% paste0("id_marker_re_", seq_len(dim(re_idm_public)[3]))

  re_id <- label_re_matrix(re_id, id_labels_public, id_term_labels)
  re_id_public <- label_re_matrix(re_id_public, id_labels_public, id_term_labels)
  re_marker <- label_re_matrix(re_marker, marker_levels, marker_term_labels)
  re_marker_public <- label_re_matrix(re_marker_public, marker_levels, marker_term_labels)
  re_idm <- label_re_array(re_idm, id_labels_public, marker_levels, idm_term_labels)
  re_idm_scaled <- label_re_array(re_idm_scaled, id_labels_public, marker_levels, idm_term_labels)
  re_idm_public <- label_re_array(re_idm_public, id_labels_public, marker_levels, idm_term_labels)

  # ---- Step 3a: Build association channels entering the survival linear predictor
  if (!is.null(assoc_from_formula)) assoc <- assoc_from_formula
  assoc <- unique(assoc)
  if (length(assoc) == 0) assoc <- "cv_total"
  assoc <- .validate_assoc_channels(assoc, context = "simulate_joinme()")

  M_corr <- .assoc_cov_feature_count(K_idm, include_diag = FALSE)
  M_vcov <- if (as.integer(indep_flags$indep_idmarker_cov %||% 0L) == 1L) {
    as.integer(K_idm)
  } else {
    .assoc_cov_feature_count(K_idm, include_diag = TRUE)
  }

  assoc_coef_scalar <- c(
    cv_total = 0,
    cv_mean = 0,
    cv_marker = 0,
    cs_total = 0,
    cs_mean = 0,
    cs_marker = 0
  )
  assoc_coef_corr <- rep(0.0, M_corr)
  assoc_coef_vcov <- rep(0.0, M_vcov)

  scalar_assoc_names <- intersect(names(assoc_coef_scalar), assoc) # active scalar association coefficients in documented formula order
  complete_assoc_names <- c(
    scalar_assoc_names,
    if ("corr" %in% assoc) paste0("corr[", seq_len(M_corr), "]") else character(0),
    if ("vcov" %in% assoc) paste0("vcov[", seq_len(M_vcov), "]") else character(0)
  ) # complete population association vector including covariance components
  raw_assoc_request <- .sim_population_component(simulation_truth_request, "assoc_coef", "slope") # association population declaration with global fallback
  resolved_assoc <- .sim_resolve_prior_or_fixed(
    raw_assoc_request$slope %||% prior_specification$assoc$slope,
    complete_assoc_names,
    "assoc$slope"
  )
  if (length(scalar_assoc_names) > 0L) assoc_coef_scalar[scalar_assoc_names] <- resolved_assoc[scalar_assoc_names]
  if (M_corr > 0L && "corr" %in% assoc) assoc_coef_corr <- unname(resolved_assoc[paste0("corr[", seq_len(M_corr), "]")])
  if (M_vcov > 0L && "vcov" %in% assoc) assoc_coef_vcov <- unname(resolved_assoc[paste0("vcov[", seq_len(M_vcov), "]")])

  weighted_terms <- intersect(c("cv_total", "cs_total", "cv_marker", "cs_marker"), assoc)
  if (length(weighted_terms) > 0) {
    assoc_coef_scalar[weighted_terms] <- abs(assoc_coef_scalar[weighted_terms])
  }

  # Use signed marker weights directly (no logistic bounding), aligned with Stan.

  assoc_coef_vec <- assoc_coef_scalar[intersect(names(assoc_coef_scalar), assoc)]
  if ("corr" %in% assoc && M_corr > 0) {
    vc_names <- paste0("corr[", seq_len(M_corr), "]")
    assoc_coef_vec <- c(assoc_coef_vec, stats::setNames(assoc_coef_corr, vc_names))
  }
  if ("vcov" %in% assoc && M_vcov > 0) {
    vc_names <- paste0("vcov[", seq_len(M_vcov), "]")
    assoc_coef_vec <- c(assoc_coef_vec, stats::setNames(assoc_coef_vcov, vc_names))
  }

  transforms <- .sim_validate_fixed_transform_shifts(.sim_normalize_transform_list(transforms))
  .sim_make_component_transform_set <- function(spec, term_name, n_components) {
    # Expand a single user-facing covariance-style transform spec into a list of
    # per-component evaluators. Each component currently shares the same runtime
    # transform definition in simulation, which mirrors the common-input user
    # layout used when building standata for Stan.
    n_components <- as.integer(n_components %||% 0L)
    if (n_components <= 0L) {
      return(list())
    }
    rep(list(.sim_make_assoc_transform(spec, term_name)), n_components)
  }
  tf_funs <- list(
    cv_total = .sim_make_assoc_transform(transforms$cv_total, "cv_total"),
    cv_mean = .sim_make_assoc_transform(transforms$cv_mean, "cv_mean"),
    cv_marker = .sim_make_assoc_transform(transforms$cv_marker, "cv_marker"),
    cs_total = .sim_make_assoc_transform(transforms$cs_total, "cs_total"),
    cs_mean = .sim_make_assoc_transform(transforms$cs_mean, "cs_mean"),
    cs_marker = .sim_make_assoc_transform(transforms$cs_marker, "cs_marker"),
    corr = .sim_make_component_transform_set(transforms$corr, "corr", M_corr),
    vcov = .sim_make_component_transform_set(transforms$vcov, "vcov", M_vcov)
  )

  has_tf_cv_marker <- !is.null(transforms$cv_marker)
  has_tf_cv_total <- !is.null(transforms$cv_total)

  # ---- Step 3b: Define latent longitudinal trajectory evaluators over continuous time
 
  assoc_row_template <- vector("list", n_id)
  marker_factor <- factor(marker_levels, levels = marker_levels)
  for (i in seq_len(n_id)) {
    row_df_i <- dataEvent[rep(i, D), , drop = FALSE]
    row_df_i[[marker_var]] <- marker_factor
    assoc_row_template[[i]] <- row_df_i
  }

  eta_components_all_markers <- function(i, t) {
    row_df <- assoc_row_template[[i]]
    row_df[[time_var]] <- t

    x_fix <- .mm(fixed_template, row_df)
    z_id <- .sim_rhs_matrix(id_templates, row_df)
    z_mk <- .sim_rhs_matrix(mk_templates, row_df)
    z_idm <- .sim_rhs_matrix(idm_templates, row_df)

    fixed_part <- if (ncol(x_fix) > 0) as.numeric(x_fix %*% beta_long_internal) else rep(0, D)
    id_part <- if (ncol(z_id) > 0) as.numeric(z_id %*% re_id[i, ]) else rep(0, D)
    mk_part <- if (ncol(z_mk) > 0) rowSums(z_mk * re_marker) else rep(0, D)
    idm_part <- if (ncol(z_idm) > 0) rowSums(z_idm * re_idm_scaled[i, , ]) else rep(0, D)

    mu_mean <- fixed_part + id_part
    mu_marker <- mk_part + idm_part
    list(
      mu_mean = mu_mean,
      mu_marker = mu_marker,
      mu_total = mu_mean + mu_marker
    )
  }

  .sim_corr_features <- function(i) {
    # Build subject-level raw correlation-factor features in the same lower-
    # triangular ordering used by Stan: (2,1), (3,1), (3,2), ... .
    #
    # Important: this returns the off-diagonal entries of the subject-specific
    # Cholesky-correlation factor K_i, not pairwise correlations from Sigma_i.
    if (K_idm < 2) return(numeric(0))
    Li <- matrix(L_i[i, , ], nrow = K_idm, ncol = K_idm)
    if (!all(is.finite(Li))) return(rep(0.0, M_corr))
    .assoc_corr_features_from_chol(Li)
  }

  .sim_vcov_features <- function(i) {
    if (K_idm < 1) return(numeric(0))
    Li <- matrix(L_i_eff[i, , ], nrow = K_idm, ncol = K_idm)
    if (!all(is.finite(Li))) {
      return(rep(0.0, M_vcov))
    }
    .assoc_vcov_features_from_chol(
      Li,
      diagonal_only = as.integer(indep_flags$indep_idmarker_cov %||% 0L) == 1L
    )
  }

  corr_by_id <- lapply(seq_len(n_id), .sim_corr_features)
  vcov_by_id <- lapply(seq_len(n_id), .sim_vcov_features)

  assoc_components <- function(i, t) {
    # Step A: build marker-resolved CV components at t and t + eps_cs.
    now <- eta_components_all_markers(i, t)
    eps <- eta_components_all_markers(i, t + eps_cs)

    w_mean <- function(term_key, x) {
      weights_term <- as.numeric(marker_weights_eff[[term_key]] %||% rep(1, D))
      sum(weights_term * x) / D
    }

    # Step B: aggregate raw CV summaries (mean, marker, total).
    cv_mean_raw <- mean(now$mu_mean)
    cv_marker_raw <- w_mean("cv_marker", now$mu_marker)
    cv_total_raw <- w_mean("cv_total", now$mu_total)

    cv_mean <- tf_funs$cv_mean(cv_mean_raw)
    cv_marker <- if (has_tf_cv_marker) w_mean("cv_marker", tf_funs$cv_marker(now$mu_marker)) else cv_marker_raw
    cv_total <- if (has_tf_cv_total) w_mean("cv_total", tf_funs$cv_total(now$mu_total)) else cv_total_raw

    # Step C: compute mean slope via finite difference.
    cs_mean_raw <- (mean(eps$mu_mean) - cv_mean_raw) / eps_cs

    # Step D: marker-level slopes for marker and total components.
    # We keep these at marker resolution so CS transforms are applied per marker
    # before weighted averaging, mirroring Stan semantics.
    cs_marker_raw_by_marker <- (eps$mu_marker - now$mu_marker) / eps_cs
    cs_total_raw_by_marker <- (eps$mu_total - now$mu_total) / eps_cs

    # Step E: CS aggregation semantics (aligned with CV aggregation semantics):
    # transform first at marker level, then weighted-average across markers.
    cs_marker <- w_mean("cs_marker", tf_funs$cs_marker(cs_marker_raw_by_marker))
    cs_total <- w_mean("cs_total", tf_funs$cs_total(cs_total_raw_by_marker))

    # Step F: keep weighted raw summaries for debugging/inspection helpers.
    cs_marker_raw <- w_mean("cs_marker", cs_marker_raw_by_marker)
    cs_total_raw <- w_mean("cs_total", cs_total_raw_by_marker)

    corr_vals <- corr_by_id[[i]]
    # CORR semantics are transform-first, then weight:
    #   corr_assoc = sum_j assoc_coef_corr[j] * tf_corr(corr_vals[j])
    # This mirrors Stan and intentionally avoids tf_corr(sum_j a_j * corr_j).
    corr_vals_tf <- if (length(corr_vals) > 0) {
      vapply(seq_along(corr_vals), function(m) {
        as.numeric(tf_funs$corr[[m]](corr_vals[m]))[1] - as.numeric(tf_funs$corr[[m]](0))[1]
      }, numeric(1))
    } else {
      numeric(0)
    }
    corr_assoc <- if (length(corr_vals_tf) > 0) sum(assoc_coef_corr * corr_vals_tf) else 0
    vcov_vals <- vcov_by_id[[i]]
    # VCOV follows the same calculation layout as CORR:
    #   Step 1: extract raw time-constant Cholesky features for subject i.
    #   Step 2: transform each component separately.
    #   Step 3: apply the corresponding association coefficient after the
    #           transform, never before it.
    vcov_vals_tf <- if (length(vcov_vals) > 0) {
      vapply(seq_along(vcov_vals), function(m) {
        as.numeric(tf_funs$vcov[[m]](vcov_vals[m]))[1] - as.numeric(tf_funs$vcov[[m]](0))[1]
      }, numeric(1))
    } else {
      numeric(0)
    }
    vcov_assoc <- if (length(vcov_vals_tf) > 0) sum(assoc_coef_vcov * vcov_vals_tf) else 0

    list(
      cv_total = cv_total,
      cv_mean = cv_mean,
      cv_marker = cv_marker,
      cs_total = cs_total,
      cs_mean = tf_funs$cs_mean(cs_mean_raw),
      cs_marker = cs_marker,
      corr = corr_assoc,
      vcov = vcov_assoc,
      raw = list(
        cv_mean = cv_mean_raw,
        cs_total = cs_total_raw,
        cs_marker = cs_marker_raw,
        cs_mean = cs_mean_raw,
        corr = corr_assoc,
        corr_vals = corr_vals,
        corr_vals_tf = corr_vals_tf,
        vcov = vcov_assoc,
        vcov_vals = vcov_vals,
        vcov_vals_tf = vcov_vals_tf,
        marker_values = now
      )
    )
  }

  h0_fn <- if (is.function(h0)) {
    h0
  } else {
    .sim_make_baseline_hazard(
      spec = baseline_hazard,
      formula_basehaz = formulaBasehaz
    )
  }
  baseline_hazard_truth <- attr(h0_fn, "basehaz_truth", exact = TRUE)
  if (is.null(baseline_hazard_truth)) {
    baseline_hazard_truth <- list(
      type = "custom_function",
      parameters = list(),
      formula = NULL,
      coefficients = NULL
    )
  }

  .eta_event_it <- function(i, t) {
    row <- .sim_event_row_at(i, t)
    for (event_clock_variable in event_clock_variables) {
      row[[event_clock_variable]] <- as.numeric(t)
    }
    W_row <- .mm_event(event_design_template, row)
    w <- rep(0, length(beta_event))
    names(w) <- names(beta_event)
    common <- intersect(colnames(W_row), names(beta_event))
    if (length(common) > 0L) {
      w[common] <- as.numeric(W_row[1, common, drop = TRUE])
    }
    sum(w * beta_event)
  }

  .sim_assoc_lp_from_components <- function(comp) {
    assoc_coef_scalar[["cv_total"]] * comp$cv_total +
      assoc_coef_scalar[["cv_mean"]] * comp$cv_mean +
      assoc_coef_scalar[["cv_marker"]] * comp$cv_marker +
      assoc_coef_scalar[["cs_total"]] * comp$cs_total +
      assoc_coef_scalar[["cs_mean"]] * comp$cs_mean +
      assoc_coef_scalar[["cs_marker"]] * comp$cs_marker +
      comp$corr +
      comp$vcov
  }

  .sim_time_key <- function(t) sprintf("%.12f", as.numeric(t))
  assoc_lp_cache <- replicate(n_id, new.env(parent = emptyenv(), hash = TRUE), simplify = FALSE)
  cumhaz_cache <- replicate(n_id, new.env(parent = emptyenv(), hash = TRUE), simplify = FALSE)

  .assoc_lp_cached <- function(i, t) {
    key <- .sim_time_key(t)
    cache_i <- assoc_lp_cache[[i]]
    if (exists(key, envir = cache_i, inherits = FALSE)) {
      return(get(key, envir = cache_i, inherits = FALSE))
    }
    lp_val <- .sim_assoc_lp_from_components(assoc_components(i, t))
    assign(key, lp_val, envir = cache_i)
    lp_val
  }

  # ---- Resolve Gauss-Kronrod nodes/weights for survival integration
  gk_spec <- gk_quadrature(nodes = quadrature_nodes)

  .sim_composite_gk <- function(f, upper, panels) {
    # Composite Gauss-Kronrod with nodes defined on [0,1] per panel.
    if (upper <= 0) return(0)
    edges <- seq(0, upper, length.out = panels + 1L)
    total <- 0
    for (p in seq_len(panels)) {
      a <- edges[p]
      b <- edges[p + 1L]
      nodes <- a + (b - a) * gk_spec$nodes
      vals <- f(nodes)
      total <- total + (b - a) * sum(gk_spec$weights * vals)
    }
    total
  }

  hazard_i <- function(i, t) {
    t <- as.numeric(t)
    if (length(t) > 1L) {
      return(vapply(t, function(tt) h0_fn(tt) * exp(.eta_event_it(i, tt) + .assoc_lp_cached(i, tt)), numeric(1)))
    }
    h0_fn(t) * exp(.eta_event_it(i, t) + .assoc_lp_cached(i, t))
  }

  cumhaz_i <- function(i, t) {
    t <- as.numeric(t)
    if (length(t) > 1L) return(vapply(t, function(tt) cumhaz_i(i, tt), numeric(1)))
    if (t <= 0) return(0)

    key <- .sim_time_key(t)
    cache_i <- cumhaz_cache[[i]]
    if (exists(key, envir = cache_i, inherits = FALSE)) {
      return(get(key, envir = cache_i, inherits = FALSE))
    }

    method <- tolower(as.character(integration_control$method %||% "integrate")[1])
    rel_tol <- as.numeric(integration_control$rel.tol %||% 1e-6)
    max_panels <- max(1L, as.integer(integration_control$subdivisions %||% 512L))
    max_refine <- max(1L, as.integer(integration_control$max_refine %||% 8L))

    gk_integral <- function() {
      panels <- 1L
      prev <- NA_real_
      curr <- NA_real_
      for (iter in seq_len(max_refine)) {
        curr <- .sim_composite_gk(function(u) hazard_i(i, u), upper = t, panels = panels)
        if (is.finite(prev) && abs(curr - prev) <= rel_tol * max(1, abs(prev))) break
        prev <- curr
        if (panels >= max_panels) break
        panels <- min(max_panels, panels * 2L)
      }
      as.numeric(curr)
    }

    value <- if (identical(method, "integrate")) {
      integrate_ctrl <- integration_control
      integrate_ctrl$method <- NULL
      integrate_ctrl$max_refine <- NULL
      out <- tryCatch(
        do.call(integrate, c(list(f = function(u) hazard_i(i, u), lower = 0, upper = t), integrate_ctrl)),
        error = function(e) NULL
      )
      if (is.null(out)) gk_integral() else as.numeric(out$value)
    } else {
      gk_integral()
    }

    assign(key, value, envir = cache_i)
    value
  }

  draw_event_time <- function(i) {
    U <- stats::runif(1)
    target_H <- -log(U)
    f_root <- function(t) cumhaz_i(i, t) - target_H
    br <- .find_bracket(
      f_root,
      lower = 0,
      upper = root_control$t_init,
      upper_max = root_control$t_max,
      expand = root_control$expand,
      max_expand = root_control$max_expand
    )
    if (is.null(br)) {
      return(list(time = time_cens, event = 0L, bracketing_failed = TRUE))
    }
    T_star <- stats::uniroot(f_root, lower = br$lower, upper = br$upper)$root
    if (T_star > time_cens) {
      list(time = time_cens, event = 0L, bracketing_failed = FALSE)
    } else {
      list(time = T_star, event = 1L, bracketing_failed = FALSE)
    }
  }

  # ---- Step 4: Draw event/censoring times from the implied hazard model
  event_draws <- .sim_parallel_lapply(seq_len(n_id), draw_event_time)

  .is_valid_event_draw <- function(x) {
    is.list(x) &&
      !is.null(x$time) && length(x$time) == 1L && is.finite(as.numeric(x$time)) &&
      !is.null(x$event) && length(x$event) == 1L && is.finite(as.numeric(x$event))
  }

  bad_idx <- which(!vapply(event_draws, .is_valid_event_draw, logical(1)))
  if (length(bad_idx) > 0L) {
    cli::cli_warn(c(
      x = "{length(bad_idx)} parallel event-draw job(s) returned malformed output.",
      i = "Recomputing those jobs deterministically in-process."
    ))

    for (idx in bad_idx) {
      set.seed(seed_base + as.integer(idx))
      event_draws[[idx]] <- draw_event_time(idx)
    }

    still_bad <- which(!vapply(event_draws, .is_valid_event_draw, logical(1)))
    if (length(still_bad) > 0L) {
      cli::cli_abort(c(
        x = "Failed to generate valid event draws for indices: {paste(still_bad, collapse = ', ')}.",
        i = "Check custom hazard/transforms for numerical failures."
      ))
    }
  }

  simulated_event_stop <- vapply(
    event_draws,
    function(x) as.numeric(x$time),
    numeric(1)
  ) # observed event or administrative-censoring time returned by numerical inversion
  simulated_event_stop <- .sim_separate_tied_event_endpoints(
    interval_start = 0,
    interval_stop = simulated_event_stop
  ) # positive event-process endpoints used consistently by observation scheduling, returned data and fitting
  dataEvent[[event_time_var]] <- simulated_event_stop
  dataEvent[[event_var]] <- vapply(event_draws, function(x) as.integer(x$event), integer(1))

  # ---- Step 5a: Sample delayed-entry times once per subject.
  # Entry is constrained to [0, stop_i) so each subject remains at risk after
  # entry and before the observed event/censoring stop time.
  entry_time_by_id <- vapply(seq_len(n_id), function(i) {
    stop_i <- as.numeric(dataEvent[[event_time_var]][i])
    entry_upper <- min(left_truncation_max, max(0, stop_i - 1e-8))
    if (entry_upper > 0) stats::runif(1, min = 0, max = entry_upper) else 0
  }, numeric(1))

  # ---- Step 5b: Build longitudinal observation schedule conditional on
  # event/censoring and optional delayed-entry truncation.
  cov_names <- setdiff(colnames(dataEvent), c(id_var, event_time_var, event_var))

  .sim_obs_rows_for_id <- function(i) {
    entry_i <- as.numeric(entry_time_by_id[i])
    obs_upper <- max(1e-8, min(dataEvent[[event_time_var]][i], time_cens))
    rows_i <- vector("list", D)
    for (d in seq_len(D)) {
      t_obs <- times_obs

      if (obs_time_noise_sd > 0) {
        noise <- stats::rnorm(length(t_obs), mean = 0, sd = obs_time_noise_sd)
        t_obs <- sort(pmin(time_cens, pmax(0, t_obs + noise)))
      }
      if (isTRUE(censor_longitudinal_after_event)) {
        t_obs <- t_obs[t_obs <= obs_upper]
      }
      if (isTRUE(truncate_longitudinal_before_entry)) {
        t_obs <- t_obs[t_obs >= entry_i]
      }
      if (!length(t_obs)) {
        rows_i[[d]] <- NULL
        next
      }

      row_df <- data.frame(
        id = rep(dataEvent[[id_var]][i], length(t_obs)),
        marker = factor(rep(marker_levels[d], length(t_obs)), levels = marker_levels),
        time = t_obs,
        stringsAsFactors = FALSE
      )
      names(row_df)[names(row_df) == "id"] <- id_var
      names(row_df)[names(row_df) == "time"] <- time_var
      names(row_df)[names(row_df) == "marker"] <- marker_var

      for (cov_nm in cov_names) {
        row_df[[cov_nm]] <- dataEvent[[cov_nm]][i]
      }
      rows_i[[d]] <- row_df
    }
    rows_i <- Filter(Negate(is.null), rows_i)
    if (!length(rows_i)) {
      row_df <- data.frame(
        id = dataEvent[[id_var]][i],
        marker = factor(marker_levels[1], levels = marker_levels),
        time = obs_upper,
        stringsAsFactors = FALSE
      )
      names(row_df)[names(row_df) == "id"] <- id_var
      names(row_df)[names(row_df) == "time"] <- time_var
      names(row_df)[names(row_df) == "marker"] <- marker_var

      for (cov_nm in cov_names) {
        row_df[[cov_nm]] <- dataEvent[[cov_nm]][i]
      }
      return(row_df)
    }
    do.call(rbind, rows_i)
  }

  obs_rows <- .sim_parallel_lapply(seq_len(n_id), .sim_obs_rows_for_id)
  obs_rows <- Filter(Negate(is.null), obs_rows)
  if (!length(obs_rows)) {
    cli::cli_abort(c(
      x = "No longitudinal observations remained after applying the schedule and event-time censoring.",
      i = "Use a denser {.arg times_obs} grid, reduce {.arg obs_time_noise_sd}, or disable {.arg censor_longitudinal_after_event}."
    ))
  }
  dataLong <- do.call(rbind, obs_rows)
  rownames(dataLong) <- NULL

  # ---- Mean structure from model matrices and sampled random effects
  X_long <- .mm(fixed_template, dataLong)
  Z_id_long <- .sim_rhs_matrix(id_templates, dataLong)
  Z_mk_long <- .sim_rhs_matrix(mk_templates, dataLong)
  Z_idm_long <- .sim_rhs_matrix(idm_templates, dataLong)

  id_index <- match(as.character(dataLong[[id_var]]), as.character(dataEvent[[id_var]]))
  marker_index <- match(as.character(dataLong[[marker_var]]), marker_levels)

  mu_long <- as.numeric(X_long %*% beta_long_internal)
  if (K_id > 0) {
    mu_long <- mu_long + rowSums(Z_id_long * re_id[id_index, , drop = FALSE])
  }
  if (K_mk > 0) {
    mu_long <- mu_long + rowSums(Z_mk_long * re_marker[marker_index, , drop = FALSE])
  }
  if (K_idm > 0) {
    re_idm_flat_scaled <- matrix(aperm(re_idm_scaled, c(2, 1, 3)), nrow = n_id * D, ncol = K_idm)
    idx_flat <- (id_index - 1L) * D + marker_index
    idm_effect <- re_idm_flat_scaled[idx_flat, , drop = FALSE]
    mu_long <- mu_long + rowSums(Z_idm_long * idm_effect)
  }

  # ---- Step: apply inverse-link to the linear predictor per marker
  # This yields response-scale means/probabilities used in sampling.
  mu_linked <- vapply(seq_len(nrow(dataLong)), function(r) {
    .sim_apply_inv_link_bc(mu_long[r], inv_link_specs[[marker_index[r]]])
  }, numeric(1))

  # ---- Distributional parameters (family defaults + formulaDist overrides)
  family_by_row <- family_names[marker_index]
  attr(dataLong, "JoiNMe_family_by_row") <- family_by_row

  sigma_vec <- vapply(family_by_row, function(f) .sim_get_family_param(f, "sigma", 1.0), numeric(1))
  nu_vec <- vapply(family_by_row, function(f) .sim_get_family_param(f, "nu", 4.0), numeric(1))
  phi_vec <- vapply(family_by_row, function(f) .sim_get_family_param(f, "phi", 2.0), numeric(1))
  alpha_vec <- vapply(family_by_row, function(f) .sim_get_family_param(f, "alpha", 0.0), numeric(1))
  kappa_vec <- vapply(family_by_row, function(f) .sim_get_family_param(f, "kappa", 10.0), numeric(1))
  tau_vec <- vapply(family_by_row, function(f) .sim_get_family_param(f, "tau", 0.5), numeric(1))
  trials_vec <- vapply(family_by_row, function(f) .sim_get_family_param(f, "trials", 10L), numeric(1))

  dist_re_effective <- list()
  dist_coef_effective <- list()
  dist_eta_effective <- list()
  dist_design_cols <- list()
  for (param_name in names(dist_formulas)) {
    X_param <- .build_dist_matrix(dist_formulas[[param_name]], dataLong, family_by_row = family_by_row)$X
    beta_param <- .sim_distributional_population(
      simulation_truth_request,
      prior_specification$distributional[[param_name]],
      param_name,
      colnames(X_param),
      marker_levels,
      family_names
    ) # distributional population coefficients drawn or fixed in model-matrix order
    eta_param <- as.numeric(X_param %*% beta_param)

    re_terms_param <- .build_dist_re_terms(dist_formulas[[param_name]], dataLong)
    if (!is.null(re_terms_param$n_re) && re_terms_param$n_re > 0) {
      cfg_dist <- (re_params$dist %||% list())[[param_name]] %||% list()
      dist_re_effective_terms <- vector("list", re_terms_param$n_re)
      for (j in seq_len(re_terms_param$n_re)) {
        cfg_j <- if (!is.null(cfg_dist$terms) && length(cfg_dist$terms) >= j) cfg_dist$terms[[j]] else cfg_dist
        resolved_j <- .sim_resolve_re_block_cfg(
          cfg = cfg_j,
          K = as.integer(re_terms_param$K[j]),
          label = paste0("dist_", param_name, "_re", j)
        )
        dist_re_effective_terms[[j]] <- resolved_j
        b_j <- .sim_draw_re_block(
          n_group = as.integer(re_terms_param$G[j]),
          K = as.integer(re_terms_param$K[j]),
          cfg = cfg_j,
          label = paste0("dist_", param_name, "_re", j),
          resolved = resolved_j
        )
        eta_param <- eta_param + rowSums(re_terms_param$Z[[j]] * b_j[re_terms_param$J[[j]], , drop = FALSE])
      }
      dist_re_effective[[param_name]] <- list(terms = dist_re_effective_terms)
    }

    dist_coef_effective[[param_name]] <- beta_param
    dist_eta_effective[[param_name]] <- eta_param
    dist_design_cols[[param_name]] <- colnames(X_param) %||% character(0)

    if (param_name == "sigma") sigma_vec <- exp(eta_param)
    if (param_name == "nu") nu_vec <- 2 + exp(eta_param)
    if (param_name == "phi") phi_vec <- exp(eta_param)
    if (param_name == "alpha") alpha_vec <- eta_param
    if (param_name == "kappa") kappa_vec <- exp(eta_param)
    if (param_name == "tau") tau_vec <- stats::plogis(eta_param)
  }

  # A fixed family quantile is the final statistical specification for its
  # marker. Apply it after distributional regression so mixed simulations
  # mirror Stan: fixed markers use their declared constants, while remaining
  # skew-Laplace markers retain their family-level or row-specific values.
  fixed_tau_row <- use_tau_fixed[marker_index] == 1L
  tau_vec[fixed_tau_row] <- tau_fixed[marker_index[fixed_tau_row]]

  # ---- Step: draw outcomes by marker-specific family
  y_out <- numeric(nrow(dataLong))
  for (r in seq_len(nrow(dataLong))) {
    fam_name <- family_by_row[r]
    fam_code <- .parse_family(fam_name)
    eta_r <- mu_long[r]
    mu_r <- mu_linked[r]

    if (fam_code == 11L) {
      cp <- family_params[[fam_name]]$cutpoints %||% c(-1, 1)
      y_out[r] <- .sim_sample_ordinal(eta = eta_r, cutpoints = cp)
    } else {
      y_out[r] <- .sample_from_family(
        n = 1,
        mu = mu_r,
        family = fam_code,
        sigma = sigma_vec[r],
        nu = nu_vec[r],
        phi = phi_vec[r],
        kappa = kappa_vec[r],
        tau = tau_vec[r],
        trials = as.integer(round(trials_vec[r])),
        skew = alpha_vec[r]
      )
    }
  }
  dataLong[[y_var]] <- y_out

  # ---- Final formatting and metadata
  dataLong[[marker_var]] <- factor(dataLong[[marker_var]], levels = marker_levels)
  ord_long <- order(dataLong[[id_var]], dataLong[[marker_var]], dataLong[[time_var]])
  dataLong <- dataLong[ord_long, , drop = FALSE]
  family_by_row <- family_by_row[ord_long]
  sigma_vec <- sigma_vec[ord_long]
  nu_vec <- nu_vec[ord_long]
  phi_vec <- phi_vec[ord_long]
  alpha_vec <- alpha_vec[ord_long]
  kappa_vec <- kappa_vec[ord_long]
  tau_vec <- tau_vec[ord_long]
  trials_vec <- trials_vec[ord_long]
  rownames(dataLong) <- NULL

  dist_param_rowwise <- data.frame(
    id = dataLong[[id_var]],
    marker = dataLong[[marker_var]],
    time = dataLong[[time_var]],
    family = family_by_row,
    sigma = sigma_vec,
    nu = nu_vec,
    phi = phi_vec,
    phi_nb = phi_vec,
    alpha = alpha_vec,
    alpha_skew = alpha_vec,
    skew = alpha_vec,
    kappa = kappa_vec,
    tau = tau_vec,
    trials = as.integer(round(trials_vec)),
    stringsAsFactors = FALSE
  )

  distributional_truth <- list(
    family_by_row = family_by_row,
    rowwise = dist_param_rowwise,
    sigma = sigma_vec,
    nu = nu_vec,
    phi = phi_vec,
    phi_nb = phi_vec,
    alpha = alpha_vec,
    alpha_skew = alpha_vec,
    skew = alpha_vec,
    kappa = kappa_vec,
    tau = tau_vec,
    trials = as.integer(round(trials_vec)),
    coef = dist_coef_effective,
    eta = dist_eta_effective,
    design_cols = dist_design_cols,
    family_defaults = family_params,
    use_tau_fixed = stats::setNames(use_tau_fixed, marker_levels),
    tau_fixed = stats::setNames(tau_fixed, marker_levels)
  )

  tmax <- max(dataEvent[[event_time_var]])
  time_scale_generation <- as.numeric(time_scale_internal %||% tmax)

  family_codes_present <- sort(unique(family_codes))
  family_names_present <- vapply(family_codes_present, .family_code_to_name, character(1))

  .sim_family_param_truth <- function(param_name, fallback, include_families = NULL) {
    fam_keep <- family_codes_present[vapply(
      family_codes_present,
      function(fc) param_name %in% .family_distrib_params(fc),
      logical(1)
    )]
    if (!is.null(include_families)) {
      fam_keep <- fam_keep[vapply(fam_keep, function(fc) .family_code_to_name(fc) %in% include_families, logical(1))]
    }
    if (length(fam_keep) == 0L) return(numeric(0))
    fam_names_keep <- vapply(fam_keep, .family_code_to_name, character(1))
    vals <- vapply(fam_names_keep, function(fam) .sim_get_family_param(fam, param_name, fallback), numeric(1))
    names(vals) <- fam_names_keep
    vals
  }

  sigma_family_truth <- .sim_family_param_truth("sigma", 1.0)
  nu_family_truth <- .sim_family_param_truth("nu", 4.0)
  phi_family_truth <- .sim_family_param_truth("phi", 2.0)
  alpha_family_truth <- .sim_family_param_truth("alpha", 0.0)
  kappa_family_truth <- .sim_family_param_truth("kappa", 10.0)
  # The fitted model declares a family-level tau only when at least one
  # skew-Laplace marker has not supplied its own fixed quantile.
  tau_family_truth <- if (any(skew_laplace_marker & use_tau_fixed == 0L)) {
    .sim_family_param_truth("tau", 0.5)
  } else {
    numeric(0)
  }
  trials_family_truth <- if (any(family_names_present == "binomial")) {
    vals <- c(binomial = .sim_get_family_param("binomial", "trials", 10L))
    as.numeric(setNames(vals, names(vals)))
  } else {
    numeric(0)
  }
  cutpoints_ord_truth <- if (any(family_names_present == "cumulative_logit")) {
    as.numeric(family_params$cumulative_logit$cutpoints %||% c(-1, 1))
  } else {
    numeric(0)
  }

  marker_to_sigma_family <- setNames(match(family_names, names(sigma_family_truth), nomatch = 0L), marker_levels)
  marker_to_nu_family <- setNames(match(family_names, names(nu_family_truth), nomatch = 0L), marker_levels)
  marker_to_phi_family <- setNames(match(family_names, names(phi_family_truth), nomatch = 0L), marker_levels)
  marker_to_alpha_family <- setNames(match(family_names, names(alpha_family_truth), nomatch = 0L), marker_levels)
  marker_to_kappa_family <- setNames(match(family_names, names(kappa_family_truth), nomatch = 0L), marker_levels)
  marker_to_tau_family <- setNames(match(family_names, names(tau_family_truth), nomatch = 0L), marker_levels)
  marker_to_tau_family[use_tau_fixed == 1L] <- 0L

  stan_fit_truth <- list(
    beta = beta_long,
    time_scale_generation = as.numeric(time_scale_generation),
    time_scale_observed_max = as.numeric(tmax),
    gamma_w = beta_event,
    bs_gamma_c = baseline_hazard_truth$coefficients,
    shrinkage = shrinkage,
    marker_weights_offset = marker_weights_offset,
    marker_weight_mean = marker_weight_means,
    z_marker_weights = marker_weight_prior_draw$raw,
    marker_weight_standardised = as.numeric(t(standardised_departure_sets)),
    marker_weight_prior_family = prior_specification$marker_weights$family$family,
    marker_weight_df = marker_weight_prior_draw$df %||% numeric(0),
    marker_weight_df_was_fitted = marker_weight_prior_draw$df_was_fitted %||% logical(0),
    marker_weight_prior_raw = marker_weight_prior_draw$raw,
    marker_weight_horseshoe_local = marker_weight_prior_draw$local_scale,
    marker_weight_horseshoe_global = marker_weight_prior_draw$global_scale,
    marker_weight_horseshoe_slab = marker_weight_prior_draw$slab_multiplier,
    marker_weights_eff = marker_weights_effective,
    tau_u = re_id_effective$sd,
    Lcorr_u = re_id_effective$Lcorr,
    Corr_u = re_id_effective$corr,
    Sigma_u = re_id_effective$cov,
    tau_v = re_marker_effective$sd,
    Lcorr_v = re_marker_effective$Lcorr,
    Corr_v = re_marker_effective$corr,
    Sigma_v = re_marker_effective$cov,
    alpha_L = alpha_cov, # fitted packed covariance-regression intercepts
    beta_L_sd = beta_cov_sd, # fitted Q_idm by K_cov_sd standard-deviation slopes
    beta_L_corr = beta_cov_corr, # fitted M_corr by K_cov_corr partial-correlation slopes
    lambda_L = lambda_cov, # fitted non-negative residual loading for each packed latent coordinate
    z_L = z_cov[, seq_len(M_cov), drop = FALSE], # subject-by-coordinate standardised covariance latents
    alpha_cv_total = unname(assoc_coef_scalar[["cv_total"]] %||% NA_real_),
    alpha_cs_total = unname(assoc_coef_scalar[["cs_total"]] %||% NA_real_),
    alpha_cv_mean = unname(assoc_coef_scalar[["cv_mean"]] %||% NA_real_),
    alpha_cs_mean = unname(assoc_coef_scalar[["cs_mean"]] %||% NA_real_),
    alpha_cv_marker = unname(assoc_coef_scalar[["cv_marker"]] %||% NA_real_),
    alpha_cs_marker = unname(assoc_coef_scalar[["cs_marker"]] %||% NA_real_),
    alpha_corr = assoc_coef_corr,
    alpha_vcov = assoc_coef_vcov,
    sigma_family = sigma_family_truth,
    nu_family = nu_family_truth,
    phi_family = phi_family_truth,
    alpha_family = alpha_family_truth,
    kappa_family = kappa_family_truth,
    tau_family = tau_family_truth,
    use_tau_fixed = stats::setNames(use_tau_fixed, marker_levels),
    tau_fixed = stats::setNames(tau_fixed, marker_levels),
    trials_family = trials_family_truth,
    cutpoints_ord = cutpoints_ord_truth,
    marker_to_sigma_family = marker_to_sigma_family,
    marker_to_nu_family = marker_to_nu_family,
    marker_to_phi_family = marker_to_phi_family,
    marker_to_alpha_family = marker_to_alpha_family,
    marker_to_kappa_family = marker_to_kappa_family,
    marker_to_tau_family = marker_to_tau_family
  )

  # Assemble a simulation-truth record in both the public statistical language
  # and the fitted Stan parameter language. The standardised draws are retained
  # because class recovery should be assessed at the level on which the mixture
  # is defined, not by classifying covariance-scaled effects after generation.
  mixture_truth <- mixture_configuration
  if (!is.null(mixture_truth)) {
    mixture_truth$formulaClass <-
      .mixture_specification$formulaClass # original membership formula request
    mixture_truth$standardised_draws <- list(
      subject = z_id,
      marker = z_marker
    ) # latent values to which component locations and scales may have been applied
    covariance_class_type <- intersect(
      c("corr", "vcov"),
      mixture_truth$class_type
    ) # selected public name for the covariance-regression latent matrix
    if (length(covariance_class_type) == 1L) {
      mixture_truth$standardised_draws[[covariance_class_type]] <-
        z_cov[, seq_len(M_cov), drop = FALSE]
    }
    stan_fit_truth$mix_probability <-
      unname(mixture_truth$probability)
    stan_fit_truth$mix_location <-
      unname(mixture_truth$location)
    stan_fit_truth$mix_scale <-
      unname(mixture_truth$scale)
    stan_fit_truth$mix_class_coefficient_subject <-
      unname(mixture_truth$coefficient$subject)
    stan_fit_truth$mix_class_coefficient_marker <-
      unname(mixture_truth$coefficient$marker)
    stan_fit_truth$posterior_class_subject <-
      unname(mixture_truth$allocation$subject)
    stan_fit_truth$posterior_class_marker <-
      unname(mixture_truth$allocation$marker)
  }

  family_truth <- list(
    by_marker = data.frame(
      marker = marker_levels,
      family = family_names,
      family_code = family_codes,
      marker_to_sigma_family = unname(marker_to_sigma_family),
      marker_to_nu_family = unname(marker_to_nu_family),
      marker_to_phi_family = unname(marker_to_phi_family),
      marker_to_alpha_family = unname(marker_to_alpha_family),
      marker_to_kappa_family = unname(marker_to_kappa_family),
      marker_to_tau_family = unname(marker_to_tau_family),
      use_tau_fixed = use_tau_fixed,
      tau_fixed = tau_fixed,
      stringsAsFactors = FALSE
    ),
    shared = stan_fit_truth[c(
      "sigma_family",
      "nu_family",
      "phi_family",
      "alpha_family",
      "kappa_family",
      "tau_family",
      "trials_family",
      "cutpoints_ord"
    )],
    defaults = family_params
  )

  truth <- list(
    beta_long = beta_long,
    beta_event = beta_event,
    alpha_cv_total = if ("cv_total" %in% names(assoc_coef_vec)) assoc_coef_vec[["cv_total"]] else NA_real_,
    alpha_cs_total = if ("cs_total" %in% names(assoc_coef_vec)) assoc_coef_vec[["cs_total"]] else NA_real_,
    alpha_cv_mean = if ("cv_mean" %in% names(assoc_coef_vec)) assoc_coef_vec[["cv_mean"]] else NA_real_,
    alpha_cs_mean = if ("cs_mean" %in% names(assoc_coef_vec)) assoc_coef_vec[["cs_mean"]] else NA_real_,
    alpha_corr = assoc_coef_corr,
    alpha_vcov = assoc_coef_vcov,
    assoc = assoc,
    assoc_coefs = assoc_coef_vec,
    # `marker_weights` remains the convenient public truth alias and denotes
    # the effective weights that actually generated the event process.
    marker_weights = marker_weights_effective,
    marker_weights_raw = marker_weights_effective,
    marker_weights_offset = marker_weights_offset,
    marker_weight_mean = marker_weight_means,
    marker_weight_mean_by_term = marker_weight_mean_by_term,
    marker_weights_latent = marker_weights_latent,
    marker_weights_eff = marker_weights_eff,
    marker_weights_by_term = marker_weights_eff,
    marker_weights_offset_by_term = marker_weights_offset_by_term,
    marker_weights_latent_by_term = marker_weights_latent_by_term,
    z_marker_weight_sets = z_marker_weight_sets,
    marker_weight_standardised_sets = standardised_departure_sets,
    marker_weight_prior_family = prior_specification$marker_weights$family$family,
    marker_weight_df = marker_weight_prior_draw$df %||% numeric(0),
    marker_weight_df_was_fitted = marker_weight_prior_draw$df_was_fitted %||% logical(0),
    marker_weight_prior_raw = marker_weight_prior_draw$raw,
    marker_weight_horseshoe_local = marker_weight_prior_draw$local_scale,
    marker_weight_horseshoe_global = marker_weight_prior_draw$global_scale,
    marker_weight_horseshoe_slab = marker_weight_prior_draw$slab_multiplier,
    shrinkage = shrinkage,
    shrinkage_distribution = c(
      `0` = "student_t(6, 0, 1)",
      `1` = "double_exponential(0, 1)",
      `2` = "normal(0, 1)"
    )[[as.character(shrinkage)]],
    link_names = link_names,
    quadrature_nodes = as.integer(gk_spec$n_gk),
    gk_nodes = gk_spec$nodes,
    gk_weights = gk_spec$weights,
    gk_rule = gk_spec$rule,
    transforms = transforms,
    declaration = simulation_truth_request,
    basehaz = h0_fn,
    baseline_hazard = baseline_hazard_truth,
    formulaVCov = formulaVCov,
    formulaBasehaz = formulaBasehaz,
    beta_basehaz = baseline_hazard_truth$coefficients,
    family = family_truth,
    distributional_params = distributional_truth,
    stan_fit = stan_fit_truth,
    dist_re_params = re_params$dist %||% list(),
    dist_re = dist_re_effective %||% list(),
    re_params = re_params,
    re_structure = list(
      id = re_id_effective,
      marker = re_marker_effective,
      id_marker_cov = list(
        standardised = re_idm_effective,
        sd = list(
          alpha = alpha_cov_sd,
          beta = beta_cov_sd,
          lambda = lambda_cov_sd,
          z = z_cov[, diagonal_positions, drop = FALSE],
          formula = formulaVCov$sd
        ),
        corr = list(
          alpha = alpha_cov_corr,
          beta = beta_cov_corr,
          lambda = lambda_cov_corr,
          z = z_cov[, correlation_positions, drop = FALSE],
          formula = formulaVCov$corr
        ),
        alpha = alpha_cov,
        beta = list(sd = beta_cov_sd, corr = beta_cov_corr),
        lambda = lambda_cov,
        z = z_cov[, seq_len(M_cov), drop = FALSE],
        lambda_sign = lambda_cov_sign,
        diag_link = diag_link_cov
      )
    ),
    re_draws = list(
      id = re_id_public,
      marker = re_marker_public,
      id_marker_cov = re_idm_public,
      id_marker_cov_standardised = re_idm,
      id_marker_cov_scaled = re_idm_scaled
    ),
    re_draws_likelihood = list(
      id = re_id,
      marker = re_marker,
      id_marker_cov = re_idm_scaled
    ),
    L_i = L_i,
    id_marker_cov_effective = list(
      standardised = re_idm_effective,
      sd = list(
        alpha = alpha_cov_sd,
        beta = beta_cov_sd,
        lambda = lambda_cov_sd,
        z = z_cov[, diagonal_positions, drop = FALSE],
        formula = formulaVCov$sd
      ),
      corr = list(
        alpha = alpha_cov_corr,
        beta = beta_cov_corr,
        lambda = lambda_cov_corr,
        z = z_cov[, correlation_positions, drop = FALSE],
        formula = formulaVCov$corr
      ),
      alpha = alpha_cov,
      beta = list(sd = beta_cov_sd, corr = beta_cov_corr),
      lambda = lambda_cov,
      z = z_cov[, seq_len(M_cov), drop = FALSE],
      lambda_sign = lambda_cov_sign,
      diag_link = diag_link_cov
    ),
    mixture = mixture_truth,
    formulaLong = formulaLong,
    formulaEvent = formulaEvent,
    formulaDist = dist_formulas
  )

  marker_info <- list(
    names = marker_levels,
    families = families,
    family_codes = family_codes,
    use_tau_fixed = stats::setNames(use_tau_fixed, marker_levels),
    tau_fixed = stats::setNames(tau_fixed, marker_levels)
  )

  helpers <- list(
    cv_total = function(i, t) assoc_components(i, t)$cv_total,
    cv_mean = function(i, t) assoc_components(i, t)$cv_mean,
    cv_marker = function(i, t) assoc_components(i, t)$cv_marker,
    cs_total = function(i, t) assoc_components(i, t)$cs_total,
    cs_mean = function(i, t) assoc_components(i, t)$cs_mean,
    cs_marker = function(i, t) assoc_components(i, t)$cs_marker,
    corr = function(i, t) assoc_components(i, t)$corr,
    vcov = function(i, t) assoc_components(i, t)$vcov,
    assoc_components_raw = function(i, t) assoc_components(i, t)$raw,
    assoc_components = assoc_components,
    baseline_hazard = h0_fn,
    hazard = hazard_i,
    cumhaz = cumhaz_i,
    survival_prob = function(i, t) exp(-cumhaz_i(i, t))
  )

  # ---- Step 5c: Build event rows automatically from the simulated risk-process
  # structure. Counting-process rows are emitted whenever split layout is
  # required by the data-generating mechanism.
  dataEvent_public <- dataEvent
  if (event_layout == "split") {
    interval_rows <- lapply(seq_len(n_id), function(i) {
      stop_i <- as.numeric(dataEvent[[event_time_var]][i])
      status_i <- as.integer(dataEvent[[event_var]][i])
      start_i <- as.numeric(entry_time_by_id[i])

      tv_internal <- numeric(0)
      if (length(tv_cov_specs) > 0L) {
        tv_internal <- unique(unlist(lapply(tv_cov_specs, function(spec_all) {
          breaks_i <- as.numeric(spec_all[[i]]$breaks %||% numeric(0))
          breaks_i[breaks_i > start_i & breaks_i < stop_i]
        }), use.names = FALSE))
      }
      internal <- times_obs[times_obs > start_i & times_obs < stop_i]
      breaks <- c(start_i, internal, tv_internal, stop_i)
      breaks <- sort(unique(as.numeric(breaks)))
      if (length(breaks) < 2L) {
        breaks <- c(start_i, stop_i)
      }
      int_start <- breaks[-length(breaks)]
      int_stop <- breaks[-1]
      int_status <- rep.int(0L, length(int_start))
      if (status_i == 1L && length(int_status) > 0L) {
        int_status[length(int_status)] <- 1L
      }

      out <- dataEvent[rep(i, length(int_start)), , drop = FALSE]
      if (length(tv_cov_specs) > 0L) {
        for (nm in names(tv_cov_specs)) {
          out[[nm]] <- vapply(int_start, function(tt) {
            .sim_time_varyring_value(tv_cov_specs[[nm]][[i]], tt)
          }, numeric(1))
        }
      }
      out$time_start <- as.numeric(int_start)
      separated_interval_stop <- .sim_separate_tied_event_endpoints(
        interval_start = out$time_start,
        interval_stop = as.numeric(int_stop)
      ) # exact numerical ties are separated without concealing a stop lying genuinely below its start
      out$time_stop <- separated_interval_stop
      out[[event_time_var]] <- separated_interval_stop
      out[[event_var]] <- as.integer(int_status)
      out
    })
    dataEvent_public <- do.call(rbind, interval_rows)
    rownames(dataEvent_public) <- NULL
  } else {
    if (length(tv_cov_specs) > 0L) {
      for (nm in names(tv_cov_specs)) {
        dataEvent_public[[nm]] <- vapply(seq_len(nrow(dataEvent_public)), function(i) {
          .sim_time_varyring_value(tv_cov_specs[[nm]][[i]], 0)
        }, numeric(1))
      }
    }
    dataEvent_public$time_start <- as.numeric(entry_time_by_id)
    separated_subject_stop <- .sim_separate_tied_event_endpoints(
      interval_start = dataEvent_public$time_start,
      interval_stop = dataEvent_public[[event_time_var]]
    ) # subject-level endpoint used by both the public event formula and its auxiliary interval columns
    dataEvent_public[[event_time_var]] <- separated_subject_stop
    dataEvent_public$time_stop <- separated_subject_stop
  }

# ---- Retain a directly reusable fitting specification.
  #
  # This record is assembled only after the public event data have reached
  # their final subject-level or counting-process layout. Every name below is
  # therefore an argument understood by the corresponding fitting entry
  # point. Generating quantities such as beta_long and the covariance
  # coefficients remain in `truth`; they are estimands, not fitting inputs.
  recovery_arguments <- list(
    formulaLong = formulaLong,
    dataLong = dataLong,
    formulaEvent = formulaEvent,
    dataEvent = dataEvent_public,
    formulaVCov = formulaVCov,
    formulaDist = dist_formulas,
    families = families,
    transforms = transforms,
    priors = prior_specification,
    control = list(
      quadrature_nodes = as.integer(gk_spec$n_gk),
      vcov_diag_link = diag_link_cov
    ),
    assoc = assoc,
    shrinkage = shrinkage,
    eps_fd = eps_cs,
    id_var = id_var,
    marker_var = marker_var,
    time_var = time_var
  ) # complete common fitting syntax for a simulation-recovery analysis
  recovery_entry_point <- "joinme" # ordinary fitting function paired with simulate_joinme()
  if (!is.null(.mixture_specification)) {
    recovery_entry_point <- "joinme_mix"
    recovery_arguments$n_classes <- .mixture_specification$n_classes
    recovery_arguments$formulaClass <- .mixture_specification$formulaClass
    recovery_arguments$class_type <- .mixture_specification$class_type
    recovery_arguments$class_dimensions <- .mixture_specification$class_dimensions
    recovery_arguments$class_ordering <- .mixture_specification$class_ordering
  }
  truth$recovery <- list(
    entry_point = recovery_entry_point,
    arguments = recovery_arguments
  ) # call with do.call(get(entry_point), arguments) to rebuild the fitted specification exactly

  if (longitudinal_only) {
    dataEvent_public <- NULL # no event data belong to the public longitudinal-only simulation
    truth$formulaEvent <- NULL # scientific truth records the absence of a survival submodel
    truth$assoc <- character(0) # association parameters are absent rather than fixed at zero
    truth$recovery$arguments$formulaEvent <- NULL # recovery invokes the longitudinal-only fitting route
    truth$recovery$arguments$dataEvent <- NULL # no neutral scaffold is exposed to fitting
    truth$recovery$arguments$assoc <- character(0) # fitted association design has zero columns
  }

  list(
    dataLong = dataLong,
    dataEvent = dataEvent_public,
    truth = truth,
    marker_info = marker_info,
    helpers = helpers,
    tmax = tmax
  )
}
