# Transform Specification Helper
#
# Provides utilities to build and manage transformation specifications
# for the JoiNMe Stan model's association term.
# File overview:
# - Map user transform specs into Stan-friendly bytecode/spline data.
# - Validate functional bytecode and spline shapes before sampling.
#' Canonicalise transform type labels
#'
#' @description
#' Map user-facing transform aliases onto the internal canonical transform names
#' used throughout standata construction, simulation, plotting, and summaries.
#'
#' @param type Character scalar naming a transform type. Accepted canonical or
#'   alias values include `"identity"`, `"functional"`, `"ispline"`,
#'   `"ispline_penalised"`, `"ispline_penalized"`, `"ispline_expit"`,
#'   `"ispline_expit_penalised"`, `"ispline_expit_penalized"`,
#'   `"ispline_exp_penalised"`, `"ispline_exp_penalized"`, and `"pwlin"`.
#'
#' @return Character scalar giving the canonical transform type.
#' @keywords internal
#' @noRd
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

#' Detect expit-domain spline transform declarations
#'
#' @description
#' Return `TRUE` when a transform specification corresponds to an expit-domain
#' I-spline family, meaning the runtime feature is mapped through `plogis()`
#' before spline evaluation.
#'
#' @param spec_or_type Either a full transform specification list or a character
#'   transform type. Accepted positive cases are `"ispline_expit"` and
#'   `"ispline_expit_penalised"` plus their canonicalised aliases.
#'
#' @return Logical scalar.
#' @keywords internal
#' @noRd
.transform_uses_expit_input <- function(spec_or_type) {
  type <- if (is.list(spec_or_type)) spec_or_type$type %||% "identity" else spec_or_type
  .canonicalise_transform_type(type) %in% c("ispline_expit", "ispline_expit_penalised")
}

#' Detect I-spline transform declarations
#'
#' @description
#' Return `TRUE` when a transform declaration belongs to the I-spline family,
#' including penalised and expit-domain variants.
#'
#' @param spec_or_type Either a full transform specification list or a character
#'   transform type. Accepted positive cases are `"ispline"`,
#'   `"ispline_penalised"`, `"pmonospline"`, `"pmono"`,
#'   `"ispline_expit"`, and `"ispline_expit_penalised"` after
#'   canonicalisation.
#'
#' @return Logical scalar.
#' @keywords internal
#' @noRd
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

#' Map runtime transform input onto its spline domain
#'
#' @description
#' Apply any transform-family-specific preprocessing required before spline
#' evaluation. At present this means applying `plogis()` for expit-domain
#' I-spline declarations and leaving all other transform families unchanged.
#'
#' @param x Numeric vector of raw association-feature values.
#' @param spec_or_type Either a transform specification list or a character
#'   transform type.
#'
#' @return Numeric vector on the domain expected by the configured transform.
#' @keywords internal
#' @noRd
.transform_input_for_spec <- function(x, spec_or_type) {
  x <- as.numeric(x)
  if (.transform_uses_expit_input(spec_or_type)) {
    return(stats::plogis(x))
  }
  x
}

#' Validate expit-domain spline inputs
#'
#' @description
#' Check that user-supplied expit-domain training values or knots are finite and
#' live in `[0, 1]`, which is the valid basis domain for
#' `ispline_expit` and `ispline_expit_penalised` declarations.
#'
#' @param values Numeric vector to validate.
#' @param arg_name Character scalar naming the argument being validated.
#'
#' @return Numeric vector identical to `values` after validation.
#' @keywords internal
#' @noRd
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
      i = "Supply values in [0, 1]; {.fn JoiNMe} applies {.fn plogis} only to the raw association feature being transformed."
    ))
  }
  values
}

#' Resolve monotone spline direction
#'
#' @description
#' Validate and canonicalise the requested monotone direction for penalised
#' I-spline transforms.
#'
#' @param direction Optional direction request. Accepted values are
#'   `"increasing"`, `"decreasing"`, `1`, and `-1`. When `NULL`, the
#'   supplied `default` is used.
#' @param default Default direction used when `direction` is `NULL`.
#'
#' @return Integer direction code: `1L` for increasing, `-1L` for decreasing.
#' @keywords internal
#' @noRd
.get_monotone_direction <- function(direction = NULL, default = "increasing") {
  direction <- direction %||% default
  if (is.numeric(direction) && length(direction) == 1L && is.finite(direction)) {
    if (direction > 0) return(1L)
    if (direction < 0) return(-1L)
  }

  direction_chr <- tolower(trimws(as.character(direction)[1]))
  if (direction_chr %in% c("increasing", "increase", "inc", "+1", "1")) {
    return(1L)
  }
  if (direction_chr %in% c("decreasing", "decrease", "dec", "-1")) {
    return(-1L)
  }

  cli::cli_abort(c(
    x = "{.arg direction} must be either {.val increasing} or {.val decreasing}.",
    i = "Numeric aliases {.val 1} and {.val -1} are also accepted."
  ))
}

#' Infer monotone spline direction from training pairs
#'
#' @description
#' Infer whether the target transform is increasing or decreasing by comparing
#' the first and last fitted targets after sorting by `x`.
#'
#' @param x Numeric input locations.
#' @param y Numeric target values.
#'
#' @return Integer direction code: `1L` for increasing, `-1L` for decreasing.
#' @keywords internal
#' @noRd
.infer_monotone_direction <- function(x, y) {
  ord <- order(x)
  x <- as.numeric(x)[ord]
  y <- as.numeric(y)[ord]
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(y) < 2L) {
    return(1L)
  }
  if (y[length(y)] < y[1]) {
    return(-1L)
  }
  1L
}

#' Format monotone spline direction labels
#'
#' @description
#' Convert an integer monotone direction code into the user-facing labels used
#' in summaries and printed transform formulas.
#'
#' @param direction Integer direction code.
#'
#' @return Character scalar.
#' @keywords internal
#' @noRd
.monotone_direction_label <- function(direction) {
  if (.get_monotone_direction(direction, default = 1L) < 0L) {
    "decreasing"
  } else {
    "increasing"
  }
}

#' Build a Stan-estimated ordered piecewise-linear specification
#'
#' @description
#' Convert a user-facing piecewise-linear declaration into the fixed knot
#' locations and estimation metadata required by Stan.  The knot ordinates are
#' deliberately not taken from an observed outcome vector.  Instead, the first
#' ordinate is anchored at zero and the remaining ordinates are reconstructed
#' from the cumulative sum of a simplex with `K - 1` elements.  The final
#' ordinate is consequently `1` for an increasing curve and `-1` for a
#' decreasing curve.  The association coefficient supplies the posterior
#' log-hazard span, while the simplex supplies the relative allocation of that
#' span between adjacent knots.
#'
#' Anchoring the first ordinate is statistically necessary: adding a common
#' constant to every ordinate would be indistinguishable from changing the log
#' baseline hazard after multiplication by the association coefficient. Thus,
#' fitted transform ordinates are anchored relative association ordinates;
#' their products with the association coefficient are relative log-hazard
#' contributions rather than separate baseline-hazard intercepts.
#'
#' The modern API uses `knots` (or the synonymous `cutpoints`) and `direction`.
#' The former `x`/`y` API remains accepted for fitted models.  In that legacy
#' form, `x` supplies the knots and `y` is used only to infer the direction when
#' `direction` is omitted; it never fixes the fitted association curve.  The
#' simulation implementation intentionally retains its former fixed `x`/`y`
#' interpolation semantics.
#'
#' @param spec Named piecewise-linear transform specification.  Recognised
#'   fields are `knots`, `cutpoints`, `x`, `y`, `direction`, and `lambda`.
#'   `lambda` is an optional non-negative second-difference penalty applied to
#'   the reconstructed ordinates; its default is zero because the simplex
#'   already regularises the finite-dimensional shape.
#'
#' @return A canonical transform specification containing sorted knots, a
#'   direction code, an anchored coefficient template, and Stan estimation
#'   metadata.
#' @keywords internal
#' @noRd
.make_stan_pwlin_transform <- function(spec) {
  knots <- spec$knots %||% spec$cutpoints %||% spec$x
  if (is.null(knots)) {
    cli::cli_abort(c(
      x = "Fitted piecewise-linear associations require {.arg knots} (or legacy {.arg x}).",
      i = "For example, use list(type = 'pwlin', knots = c(-2, -1, 0, 1, 2), direction = 'increasing')."
    ))
  }
  knots <- as.numeric(knots)
  if (length(knots) < 2L || any(!is.finite(knots)) || any(diff(knots) <= 0)) {
    cli::cli_abort(c(
      x = "{.arg knots} must be a strictly increasing finite vector with at least two values.",
      i = "Each adjacent pair defines one linear association interval."
    ))
  }

  legacy_y <- spec$y
  if (!is.null(legacy_y)) {
    legacy_y <- as.numeric(legacy_y)
    if (length(legacy_y) != length(knots) || any(!is.finite(legacy_y))) {
      cli::cli_abort(c(
        x = "Legacy {.arg y} must be finite and have one value per knot.",
        i = "The values no longer fix a fitted curve; omit y and state direction explicitly for new analyses."
      ))
    }
  }

  inferred_direction <- if (is.null(legacy_y)) "increasing" else {
    .infer_monotone_direction(knots, legacy_y)
  }
  direction <- .get_monotone_direction(spec$direction, default = inferred_direction)

  lambda <- spec$lambda %||% 0
  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda < 0) {
    cli::cli_abort(c(
      x = "{.arg lambda} must be a non-negative numeric scalar.",
      i = "Use zero to omit the optional second-difference penalty."
    ))
  }
  lambda <- as.numeric(lambda)

  n_knots <- length(knots)
  list(
    type = "pwlin",
    knots = knots,
    x = knots,
    legacy_y = legacy_y,
    direction = .monotone_direction_label(direction),
    spline_direction = direction,
    coeff = seq(0, direction, length.out = n_knots),
    degree = 1L,
    estimate_spline = 1L,
    lambda_spline = lambda,
    n_free_spline = as.integer(n_knots - 1L)
  )
}

#' Transform Specification Helper
#'
#' Provides utilities to build and manage transformation specifications
#' for the JoiNMe Stan model's association term.
#'
#' @description
#' Transform can be specified in four modes:
#'   - Mode 0: Identity (no transformation)
#'   - Mode 1: Functional (arbitrary nested functions via parser)
#'   - Mode 2: I-spline basis (monotonic spline)
#'   - Mode 2 (penalised): I-spline coefficients estimated with smoothness penalty
#'   - Mode 4: I-spline basis evaluated on `expit(x)` for a bounded, numerically
#'     stable spline input domain
#'   - Mode 3: Increasing ordered piecewise-linear association
#'   - Mode 7: Decreasing ordered piecewise-linear association
#'
#' Each transformation term (CV_total, CS_total, CV_mean, CS_mean, CV_marker,
#' CS_marker, corr, vcov) can have its own independent specification, enabling
#' flexibility in model building. For covariance-style channels, `corr`
#' transformations act on off-diagonal entries of the subject-specific
#' Cholesky-correlation factor `K`, while `vcov` transformations act on those
#' same off-diagonal `K` entries together with the subject-specific standard
#' deviations.
#'
#' @param transform_list List with elements cv_total, cs_total, cv_mean, cs_mean,
#'   cv_marker, cs_marker, corr, vcov, each specifying a transformation. See details.
#' @param default_mode Default transformation mode if not specified.
#' @param n_corr_components Optional non-negative number of correlation
#'   components. It determines the row count of component-specific coefficient
#'   matrices.
#' @param n_vcov_components Optional non-negative number of covariance
#'   components. It determines the row count of component-specific coefficient
#'   matrices.
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
#'       * if y is supplied, JoiNMe fits the monotone spline to the training pairs (x, y) in R
#'         using the same anchored endpoint convention as the Stan-estimated path
#'         (first coefficient = 0, last coefficient = 1)
#'       * if y is omitted, Stan estimates the monotone spline coefficients directly and lambda controls smoothness
#'     - ispline_expit / ispline_expit_penalised: identical to the I-spline
#'       variants above, except the raw association feature is first mapped
#'       through `plogis(x)` and the spline basis is then evaluated on that
#'       bounded expit-scale input. User-supplied training `x` values and
#'       explicit `knots` for these transform types must therefore already be
#'       specified on the expit scale in `[0, 1]`.
#'     - pwlin: knots (or cutpoints/x) and direction; legacy y is accepted but
#'       does not determine the fitted ordinates
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
#'     `list(type = "pwlin", knots = c(-2, -1, 0, 1, 2), direction = "increasing")`
#'
#' Default values used internally:
#'   - omitted term -> identity mode (`default_mode = 0`),
#'   - `ispline`: `degree = 3` if omitted,
#'   - `ispline_penalised`: `n_knots = 6`, `degree = 3`, `lambda = 1` if omitted,
#'   - `pwlin`: `direction = "increasing"` if omitted,
#'   - `functional`: no additional defaults beyond its required fields.
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
#'   - estimate_iota_intercept_* and estimate_iota_slope_* (fit-only functional affine-shift counts)
#'   - functional_iota_intercept_idx_* and functional_iota_slope_idx_* (per-op affine-shift indices)
#'   - knots_* and coeff_* and spline_degree_* (if mode 2, 3, or 4)
#'
#' @section Usage:
#' Use `build_standata_transforms()` to construct the data list entries,
#' and `validate_transforms()` to verify consistency before sampling.
build_standata_transforms <- function(
  transform_list = NULL,
  default_mode = 0,  # 0 = identity
  n_corr_components = NULL,
  n_vcov_components = NULL
) {
  # Output list to be merged into the main standata object
  standata <- list()
  # Covariance-style channels can carry one transform per component. The same
  # user-facing spec is expanded across those components here so Stan receives a
  # rectangular coefficient matrix even when the user wrote only one spec.
  component_counts <- list(
    corr = as.integer(n_corr_components %||% 1L),
    vcov = as.integer(n_vcov_components %||% 1L)
  )
  
  term_defaults <- list(
    cv_tot = list(mode_suffix = "cv_tot", short_suffix = "cv"),
    cs_tot = list(mode_suffix = "cs_tot", short_suffix = "cs"),
    cv_mean = list(mode_suffix = "cv_mean", short_suffix = "cv_mean"),
    cs_mean = list(mode_suffix = "cs_mean", short_suffix = "cs_mean"),
    cv_marker = list(mode_suffix = "cv_marker", short_suffix = "cv_marker"),
    cs_marker = list(mode_suffix = "cs_marker", short_suffix = "cs_marker"),
    corr = list(mode_suffix = "corr", short_suffix = "corr"),
    vcov = list(mode_suffix = "vcov", short_suffix = "vcov")
  )
  
  # Default empty specs for all terms (identity transformation)
  for (term in names(term_defaults)) {
    mode_suffix <- term_defaults[[term]]$mode_suffix
    short_suffix <- term_defaults[[term]]$short_suffix
    n_components <- if (term %in% c("corr", "vcov")) max(0L, component_counts[[term]]) else 1L
    standata[[paste0("tf_mode_", mode_suffix)]] <- default_mode
    standata[[paste0("n_functional_ops_", short_suffix)]] <- 0
    standata[[paste0("functional_ops_", short_suffix)]] <- integer()
    standata[[paste0("n_const_", short_suffix)]] <- 0
    standata[[paste0("const_data_", short_suffix)]] <- numeric()
    standata[[paste0("functional_iota_intercept_idx_", short_suffix)]] <- integer()
    standata[[paste0("functional_iota_slope_idx_", short_suffix)]] <- integer()
    standata[[paste0("n_knots_", short_suffix)]] <- 0
    standata[[paste0("knots_", short_suffix)]] <- numeric()
    standata[[paste0("n_coeff_", short_suffix)]] <- 0
    standata[[paste0("coeff_", short_suffix)]] <- if (term %in% c("corr", "vcov")) {
      matrix(numeric(0), nrow = n_components, ncol = 0L)
    } else {
      numeric()
    }
    standata[[paste0("spline_degree_", short_suffix)]] <- 1L
    standata[[paste0("estimate_spline_", short_suffix)]] <- 0L
    standata[[paste0("lambda_spline_", short_suffix)]] <- 0.0
    standata[[paste0("n_free_spline_", short_suffix)]] <- 0L
    standata[[paste0("estimate_iota_intercept_", short_suffix)]] <- 0L
    standata[[paste0("estimate_iota_slope_", short_suffix)]] <- 0L
  }
  
  # Override with user specifications
  # Algorithm:
  #   Step 1: normalise the user-facing term name and transform type.
  #   Step 2: convert the spec into a Stan mode code plus the required inputs.
  #   Step 3: for corr/vcov, replicate the coefficient matrix across
  #           components so each component can later evolve independently if the
  #           spline is estimated in Stan.
  #   Step 4: store the serialised result under the common standata schema.
  if (!is.null(transform_list)) {
    for (term_name in names(transform_list)) {
      spec <- transform_list[[term_name]]
      term_info <- .get_transform_term(term_name)
      mode_suffix <- term_info$mode_suffix
      short_suffix <- term_info$short_suffix
      n_components <- if (term_name %in% c("corr", "vcov")) max(0L, component_counts[[term_name]]) else 1L
      spec$type <- .canonicalise_transform_type(spec$type)
      
      if (is.null(spec) || spec$type == "identity") {
        standata[[paste0("tf_mode_", mode_suffix)]] <- 0
      } else if (spec$type == "functional") {
        bc <- parse_transform_expr(spec$expr, iota_nodes = .transform_iota_nodes(spec))
        ops <- bc$bytecode
        standata[[paste0("tf_mode_", mode_suffix)]] <- 1
        standata[[paste0("n_functional_ops_", short_suffix)]] <- length(ops)
        standata[[paste0("functional_ops_", short_suffix)]] <- ops
        standata[[paste0("n_const_", short_suffix)]] <- bc$n_const
        standata[[paste0("const_data_", short_suffix)]] <- bc$const_data
        standata[[paste0("functional_iota_intercept_idx_", short_suffix)]] <- as.integer(bc$op_iota_intercept_idx %||% rep(0L, length(ops)))
        standata[[paste0("functional_iota_slope_idx_", short_suffix)]] <- as.integer(bc$op_iota_slope_idx %||% rep(0L, length(ops)))
        standata[[paste0("estimate_iota_intercept_", short_suffix)]] <- as.integer(spec$n_iota_intercept %||% bc$n_iota_intercept %||% 0L)
        standata[[paste0("estimate_iota_slope_", short_suffix)]] <- as.integer(spec$n_iota_slope %||% bc$n_iota_slope %||% 0L)
      } else if (.is_ispline_transform_type(spec$type)) {
        if (spec$type %in% c("ispline_penalised", "pmonospline", "pmono", "ispline_expit_penalised")) {
          # Two supported semantics:
          # 1) plug-in mode: user supplies x/y and we fit coefficients in R,
          # 2) Stan-estimated mode: user omits y and Stan estimates spline shape
          #    with a smoothness penalty controlled by lambda.
          # The branch taken here determines whether alpha recovery should be
          # interpreted as "given a fixed transform" (plug-in mode) or as part
          # of a jointly estimated transform-plus-association system.
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
        spline_direction <- .get_monotone_direction(spec$direction %||% spec$spline_direction %||% 1L, default = 1L)
        standata[[paste0("tf_mode_", mode_suffix)]] <- if (.transform_uses_expit_input(spec$type)) {
          if (spline_direction < 0L) 6L else 4L
        } else {
          if (spline_direction < 0L) 5L else 2L
        }
        standata[[paste0("n_knots_", short_suffix)]] <- length(spec$knots)
        standata[[paste0("knots_", short_suffix)]] <- spec$knots
        standata[[paste0("n_coeff_", short_suffix)]] <- length(spec$coeff)
        # For corr/vcov, the same initial spline shape is copied to every
        # component. Stan later treats each component row as its own curve when
        # estimate_spline_* = 1, but this common initialisation keeps the R-side
        # layout concise.
        standata[[paste0("coeff_", short_suffix)]] <- if (term_name %in% c("corr", "vcov")) {
          matrix(
            rep(as.numeric(spec$coeff), times = n_components),
            nrow = n_components,
            ncol = length(spec$coeff),
            byrow = TRUE
          )
        } else {
          spec$coeff
        }
        standata[[paste0("spline_degree_", short_suffix)]] <- spec$degree %||% 3L
        standata[[paste0("estimate_spline_", short_suffix)]] <- as.integer(spec$estimate_spline %||% 0L)
        standata[[paste0("lambda_spline_", short_suffix)]] <- as.numeric(spec$lambda_spline %||% 0.0)
        standata[[paste0("n_free_spline_", short_suffix)]] <- as.integer(spec$n_free_spline %||% 0L)
      } else if (spec$type == "pwlin") {
        # A fitted piecewise-linear curve learns its ordered ordinates in Stan.
        # User-supplied y values from the former API are never treated as
        # observed responses.  They are retained only for direction inference
        # by .make_stan_pwlin_transform().
        spec <- .make_stan_pwlin_transform(spec)
        standata[[paste0("tf_mode_", mode_suffix)]] <- if (spec$spline_direction < 0L) 7 else 3
        standata[[paste0("n_knots_", short_suffix)]] <- length(spec$knots)
        standata[[paste0("knots_", short_suffix)]] <- spec$knots
        standata[[paste0("n_coeff_", short_suffix)]] <- length(spec$coeff)
        standata[[paste0("coeff_", short_suffix)]] <- if (term_name %in% c("corr", "vcov")) {
          matrix(
            rep(as.numeric(spec$coeff), times = n_components),
            nrow = n_components,
            ncol = length(spec$coeff),
            byrow = TRUE
          )
        } else {
          spec$coeff
        }
        standata[[paste0("spline_degree_", short_suffix)]] <- 1L
        standata[[paste0("estimate_spline_", short_suffix)]] <- 1L
        standata[[paste0("lambda_spline_", short_suffix)]] <- spec$lambda_spline
        standata[[paste0("n_free_spline_", short_suffix)]] <- spec$n_free_spline
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

#' @keywords internal
.get_transform_term <- function(term_name) {
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
    vcov = list(mode_suffix = "vcov", short_suffix = "vcov"),
    cli::cli_abort(c(
      x = "Unknown transform term: {term_name}.",
      i = "Use cv_total, cs_total, cv_mean, cs_mean, cv_marker, cs_marker, corr, or vcov."
    ))
  )
}

#' Validate Transformation Specifications
#'
#' Check functional bytecode integrity and spline coefficient shapes before sampling.
#'
#' @param standata Prepared standata list with transformation entries
#' @return TRUE if valid; stops with error message if invalid
#'
#' @keywords internal
#' @noRd
validate_transforms <- function(standata) {
  term_map <- list(
    cv = list(mode_suffix = "cv_tot", short_suffix = "cv"),
    cs = list(mode_suffix = "cs_tot", short_suffix = "cs"),
    cv_mean = list(mode_suffix = "cv_mean", short_suffix = "cv_mean"),
    cs_mean = list(mode_suffix = "cs_mean", short_suffix = "cs_mean"),
    cv_marker = list(mode_suffix = "cv_marker", short_suffix = "cv_marker"),
    cs_marker = list(mode_suffix = "cs_marker", short_suffix = "cs_marker"),
    corr = list(mode_suffix = "corr", short_suffix = "corr"),
    vcov = list(mode_suffix = "vcov", short_suffix = "vcov")
  )
  for (term_name in names(term_map)) {
    mode_suffix <- term_map[[term_name]]$mode_suffix
    short_suffix <- term_map[[term_name]]$short_suffix
    mode <- standata[[paste0("tf_mode_", mode_suffix)]]
    
    if (mode == 1) {
      # Functional bytecode validation
      functional_ops <- standata[[paste0("functional_ops_", short_suffix)]]
      const_data <- standata[[paste0("const_data_", short_suffix)]]
      verify_bytecode(functional_ops, const_data)
    } else if (mode %in% c(2, 3, 4, 5, 6, 7)) {
      # Spline validation
      knots <- standata[[paste0("knots_", short_suffix)]]
      coeff <- standata[[paste0("coeff_", short_suffix)]]
      coeff_len <- if (is.matrix(coeff)) ncol(coeff) else length(coeff)
      
      if (length(knots) == 0) {
        cli::cli_abort(c(x = "Spline knots empty for {short_suffix}.", i = "Provide at least two knots."))
      }
      if (coeff_len == 0) {
        cli::cli_abort(c(x = "Spline coefficients empty for {short_suffix}.", i = "Provide coefficients matching the knots."))
      }
      
      if (any(diff(knots) <= 0)) {
        cli::cli_abort(c(x = "Knots not strictly increasing for {short_suffix}.", i = "Ensure knots are sorted and unique."))
      }
      if (mode %in% c(4, 6) && any(knots < 0 | knots > 1)) {
        cli::cli_abort(c(
          x = "Expit-spline knots must lie in [0, 1] for {short_suffix}.",
          i = "Provide explicit expit-scale knots, or provide expit-scale x values so knots can be derived from them."
        ))
      }
      
      if (mode %in% c(2, 4, 5, 6)) {
        # I-spline: coeff length should match knot/degree conventions.
        degree <- standata[[paste0("spline_degree_", short_suffix)]]
        expected_len_a <- length(knots) + degree
        expected_len_b <- length(knots) + degree - 1
        if (!(coeff_len %in% c(expected_len_a, expected_len_b))) {
          cli::cli_abort(c(
            x = "I-spline coeff length mismatch for {short_suffix}.",
            i = "Expected {expected_len_a} or {expected_len_b} coefficients, got {coeff_len}."
          ))
        }
      }
      if (mode %in% c(3, 7) && coeff_len != length(knots)) {
        cli::cli_abort(c(
          x = "Piecewise-linear coefficient count mismatch for {short_suffix}.",
          i = "Expected one relative association ordinate per knot."
        ))
      }
    }
  }
  
  return(TRUE)
}

#' Example: Build Transform Specification from Formula String
#'
#' Dummy function to print help this help page.
#' See `build_standata_transforms()` for details on how to construct a transform
#' for the association term.
#'
#' @export
#'
#' @examples
#' # Example 1: Functional transformation (bytecode-backed)
#' spec_cv <- list(
#'   type = "functional",
#'   expr = ~ (log(sqrt(x + 1/inv_logit(3*x - 3))))^2
#' )
#'
#' # Example 2: Piecewise-linear
#' spec_cs <- list(
#'   type = "pwlin",
#'   knots = c(-2, -1, 0, 1, 2),
#'   direction = "increasing"
#' )
#'
#' # Example 3: I-spline (monotonic)
#' # - knots define basis locations (first/last are boundary knots)
#' # - coeff defines the monotone shape directly
#' spec_corr <- list(
#'   type = "ispline",
#'   knots = c(-1, 0, 1),
#'   coeff = c(0, 0.5, 1, 1.2, 1.5),
#'   degree = 3
#' )
#'
#' # Example 4: Penalised monotone I-spline (fit in R)
#' # - x is the input scale of the raw association feature
#' # - y is the desired transformed output at each x
#' # - lambda controls smoothness (higher = smoother)
#' spec_corr_pen <- list(
#'   type = "ispline_penalised",
#'   x = seq(-2, 2, length.out = 50),
#'   y = exp(seq(-2, 2, length.out = 50)),
#'   n_knots = 6,
#'   degree = 3,
#'   lambda = 1.0,
#'   direction = "increasing"
#' )
#'
#' # Example 5: Penalised monotone I-spline (Stan-estimated coefficients)
#' # - omit y so Stan learns the monotone shape directly
#' # - if knots are omitted, n_knots = 6 and x quantiles define them
#' spec_corr_pen_stan <- list(
#'   type = "ispline_penalised",
#'   x = seq(-2, 2, length.out = 50),
#'   n_knots = 6,
#'   degree = 3,
#'   lambda = 1.0
#' )
#'
#' # Example 6: Penalised monotone I-spline on expit(x)
#' spec_corr_expit <- list(
#'   type = "ispline_expit_penalised",
#'   x = seq(0.02, 0.98, length.out = 50),
#'   y = seq(0.02, 0.98, length.out = 50)^0.75,
#'   n_knots = 6,
#'   degree = 3,
#'   lambda = 1.0,
#'   direction = "increasing"
#' )
#'
#' # Combine and build standata
#' transforms <- joinme_tf(
#'   cv_total = spec_cv,
#'   cs_total = spec_cs,
#'   corr = spec_corr
#' )
#' standata_tf <- joinme:::build_standata_transforms(transforms)
example_transform_spec <- function() {
  help('example_transform_spec')
}

#' Penalised Monotone I-spline Transform Spec
#'
#' @description
#' Fits a monotone I-spline transformation with a smoothness penalty and
#' returns a transform specification compatible with `build_standata_transforms()`.
#'
#' This helper is the legacy plug-in constructor for `type = "ispline_penalised"`.
#' Its defaults are `n_knots = 6`, `degree = 3`, `lambda = 1.0`,
#' `weights = NULL` (equal weights), and `diff_order = 2`.
#' Returned coefficients follow the same anchored convention as the Stan-
#' estimated path: increasing splines run from `0` to `1`, decreasing splines
#' run from `1` to `0`, and interior coefficients stay monotone in the chosen
#' direction.
#' Example:
#' `penalised_ispline_transform(x = seq(-2, 2, length.out = 50), y = exp(seq(-2, 2, length.out = 50)))`
#'
#' @param x Numeric vector of input values on the raw feature scale to transform.
#' @param y Numeric vector of target transformed values at `x` (same length as `x`).
#' @param knots Optional numeric vector of knots (including boundary knots).
#' @param n_knots Integer. If knots are not provided, number of knots to use
#'   (including boundary knots). Default 6.
#' @param degree Integer spline degree (default 3).
#' @param lambda Non-negative smoothness penalty weight (default 1.0).
#'   Larger values produce smoother fitted transforms; smaller values allow
#'   more local curvature.
#' @param direction Monotone direction. Accepted values are `"increasing"`
#'   and `"decreasing"`. When omitted, the direction is inferred from the
#'   supplied `(x, y)` pairs.
#' @param weights Optional non-negative weights (same length as x).
#' @param diff_order Integer difference order for penalty (default 2).
#'
#' @return A list suitable for `build_standata_transforms()`.
#'
#' @export
penalised_ispline_transform <- function(
  x,
  y,
  knots = NULL,
  n_knots = 6,
  degree = 3,
  lambda = 1.0,
  direction = NULL,
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
    direction = direction,
    weights = weights,
    diff_order = diff_order
  )
  .make_penalised_ispline_transform(spec)
}

#' @rdname penalised_ispline_transform
#' @param ... Arguments passed to [penalised_ispline_transform()].
#' @export
penalized_ispline_transform <- function(...) {
  penalised_ispline_transform(...)
}

#' Build Stan-estimated penalised I-spline specification
#'
#' @description
#' Construct the standata payload for a penalised I-spline whose coefficients
#' will be estimated inside Stan. This path keeps the knot sequence fixed while
#' initializing a monotone coefficient vector and the associated penalty
#' metadata.
#'
#' @param spec Named list describing the transform. Accepted fields include
#'   `type`, `degree`, `lambda`, `direction`, `knots`, `n_knots`, `x`, and
#'   `n_coeff`.
#'   For expit-domain variants, any supplied `x` or `knots` values must already
#'   lie in `[0, 1]`.
#'
#' @return A canonicalised transform spec list suitable for
#'   `build_standata_transforms()`.
#' @keywords internal
#' @noRd
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

  spline_direction <- .get_monotone_direction(spec$direction, default = "increasing")

  list(
    type = if (.transform_uses_expit_input(spec)) "ispline_expit" else "ispline",
    knots = knots,
    raw_knots = if (is.null(knots_raw)) NULL else as.numeric(knots_raw),
    coeff = rep(0, n_coeff),
    degree = degree,
    direction = .monotone_direction_label(spline_direction),
    estimate_spline = 1L,
    lambda_spline = as.numeric(lambda),
    n_free_spline = as.integer(max(0L, n_coeff - 1L)),
    spline_direction = spline_direction
  )
}

#' Fit penalised I-spline coefficients in R
#'
#' @description
#' Fit a monotone penalised I-spline to training pairs `(x, y)` in R and return
#' the fixed-knot, fixed-coefficient specification that Stan later evaluates at
#' runtime. This is the plug-in path used when `y` is supplied.
#'
#' @param spec Named list describing the transform. Accepted fields include
#'   `type`, `x`, `y`, `knots`, `n_knots`, `degree`, `lambda`, `direction`,
#'   `weights`, and `diff_order`. For expit-domain variants, any supplied `x`
#'   or `knots` values
#'   must already lie in `[0, 1]`.
#'
#' @return A canonicalised fixed I-spline transform spec list.
#' @keywords internal
#' @noRd
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
  spline_direction <- .get_monotone_direction(spec$direction, default = .infer_monotone_direction(x, y))
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
  target_y <- if (spline_direction < 0L) max(y, na.rm = TRUE) - y else y - min(y, na.rm = TRUE)
  target_y <- as.numeric(target_y)
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
    b <- .anchored_coeff_from_z(z)
    fit_vals <- as.numeric(B_w %*% b)
    if (spline_direction < 0L) {
      fit_vals <- w_sqrt * (sum(b) - as.numeric(basis %*% b))
    }
    r <- y_w - fit_vals
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
    degree = degree,
    direction = .monotone_direction_label(spline_direction),
    spline_direction = spline_direction
  )
}

#' Default Backward-Compatible Transforms
#'
#' If no transforms specified, create identity specs so existing models
#' continue to work unchanged.
#'
#' @keywords internal
#' @noRd
default_identity_transforms <- function() {
  list(
    cv_tot = list(type = "identity"),
    cs_tot = list(type = "identity"),
    corr = list(type = "identity"),
    vcov = list(type = "identity")
  )
}
