#' Declare association transformations with validation
#'
#' @description
#' `joinme_tf()` is a user-facing wrapper for transformation declarations used by
#' `joinme()` and `simulate_joinme()`. It validates channel names, normalises
#' common shorthands, and returns a structured object that can be passed directly
#' as the `transforms` argument.
#'
#' Supported channels are `cv_total`, `cs_total`, `cv_mean`, `cs_mean`,
#' `cv_marker`, `cs_marker`, `corr`, and `vcov`. The `corr` channel refers to
#' off-diagonal entries of the subject-specific Cholesky-correlation factor `K`,
#' while `vcov` refers to those same off-diagonal `K` entries together with the
#' subject-specific standard deviations.
#'
#' Supported shorthands per channel:
#' - `"identity"`, `"ispline"`, `"ispline_penalised"`, `"ispline_penalized"`,
#'   `"ispline_expit"`, `"ispline_expit_penalised"`, `"ispline_expit_penalized"`,
#'   `"ispline_exp_penalised"`, `"ispline_exp_penalized"`, `"pwlin"`,
#'   `"functional"`
#' - a formula such as `~ log1p(x)`, which is interpreted as
#'   `list(type = "functional", expr = ~ log1p(x))`
#' - a full named list specification such as
#'   `list(type = "ispline", knots = ..., coeff = ..., degree = 3)`
#'
#' @param ... Named channel specifications.
#' @param .validate Logical; if `TRUE` (default), validate the resulting
#'   transform declarations immediately.
#'
#' @return An object of class `joinme_tf`.
#' @aliases jm_tf
#' @export
#'
#' @examples
#' tf <- joinme_tf(
#'   cv_total = "identity",
#'   corr = ~ -x,
#'   vcov = "identity",
#'   cv_marker = list(
#'     type = "pwlin",
#'     knots = c(-1, 0, 1),
#'     direction = "increasing"
#'   )
#' )
#' print(tf)
joinme_tf <- function(..., .validate = TRUE) {
  args <- list(...)
  if (length(args) == 1L && is.null(names(args)) &&
      (inherits(args[[1]], "joinme_tf") || is.list(args[[1]]))) {
    args <- args[[1]]
  }
  .normalise_joinme_tf_input(args, validate = .validate)
}

#' @rdname joinme_tf
#' @export
jm_tf <- joinme_tf

#' @export
print.joinme_tf <- function(x, ...) {
  specs <- unclass(x)
  if (length(specs) == 0) {
    cat("<joinme_tf: no explicit transforms; identity defaults apply>\n")
    return(invisible(x))
  }

  describe_spec <- function(spec) {
    type <- spec$type %||% "identity"
    if (identical(type, "functional")) {
      expr_txt <- tryCatch(paste(deparse(spec$expr), collapse = ""), error = function(e) "<expr>")
      return(paste0(type, " (", expr_txt, ")"))
    }
    if (type %in% c("ispline", "ispline_penalised", "ispline_expit", "ispline_expit_penalised")) {
      knot_txt <- if (!is.null(spec$knots)) paste0("knots=", length(spec$knots)) else paste0("n_knots=", spec$n_knots %||% "?")
      return(paste(type, knot_txt, paste0("degree=", spec$degree %||% 3L), sep = ", "))
    }
    if (identical(type, "pwlin")) {
      direction <- .monotone_direction_label(spec$direction %||% 1L)
      return(paste0(type, ", knots=", length(spec$knots %||% spec$cutpoints %||% spec$x %||% numeric(0)), ", ", direction))
    }
    type
  }

  tbl <- data.frame(
    channel = names(specs),
    specification = vapply(specs, describe_spec, character(1)),
    stringsAsFactors = FALSE
  )
  cat("Transformation specification for Joint Mixed Effects model\n")
  print(tbl, row.names = FALSE)
  invisible(x)
}

#' Declare a Normal coefficient prior
#'
#' @description
#' Creates a JoiNMe prior declaration without evaluating a density in R. The
#' location and scale may be scalars or vectors; scalar values are recycled to
#' the number of coefficients in the selected model block.
#'
#' The `prior_` prefix is intentional. It avoids masking the short constructor
#' names exported by other Bayesian modelling packages.
#'
#' @param mu Prior location or vector of coefficient-specific locations.
#' @param scale Positive prior scale or vector of coefficient-specific scales.
#'   Numeric lists, such as `list(1, 1, 2)`, are accepted and flattened in
#'   their supplied order. The same convention applies to `mu`.
#'
#' @return A `joinme_prior_spec` object for use inside [jm_prior()].
#' @export
prior_normal <- function(mu = 0, scale = 1) {
  make_prior_dist("normal", mu = mu, scale = scale)
}

#' Declare a Student-t coefficient prior
#'
#' @inheritParams prior_normal
#' @param df Positive fixed degrees of freedom. A single value is shared by the
#'   complete parameter block.
#'
#' @return A `joinme_prior_spec` object for use inside [jm_prior()].
#' @export
prior_student_t <- function(df = 6, mu = 0, scale = 1) {
  make_prior_dist(
    "student_t",
    mu = mu,
    scale = scale,
    df = df
  )
}

#' Declare a Laplace coefficient prior
#'
#' @description
#' The Laplace distribution is Stan's double-exponential distribution. Its
#' `scale` is the exponential-decay scale, not its standard deviation; the
#' standard deviation is `sqrt(2) * scale`.
#'
#' @inheritParams prior_normal
#'
#' @return A `joinme_prior_spec` object for use inside [jm_prior()].
#' @export
prior_laplace <- function(mu = 0, scale = 1) {
  make_prior_dist("laplace", mu = mu, scale = scale)
}

#' Declare a regularised horseshoe prior
#'
#' @description
#' Declares the regularised horseshoe hierarchy used by brms and rstanarm. The
#' coefficient is conditionally Normal with local and global half-Student-t
#' scales and a finite Student-t slab. The slab regularises very large signals,
#' whilst the local scales permit coefficient-specific escape from shrinkage.
#'
#' @param df Positive fixed degrees of freedom for local shrinkage scales.
#' @param global_df Positive fixed degrees of freedom for the global scale.
#' @param global_scale Positive global scale. Smaller values express stronger
#'   prior sparsity.
#' @param slab_df Positive fixed degrees of freedom for the regularising slab.
#' @param slab_scale Positive scale of the regularising slab.
#'
#' @return A `joinme_prior_spec` object for use inside [jm_prior()].
#' @export
prior_horseshoe <- function(
  df = 1,
  global_df = 1,
  global_scale = 1,
  slab_df = 4,
  slab_scale = 2
) {
  make_prior_dist(
    "horseshoe",
    mu = 0,
    scale = 1,
    df = df,
    global_df = global_df,
    global_scale = global_scale,
    slab_df = slab_df,
    slab_scale = slab_scale
  )
}

#' Declare an LKJ correlation prior
#'
#' @param eta Positive LKJ concentration. `eta = 1` is uniform over correlation
#'   matrices; values above one favour correlations nearer zero.
#'
#' @return A `joinme_lkj_prior` object for use inside [jm_prior()].
#' @export
prior_lkj <- function(eta = 1) {
  if (!is.numeric(eta) || length(eta) != 1L || !is.finite(eta) || eta <= 0) {
    cli::cli_abort("{.arg eta} must be one positive finite number.")
  }
  structure(
    list(family = "lkj", eta = as.numeric(eta)),
    class = c("joinme_lkj_prior", "list")
  )
}

#' Build one validated JoiNMe coefficient-prior declaration
#'
#' @keywords internal
#' @noRd
make_prior_dist <- function(
  family,
  mu,
  scale,
  df = NA_real_,
  global_df = NA_real_,
  global_scale = NA_real_,
  slab_df = NA_real_,
  slab_scale = NA_real_
) {
  normalise_numeric_values <- function(value, argument, positive = FALSE) {
    # Lists make long coefficient-wise declarations easier to read. Requiring
    # one numeric value per element preserves an unambiguous parameter order.
    if (is.list(value)) {
      valid_elements <- vapply(
        value,
        function(element) is.numeric(element) && length(element) == 1L,
        logical(1)
      )
      if (!all(valid_elements)) {
        cli::cli_abort("{.arg {argument}} must be a numeric vector or a list of single numeric values.")
      }
      value <- unlist(value, recursive = FALSE, use.names = TRUE)
    }
    if (!is.numeric(value) || length(value) < 1L || any(!is.finite(value))) {
      cli::cli_abort("{.arg {argument}} must contain finite numeric values.")
    }
    if (isTRUE(positive) && any(value <= 0)) {
      cli::cli_abort("{.arg {argument}} must contain positive finite values.")
    }
    as.numeric(value)
  }
  mu <- normalise_numeric_values(mu, "mu") # coefficient locations in declared order
  scale <- normalise_numeric_values(scale, "scale", positive = TRUE) # positive coefficient scales in declared order
  scalar_hyperparameters <- list(
    df = df,
    global_df = global_df,
    global_scale = global_scale,
    slab_df = slab_df,
    slab_scale = slab_scale
  ) # fixed scalar degrees of freedom and horseshoe hyper-scales
  for (hyperparameter_name in names(scalar_hyperparameters)) {
    hyperparameter_value <- scalar_hyperparameters[[hyperparameter_name]] # one fixed hyperparameter supplied to Stan data
    if (
      !is.numeric(hyperparameter_value) || length(hyperparameter_value) != 1L ||
        (!is.na(hyperparameter_value) &&
          (!is.finite(hyperparameter_value) || hyperparameter_value <= 0))
    ) {
      cli::cli_abort(
        "Prior hyperparameter {.arg {hyperparameter_name}} must be one positive finite number."
      )
    }
  }
  structure(
    list(
      family = family,
      mu = mu,
      scale = scale,
      df = as.numeric(df),
      global_df = as.numeric(global_df),
      global_scale = as.numeric(global_scale),
      slab_df = as.numeric(slab_df),
      slab_scale = as.numeric(slab_scale)
    ),
    class = c("joinme_prior_spec", "list")
  )
}

#' Declare priors with validation
#'
#' @description
#' `joinme_priors()` is a user-facing wrapper for prior declarations passed to
#' `joinme()`. It validates the exposed prior components and returns a
#' structured object that can be supplied directly as the `priors` argument.
#'
#' Exposed components are:
#' - `beta`: longitudinal fixed-effect coefficients;
#' - `alpha`: coefficients multiplying transformed association summaries;
#' - `iota`: fitted affine shifts for functional association transforms;
#' - `marker`: standardised marker-level random effects;
#' - `marker_weight`: estimated marker-weight perturbations;
#' - `vcov_sd`: intercepts and slopes of the standard-deviation regression;
#' - `vcov_corr`: intercepts and slopes of the off-diagonal correlation regression;
#' - `lkj`: random-effect correlation matrices;
#' - `class_probability`: Dirichlet concentration for baseline class probabilities
#' - `class_regression`: class-membership regression coefficients.
#'
#' @param beta,alpha,iota,vcov_sd,vcov_corr Prior declarations made with [prior_normal()],
#'   [prior_student_t()], [prior_laplace()] or [prior_horseshoe()]. Numeric
#'   values remain supported as deprecated shorthand for Student-t(6) scales.
#'   `vcov_sd` is packed as all standard-deviation intercepts followed by their
#'   row-major slope matrix. `vcov_corr` uses the same order for off-diagonal
#'   partial-correlation coordinates `(2,1), (3,1), (3,2), ...`.
#' @param marker Prior family for standardised marker-level random effects.
#'   Its location and scale are fixed at zero and one. Consequently a supplied
#'   Normal, Student-t or Laplace declaration must retain `mu = 0` and
#'   `scale = 1`; horseshoe hyper-scales must retain their unit defaults.
#' @param marker_weight Prior family for estimated marker-weight perturbations.
#'   The location and scale are likewise fixed at zero and one.
#' @param lkj An object from [prior_lkj()] or a positive numeric concentration
#'   retained for compatibility.
#' @param class_probability Positive Dirichlet concentration for the baseline
#'   class probabilities. A scalar is repeated over classes; a vector may give
#'   one concentration per class.
#' @param class_regression Prior declaration for coefficients from
#'   `formulaClass`, made with the same constructors accepted by `beta`,
#'   `alpha`, and `iota`. A positive numeric value remains accepted as shorthand
#'   for `prior_normal(scale = value)`.
#' @param .validate Logical; if `TRUE` (default), validate the resulting prior
#'   declarations immediately.
#'
#' @return An object of class `joinme_priors`.
#' @aliases joinme_prior jm_priors jm_prior
#' @export
#'
#' @examples
#' pri <- jm_prior(
#'   beta = prior_normal(mu = c(0, 1), scale = c(1, 0.5)),
#'   alpha = prior_student_t(df = 4, scale = 1),
#'   iota = prior_laplace(scale = 0.75),
#'   marker = prior_normal(),
#'   marker_weight = prior_student_t(df = 6),
#'   vcov_sd = prior_normal(scale = 1),
#'   vcov_corr = prior_student_t(df = 4, scale = 0.75),
#'   lkj = prior_lkj(2),
#'   class_probability = c(2, 2, 2),
#'   class_regression = prior_normal(scale = 1)
#' )
#' print(pri)
joinme_priors <- function(
  beta = NULL, # prior scale for fixed-effect coefficients
  alpha = NULL, # prior scale for association coefficients
  iota = NULL, # prior scale for fitted transformation shifts
  marker = NULL, # unit-scale family for standardised marker random effects
  marker_weight = NULL, # unit-scale family for marker-weight perturbations
  vcov_sd = NULL, # prior for standard-deviation covariance-regression coefficients
  vcov_corr = NULL, # prior for off-diagonal correlation-regression coefficients
  lkj = NULL, # LKJ concentration for correlation matrices
  class_probability = 1, # Dirichlet concentration for baseline class weights
  class_regression = NULL, # family and hyperparameters for class-regression coefficients
  .validate = TRUE # whether to validate the assembled prior specification
) {
  .joinme_priors_(
    list(
      beta = beta,
      alpha = alpha,
      iota = iota,
      marker = marker,
      marker_weight = marker_weight,
      vcov_sd = vcov_sd,
      vcov_corr = vcov_corr,
      lkj = lkj,
      class_probability = class_probability,
      class_regression = class_regression
    ),
    validate = .validate
  )
}

#' @rdname joinme_priors
#' @export
jm_priors <- joinme_priors

#' @rdname joinme_priors
#' @export
jm_prior <- joinme_priors

#' @rdname joinme_priors
#' @export
joinme_prior <- joinme_priors

#' @export
print.joinme_priors <- function(x, ...) {
  priors <- unclass(x)
  describe_prior <- function(prior) {
    if (is.null(prior)) return("default")
    if (inherits(prior, "joinme_lkj_prior")) {
      return(paste0("lkj(eta = ", prior$eta, ")"))
    }
    if (inherits(prior, "joinme_prior_spec")) {
      if (identical(prior$family, "horseshoe")) {
        return(paste0(
          "regularised horseshoe(global_scale = ", prior$global_scale,
          ", slab_scale = ", prior$slab_scale, ")"
        ))
      }
      return(paste0(
        prior$family,
        "(mu = ", paste(prior$mu, collapse = ", "),
        "; scale = ", paste(prior$scale, collapse = ", "),
        if (identical(prior$family, "student_t")) {
          paste0("; df = ", prior$df)
        } else "",
        ")"
      ))
    }
    paste(prior, collapse = ", ")
  }
  class_probability_txt <- paste(priors$class_probability, collapse = ", ")

  tbl <- data.frame(
    component = c(
      "beta",
      "alpha",
      "iota",
      "marker",
      "marker_weight",
      "vcov_sd",
      "vcov_corr",
      "lkj",
      "class_probability",
      "class_regression"
    ),
    value = c(
      describe_prior(priors$beta),
      describe_prior(priors$alpha),
      describe_prior(priors$iota),
      describe_prior(priors$marker),
      describe_prior(priors$marker_weight),
      describe_prior(priors$vcov_sd),
      describe_prior(priors$vcov_corr),
      describe_prior(priors$lkj),
      class_probability_txt,
      describe_prior(priors$class_regression)
    ),
    stringsAsFactors = FALSE
  )
  cat("Prior specification for Joint Mixed Effects model\n")
  print(tbl, row.names = FALSE)
  invisible(x)
}

#' Build specs for baseline hazard
#'
#' @description
#' This function provides a formal mean to specify the baseline hazard specification for the survival submodel.
#' It evaluates nothing, but it checks the validity of the input and returns a structured object that can be passed to `joinme()`.
#' @param type Baseline hazard type: "bs", "ns", or "formula".
#' @param n_knots Number of internal knots for spline baseline hazards.
#' @param knots Optional numeric vector of internal knots for spline baseline hazards.
#' @param degree Degree of spline basis for baseline hazard.
#' @param formula Formula for baseline hazard when `type = "formula"`. Ignored otherwise.
#' @return An object of class `joinme_basehaz`.
#' @aliases jm_basehaz
#' @export
#' @examples
#' basehaz_spec <- joinme_basehaz(
#'  type = "bs",
#'  n_knots = 5,
#'  degree = 3
#' )
joinme_basehaz <- function(
  type = c('bs', 'ns', 'formula'),
  n_knots = 5L,
  knots = NULL,
  degree = 3L,
  formula = ~ 1 + time
){
  type <- match.arg(type)
  if (type %in% c('bs', 'ns')) {
    if (!is.null(knots)) {
      if (!is.numeric(knots) || length(knots) < 1) {
        cli::cli_abort(c(
          x = "When specifying baseline hazard knots, {.arg knots} must be a numeric vector of length >= 1."
        ))
      }
    }
    if (!is.numeric(n_knots) || length(n_knots) != 1 || n_knots < 1) {
      cli::cli_abort(c(
        x = "When specifying baseline hazard spline, {.arg n_knots} must be a single positive integer."
      ))
    }
    if (!is.numeric(degree) || length(degree) != 1 || degree < 1) {
      cli::cli_abort(c(
        x = "When specifying baseline hazard spline, {.arg degree} must be a single positive integer."
      ))
    }
  } else if (type == 'formula') {
    formula <- as.formula(formula)
  }

  structure(
    list(
      type = type,
      n_knots = n_knots,
      knots = knots,
      degree = degree,
      formula = formula
    ),
    class = "joinme_basehaz"
  )
}

#' @rdname joinme_basehaz
#' @export
jm_basehaz <- joinme_basehaz

#' Build named conditional-effects profiles
#'
#' @description
#' `make_conditions()` is a thin wrapper around [brms::make_conditions()] so the
#' resulting labelled condition tables can be passed directly to
#' [conditional_effects.JoiNMeFit()]. Conditioning belongs to the conditional-
#' effects estimand; neither [predict.JoiNMeFit()] nor `plot.JoiNMeFit()` accepts
#' a `condition` argument.
#'
#' @param x A data frame containing the baseline covariates used to define the
#'   conditioning rows.
#' @param ... Additional arguments passed to [brms::make_conditions()].
#'
#' @return A data frame with one row per requested condition and a `cond__`
#'   column containing the profile labels used in conditional-effects facets.
#'
#' @seealso [conditional_effects.JoiNMeFit()],
#'   [brms::make_conditions()]
#' @export
#' @importFrom brms make_conditions
make_conditions <- function(x, ...) {
  brms::make_conditions(x = x, ...)
}

#' Functions admitting fitted affine shifts in functional transforms
#'
#' @description
#' Lists nonlinear operations for which `joinme_tf()` may estimate an
#' operation-specific intercept or slope before applying the transformation.
#' Both directions of the standard normal transformation are included:
#' `Phi`/`pnorm` for the CDF and `inv_Phi`/`qnorm`/`probit` for the quantile.
#'
#' @return Character vector of accepted, lower-case parser names.
#' @keywords internal
#' @noRd
.fit_affine_shift_supported_functions <- function() {
  c(
    "log", "exp", "sqrt", "inv_logit", "sigmoid", "expit", "softmax",
    "logit", "rec", "sin", "cos", "tan", "abs", "sinh", "cosh", "tanh",
    "asinh", "acosh", "atanh", "softplus", "log1p_exp", "cbrt", "phi",
    "pnorm", "inv_phi", "qnorm", "probit", "power"
  )
}

#' @keywords internal
.parse_fit_affine_shift_flag <- function(value, arg_name, term_name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    cli::cli_abort(c(
      x = "Functional transform option {.arg {arg_name}} for {.val {term_name}} must be TRUE or FALSE.",
      i = "Use syntax such as {.code joinme_tf(cv_total = ~ expit(x, intercept = TRUE, slope = TRUE))}."
    ))
  }
  isTRUE(value)
}

#' @keywords internal
.extract_fit_affine_shift_node <- function(
  node,
  term_name = "transform",
  path = "root",
  intercept_count = 0L,
  slope_count = 0L
) {
  if (!is.call(node)) {
    return(list(
      node = node,
      iota_nodes = list(),
      n_iota_intercept = as.integer(intercept_count),
      n_iota_slope = as.integer(slope_count)
    ))
  }

  op <- as.character(node[[1]])
  op_lower <- if (grepl("^[[:alpha:].][[:alnum:]_.]*$", op)) tolower(op) else op
  args <- as.list(node[-1])
  arg_names <- names(args) %||% rep("", length(args))
  intercept_idx <- which(arg_names == "intercept")
  slope_idx <- which(arg_names == "slope")

  if (length(intercept_idx) > 1L || length(slope_idx) > 1L) {
    cli::cli_abort(c(
      x = "Functional transform {.val {term_name}} may declare {.arg intercept} and {.arg slope} at most once.",
      i = "Use each fit-only affine-shift option no more than once."
    ))
  }

  has_affine_shift <- length(intercept_idx) > 0L || length(slope_idx) > 0L
  shift_intercept <- FALSE
  shift_slope <- FALSE
  if (has_affine_shift) {
    if (!(op_lower %in% .fit_affine_shift_supported_functions())) {
      cli::cli_abort(c(
        x = "Functional transform {.val {term_name}} uses fit-only affine-shift options with unsupported function {.val {op}}.",
        i = "Use fit-only {.arg intercept}/{.arg slope} with supported nonlinear functional transforms such as expit, softplus, sinh, cosh, tanh, or power."
      ))
    }
    if (length(intercept_idx) > 0L) {
      shift_intercept <- .parse_fit_affine_shift_flag(args[[intercept_idx]], "intercept", term_name)
    }
    if (length(slope_idx) > 0L) {
      shift_slope <- .parse_fit_affine_shift_flag(args[[slope_idx]], "slope", term_name)
    }
  }

  keep_idx <- setdiff(seq_along(args), c(intercept_idx, slope_idx))
  args_clean <- args[keep_idx]
  arg_names_clean <- arg_names[keep_idx]
  iota_nodes <- list()

  if (length(args_clean) > 0L) {
    cleaned_args <- vector("list", length(args_clean))
    for (idx in seq_along(args_clean)) {
      child <- .extract_fit_affine_shift_node(
        args_clean[[idx]],
        term_name = term_name,
        path = paste0(path, "/", idx),
        intercept_count = intercept_count,
        slope_count = slope_count
      )
      cleaned_args[[idx]] <- child$node
      intercept_count <- as.integer(child$n_iota_intercept)
      slope_count <- as.integer(child$n_iota_slope)
      iota_nodes <- c(iota_nodes, child$iota_nodes)
    }
  } else {
    cleaned_args <- list()
  }

  intercept_index <- 0L
  slope_index <- 0L
  if (isTRUE(shift_intercept)) {
    intercept_count <- intercept_count + 1L
    intercept_index <- intercept_count
  }
  if (isTRUE(shift_slope)) {
    slope_count <- slope_count + 1L
    slope_index <- slope_count
  }
  if (isTRUE(shift_intercept) || isTRUE(shift_slope)) {
    iota_nodes <- c(iota_nodes, list(list(
      path = path,
      function_name = op_lower,
      intercept_index = intercept_index,
      slope_index = slope_index,
      estimate_iota_intercept = shift_intercept,
      estimate_iota_slope = shift_slope
    )))
  }

  parts <- c(list(as.name(op_lower)), cleaned_args)
  names(parts) <- c("", arg_names_clean)

  list(
    node = as.call(parts),
    iota_nodes = iota_nodes,
    n_iota_intercept = as.integer(intercept_count),
    n_iota_slope = as.integer(slope_count)
  )
}

#' @keywords internal
.extract_fit_affine_shift_spec <- function(expr, term_name = "transform") {
  expr_call <- .coerce_transform_expr(expr)
  extracted <- .extract_fit_affine_shift_node(expr_call, term_name = term_name)
  expr_call <- extracted$node

  list(
    expr = stats::as.formula(paste0("~ ", paste(deparse(expr_call, width.cutoff = 500L), collapse = ""))),
    estimate_iota_intercept = isTRUE(extracted$n_iota_intercept > 0L),
    estimate_iota_slope = isTRUE(extracted$n_iota_slope > 0L),
    n_iota_intercept = as.integer(extracted$n_iota_intercept),
    n_iota_slope = as.integer(extracted$n_iota_slope),
    iota_nodes = extracted$iota_nodes
  )
}

#' @keywords internal
.transform_iota_nodes <- function(spec) {
  nodes <- spec$iota_nodes %||% list()
  if (length(nodes) > 0L) {
    return(nodes)
  }
  if (is.null(spec) || !identical(spec$type %||% NULL, "functional")) {
    return(list())
  }

  n_iota_intercept <- as.integer(spec$n_iota_intercept %||% if (isTRUE(spec$estimate_iota_intercept)) 1L else 0L)
  n_iota_slope <- as.integer(spec$n_iota_slope %||% if (isTRUE(spec$estimate_iota_slope)) 1L else 0L)
  if (n_iota_intercept < 1L && n_iota_slope < 1L) {
    return(list())
  }

  expr_call <- tryCatch(.coerce_transform_expr(spec$expr), error = function(e) NULL)
  if (!is.call(expr_call)) {
    return(list())
  }

  list(list(
    path = "root",
    function_name = tolower(as.character(expr_call[[1]])),
    intercept_index = if (n_iota_intercept > 0L) 1L else 0L,
    slope_index = if (n_iota_slope > 0L) 1L else 0L,
    estimate_iota_intercept = n_iota_intercept > 0L,
    estimate_iota_slope = n_iota_slope > 0L
  ))
}

#' @keywords internal
.normalise_joinme_tf_input <- function(transforms = NULL, validate = TRUE) {
  if (is.null(transforms)) {
    return(structure(list(), class = c("joinme_tf", "list")))
  }

  specs <- if (inherits(transforms, "joinme_tf")) unclass(transforms) else transforms
  if (length(specs) == 1L && is.null(names(specs)) && is.list(specs[[1]]) &&
      !is.null(names(specs[[1]])) && !"type" %in% names(specs[[1]])) {
    specs <- specs[[1]]
  }

  if (!is.list(specs)) {
    cli::cli_abort(c(
      x = "{.arg transforms} must be a named list or a {.fn joinme_tf} object.",
      i = "Example: joinme_tf(cv_total = 'identity', corr = ~ -x, vcov = 'identity')."
    ))
  }
  if (length(specs) > 0 && (is.null(names(specs)) || any(names(specs) %in% c("", NA_character_)))) {
    cli::cli_abort(c(
      x = "{.arg transforms} must be named by association channel.",
      i = "Use names such as cv_total, cs_total, cv_mean, cs_mean, cv_marker, cs_marker, corr, vcov."
    ))
  }

  allowed <- c("cv_total", "cs_total", "cv_mean", "cs_mean", "cv_marker", "cs_marker", "corr", "vcov")
  bad <- setdiff(names(specs), allowed)
  if (length(bad) > 0) {
    cli::cli_abort(c(
      x = "Unknown transform channel(s): {paste(bad, collapse = ', ')}.",
      i = "Allowed channels: {paste(allowed, collapse = ', ')}."
    ))
  }

  specs <- stats::setNames(lapply(names(specs), function(nm) {
    .normalise_joinme_tf_spec(specs[[nm]], term_name = nm)
  }), names(specs))

  if (isTRUE(validate) && length(specs) > 0) {
    build_standata_transforms(specs)
  }

  structure(specs, class = c("joinme_tf", "list"))
}

#' @keywords internal
.normalise_joinme_tf_spec <- function(spec, term_name = "transform") {
  allowed_types <- c(
    "identity",
    "functional",
    "ispline",
    "ispline_penalised",
    "ispline_penalized",
    "ispline_expit",
    "ispline_expit_penalised",
    "ispline_expit_penalized",
    "ispline_exp_penalised",
    "ispline_exp_penalized",
    "pwlin"
  )

  if (is.null(spec)) {
    return(list(
      type = "identity",
      estimate_iota_intercept = FALSE,
      estimate_iota_slope = FALSE,
      n_iota_intercept = 0L,
      n_iota_slope = 0L,
      iota_nodes = list()
    ))
  }
  if (inherits(spec, "formula")) {
    parsed <- .extract_fit_affine_shift_spec(spec, term_name = term_name)
    return(c(list(type = "functional"), parsed))
  }
  if (is.character(spec) && length(spec) == 1L) {
    if (!(spec %in% allowed_types)) {
      cli::cli_abort(c(
        x = "Unknown transform shorthand {.val {spec}} for {.val {term_name}}.",
        i = "Use one of: {paste(allowed_types, collapse = ', ')} or a formula such as ~ log1p(x)."
      ))
    }
    return(list(
      type = .canonicalise_transform_type(spec),
      estimate_iota_intercept = FALSE,
      estimate_iota_slope = FALSE,
      n_iota_intercept = 0L,
      n_iota_slope = 0L,
      iota_nodes = list()
    ))
  }
  if (!is.list(spec)) {
    cli::cli_abort(c(
      x = "Transform {.val {term_name}} must be declared as a list, a formula, or a supported shorthand string.",
      i = "Example: joinme_tf({term_name} = list(type = 'functional', expr = ~ log1p(x)))."
    ))
  }

  spec <- as.list(spec)
  if (is.null(spec$type)) {
    if (!is.null(spec$expr)) {
      spec$type <- "functional"
    } else if (!is.null(spec$coeff) && !is.null(spec$knots)) {
      spec$type <- "ispline"
    } else {
      cli::cli_abort(c(
        x = "Transform {.val {term_name}} is missing {.arg type}.",
        i = "Use a full specification like list(type = 'ispline', knots = ..., coeff = ...) or a shorthand such as ~ log1p(x)."
      ))
    }
  }
  spec$type <- .canonicalise_transform_type(spec$type)
  if (identical(spec$type, "functional")) {
    parsed <- .extract_fit_affine_shift_spec(spec$expr, term_name = term_name)
    spec$expr <- parsed$expr
    spec$estimate_iota_intercept <- parsed$estimate_iota_intercept
    spec$estimate_iota_slope <- parsed$estimate_iota_slope
    spec$n_iota_intercept <- parsed$n_iota_intercept
    spec$n_iota_slope <- parsed$n_iota_slope
    spec$iota_nodes <- parsed$iota_nodes
  } else {
    if (isTRUE(spec$estimate_iota_intercept) || isTRUE(spec$estimate_iota_slope) ||
        !is.null(spec$intercept) || !is.null(spec$slope)) {
      cli::cli_abort(c(
        x = "Fit-only affine-shift options are supported only for functional association transforms.",
        i = "Use {.arg intercept}/{.arg slope} only inside formulas such as {.code ~ expit(x, intercept = TRUE, slope = TRUE)}."
      ))
    }
    spec$estimate_iota_intercept <- FALSE
    spec$estimate_iota_slope <- FALSE
    spec$n_iota_intercept <- 0L
    spec$n_iota_slope <- 0L
    spec$iota_nodes <- list()
  }
  spec
}

#' @keywords internal
.joinme_priors_ <- function(priors = NULL, validate = TRUE) {
  if (is.null(priors)) {
    priors <- list(
      beta = NULL,
      alpha = NULL,
      iota = NULL,
      marker = NULL,
      marker_weight = NULL,
      vcov_sd = NULL,
      vcov_corr = NULL,
      lkj = NULL,
      class_probability = 1,
      class_regression = NULL
    )
  }
  priors <- if (inherits(priors, "joinme_priors")) unclass(priors) else priors

  if (!is.list(priors)) {
    cli::cli_abort(c(
      x = "{.arg priors} must be a named list or a {.fn joinme_priors} object.",
      i = "Example: jm_prior(beta = prior_normal(scale = 2.5), alpha = prior_student_t(df = 6), lkj = prior_lkj(2))."
    ))
  }
  if (length(priors) > 0 && (is.null(names(priors)) || any(names(priors) %in% c("", NA_character_)))) {
    cli::cli_abort(c(
      x = "{.arg priors} must be a named list.",
      i = "Allowed components are beta, alpha, iota, marker, marker_weight, vcov_sd, vcov_corr, lkj, class_probability, and class_regression."
    ))
  }

  allowed <- c(
    "beta",
    "alpha",
    "iota",
    "marker",
    "marker_weight",
    "vcov_sd",
    "vcov_corr",
    "lkj",
    "class_probability",
    "class_regression"
  )
  bad <- setdiff(names(priors), allowed)
  if (length(bad) > 0) {
    cli::cli_abort(c(
      x = "Unknown prior component(s): {paste(bad, collapse = ', ')}.",
      i = "Allowed components: beta, alpha, iota, marker, marker_weight, vcov_sd, vcov_corr, lkj, class_probability, class_regression."
    ))
  }

  out <- list(
    beta = .normalise_joinme_prior_component(
      priors$beta,
      name = "beta",
      default = prior_student_t(df = 6, scale = 2)
    ),
    alpha = .normalise_joinme_prior_component(
      priors$alpha,
      name = "alpha",
      default = prior_student_t(df = 6, scale = 1)
    ),
    iota = .normalise_joinme_prior_component(
      priors$iota,
      name = "iota",
      default = prior_student_t(df = 6, scale = 1)
    ),
    marker = .normalise_joinme_prior_component(
      priors$marker,
      name = "marker",
      default = prior_normal(),
      fixed_unit_scale = TRUE
    ),
    marker_weight = .normalise_joinme_prior_component(
      priors$marker_weight,
      name = "marker_weight",
      default = prior_student_t(df = 6),
      fixed_unit_scale = TRUE
    ),
    vcov_sd = .normalise_joinme_prior_component(
      priors$vcov_sd,
      name = "vcov_sd",
      default = prior_student_t(df = 6, scale = 1)
    ),
    vcov_corr = .normalise_joinme_prior_component(
      priors$vcov_corr,
      name = "vcov_corr",
      default = prior_student_t(df = 6, scale = 1)
    ),
    lkj = .normalise_joinme_lkj_prior(priors$lkj),
    class_probability = priors$class_probability %||% 1,
    class_regression = if (is.numeric(priors$class_regression)) {
      prior_normal(scale = priors$class_regression)
    } else {
      .normalise_joinme_prior_component(
        priors$class_regression,
        name = "class_regression",
        default = prior_normal()
      )
    }
  )

  if (
    !is.numeric(out$class_probability) ||
      length(out$class_probability) < 1L ||
      any(!is.finite(out$class_probability)) ||
      any(out$class_probability <= 0)
  ) {
    cli::cli_abort(
      "{.arg class_probability} must contain positive finite Dirichlet concentrations."
    )
  }
  out$class_probability <- as.numeric(out$class_probability)
  if (isTRUE(validate)) {
    .build_priors(
      beta_prior = out$beta,
      alpha_prior = out$alpha,
      iota_prior = out$iota,
      marker_prior = out$marker,
      marker_weight_prior = out$marker_weight,
      vcov_sd_prior = out$vcov_sd,
      vcov_corr_prior = out$vcov_corr,
      lkj_prior = out$lkj
    )
  }

  structure(out, class = c("joinme_priors", "list"))
}

#' @keywords internal
.normalise_joinme_prior_component <- function(
  x,
  name,
  default,
  fixed_unit_scale = FALSE
) {
  if (is.null(x)) return(default)
  if (inherits(x, "joinme_prior_spec")) {
    output <- x
  } else if (is.numeric(x)) {
    if (length(x) < 1L || any(!is.finite(x)) || any(x <= 0)) {
      cli::cli_abort(c(
        x = "{.arg {name}} must contain positive finite scale values.",
        i = "Prefer {.code {name} = prior_student_t(df = 6, scale = 1)}."
      ))
    }
    output <- prior_student_t(df = 6, scale = x)
  } else if (is.list(x) && !is.null(x$scale %||% x$sd)) {
    scale <- x$scale %||% x$sd
    mu <- x$mu %||% x$location %||% 0
    output <- prior_student_t(df = x$df %||% 6, mu = mu, scale = scale)
  } else {
    cli::cli_abort(c(
      x = "Prior component {.arg {name}} is not a recognised JoiNMe prior declaration.",
      i = "Use {.fn prior_normal}, {.fn prior_student_t}, {.fn prior_laplace}, or {.fn prior_horseshoe}."
    ))
  }

  if (
    isTRUE(fixed_unit_scale) &&
      (
        !identical(as.numeric(output$mu), 0) ||
          !identical(as.numeric(output$scale), 1) ||
          (
            identical(output$family, "horseshoe") &&
              (!identical(output$global_scale, 1) || !identical(output$slab_scale, 2))
          )
      )
  ) {
    cli::cli_abort(c(
      x = "{.arg {name}} is a standardised latent block and cannot accept a location or scale.",
      i = "Choose only its family, for example {.code {name} = prior_laplace()} or {.code {name} = prior_horseshoe()}.",
      i = "Its location is fixed at zero and its ordinary scale is fixed at one."
    ))
  }
  output
}

#' @keywords internal
#' @noRd
.normalise_joinme_lkj_prior <- function(x) {
  if (is.null(x)) return(prior_lkj(1))
  if (is.numeric(x)) {
    return(prior_lkj(x))
  }
  if (inherits(x, "joinme_lkj_prior")) {
    return(x)
  }
  cli::cli_abort(c(
    x = "{.arg lkj} must be created with {.fn prior_lkj}.",
    i = "A positive numeric concentration remains accepted for compatibility."
  ))
}
