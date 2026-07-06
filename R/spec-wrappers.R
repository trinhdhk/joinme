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
#' @export
#'
#' @examples
#' tf <- joinme_tf(
#'   cv_total = "identity",
#'   corr = ~ -x,
#'   vcov = "identity",
#'   cv_marker = list(type = "pwlin", x = c(-1, 0, 1), y = c(0.2, 1, 0.2))
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
      return(paste0(type, ", points=", length(spec$x %||% numeric(0))))
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

#' Declare priors with validation
#'
#' @description
#' `joinme_priors()` is a user-facing wrapper for prior declarations passed to
#' `joinme()`. It validates the exposed prior components and returns a
#' structured object that can be supplied directly as the `priors` argument.
#'
#' Exposed components are:
#' - `beta`: longitudinal/survival regression prior scale(s)
#' - `alpha`: association prior scale
#' - `iota`: fit-only affine-shift prior scale for functional association transforms
#' - `lkj`: LKJ concentration parameter
#'
#' @param beta Beta prior specification. Use either a numeric scale (or vector of
#'   scales) or a list with `scale`/`sd`.
#' @param alpha Alpha prior specification. Use either a numeric scale or a list
#'   with `scale`/`sd`.
#' @param iota Iota prior specification for fit-only affine-shift intercept and
#'   slope parameters in functional association transforms. Use either a numeric
#'   scale or a list with `scale`/`sd`.
#' @param lkj LKJ concentration parameter.
#' @param .validate Logical; if `TRUE` (default), validate the resulting prior
#'   declarations immediately.
#'
#' @return An object of class `joinme_priors`.
#' @export
#'
#' @examples
#' pri <- joinme_priors(
#'   beta = list(scale = 2.5),
#'   alpha = list(scale = 1.0),
#'   iota = list(scale = 1.0),
#'   lkj = 2
#' )
#' print(pri)
joinme_priors <- function(beta = NULL, alpha = NULL, iota = NULL, lkj = NULL, .validate = TRUE) {
  .normalise_joinme_priors_input(
    list(beta = beta, alpha = alpha, iota = iota, lkj = lkj),
    validate = .validate
  )
}

#' @export
print.joinme_priors <- function(x, ...) {
  priors <- unclass(x)
  beta_txt <- if (is.null(priors$beta)) "default" else {
    beta_scale <- if (is.list(priors$beta)) priors$beta$scale %||% priors$beta$sd else priors$beta
    paste(beta_scale, collapse = ", ")
  }
  alpha_txt <- if (is.null(priors$alpha)) "default" else {
    alpha_scale <- if (is.list(priors$alpha)) priors$alpha$scale %||% priors$alpha$sd else priors$alpha
    paste(alpha_scale, collapse = ", ")
  }
  iota_txt <- if (is.null(priors$iota)) "default" else {
    iota_scale <- if (is.list(priors$iota)) priors$iota$scale %||% priors$iota$sd else priors$iota
    paste(iota_scale, collapse = ", ")
  }
  lkj_txt <- if (is.null(priors$lkj)) "default" else paste(priors$lkj, collapse = ", ")

  tbl <- data.frame(
    component = c("beta", "alpha", "iota", "lkj"),
    value = c(beta_txt, alpha_txt, iota_txt, lkj_txt),
    stringsAsFactors = FALSE
  )
  cat("Prior specification for Joint Mixed Effects model\n")
  print(tbl, row.names = FALSE)
  invisible(x)
}

#' Build named plotting conditions
#'
#' @description
#' `make_conditions()` is a thin wrapper around [brms::make_conditions()] so the
#' resulting condition tables can be passed directly to `plot.JoiNMeFit()`
#' through its `condition` argument. This keeps the conditioning workflow close
#' to `brms::conditional_effects()` while staying within the JoiNMe interface.
#'
#' @param x A data frame containing the baseline covariates used to define the
#'   conditioning rows.
#' @param ... Additional arguments passed to [brms::make_conditions()].
#'
#' @return A data frame with one row per requested condition.
#' @export
#' @importFrom brms make_conditions
make_conditions <- function(x, ...) {
  brms::make_conditions(x = x, ...)
}

.fit_affine_shift_supported_functions <- function() {
  c(
    "log", "exp", "sqrt", "inv_logit", "sigmoid", "expit", "softmax",
    "logit", "rec", "sin", "cos", "tan", "abs", "sinh", "cosh", "tanh",
    "asinh", "acosh", "atanh", "softplus", "log1p_exp", "cbrt", "probit", "power"
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
.normalise_joinme_priors_input <- function(priors = NULL, validate = TRUE) {
  if (is.null(priors)) {
    priors <- list(beta = NULL, alpha = NULL, iota = NULL, lkj = NULL)
  }
  priors <- if (inherits(priors, "joinme_priors")) unclass(priors) else priors

  if (!is.list(priors)) {
    cli::cli_abort(c(
      x = "{.arg priors} must be a named list or a {.fn joinme_priors} object.",
      i = "Example: joinme_priors(beta = list(scale = 2.5), alpha = list(scale = 1), lkj = 2)."
    ))
  }
  if (length(priors) > 0 && (is.null(names(priors)) || any(names(priors) %in% c("", NA_character_)))) {
    cli::cli_abort(c(
      x = "{.arg priors} must be a named list.",
      i = "Allowed components are beta, alpha, and lkj."
    ))
  }

  allowed <- c("beta", "alpha", "iota", "lkj")
  bad <- setdiff(names(priors), allowed)
  if (length(bad) > 0) {
    cli::cli_abort(c(
      x = "Unknown prior component(s): {paste(bad, collapse = ', ')}.",
      i = "Allowed components: beta, alpha, lkj."
    ))
  }

  out <- list(
    beta = priors$beta %||% NULL,
    alpha = priors$alpha %||% NULL,
    iota = priors$iota %||% NULL,
    lkj = priors$lkj %||% NULL
  )

  .validate_JoiNMe_prior_component(out$beta, "beta")
  .validate_JoiNMe_prior_component(out$alpha, "alpha")
  .validate_JoiNMe_prior_component(out$iota, "iota")
  if (!is.null(out$lkj)) {
    if (!is.numeric(out$lkj) || length(out$lkj) != 1L || !is.finite(out$lkj) || out$lkj <= 0) {
      cli::cli_abort(c(
        x = "{.arg lkj} must be a positive numeric scalar.",
        i = "Example: joinme_priors(lkj = 2)."
      ))
    }
    out$lkj <- as.numeric(out$lkj)
  }

  if (isTRUE(validate)) {
    .build_priors(beta_prior = out$beta, alpha_prior = out$alpha, iota_prior = out$iota, lkj_prior = out$lkj)
  }

  structure(out, class = c("joinme_priors", "list"))
}

#' @keywords internal
.validate_JoiNMe_prior_component <- function(x, name) {
  if (is.null(x)) return(invisible(NULL))
  if (is.numeric(x)) {
    if (length(x) < 1L || any(!is.finite(x)) || any(x <= 0)) {
      cli::cli_abort(c(
        x = "{.arg {name}} must contain positive finite scale values.",
        i = "Example: joinme_priors({name} = list(scale = 2.5))."
      ))
    }
    return(invisible(NULL))
  }
  if (is.list(x)) {
    scale <- x$scale %||% x$sd
    if (is.null(scale)) {
      cli::cli_abort(c(
        x = "Prior component {.arg {name}} must provide {.arg scale} or {.arg sd} when declared as a list.",
        i = "Example: joinme_priors({name} = list(scale = 1.5))."
      ))
    }
    if (!is.numeric(scale) || length(scale) < 1L || any(!is.finite(scale)) || any(scale <= 0)) {
      cli::cli_abort(c(
        x = "Prior component {.arg {name}} must use positive finite scale values.",
        i = "Example: joinme_priors({name} = list(scale = 1.5))."
      ))
    }
    return(invisible(NULL))
  }

  cli::cli_abort(c(
    x = "Prior component {.arg {name}} must be numeric or a list with {.arg scale}/{.arg sd}.",
    i = "Example: joinme_priors({name} = list(scale = 1.5))."
  ))
}
