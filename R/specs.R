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
#' @return A `joinme_prior_spec` object for use inside [jm_priors()].
#' @export
prior_normal <- function(mu = 0, scale = 1) {
  make_prior_dist("normal", mu = mu, scale = scale)
}

#' Declare a Student-t coefficient prior
#'
#' @inheritParams prior_normal
#' @param df Positive fixed degrees of freedom. A single value is shared by the
#'   complete parameter block. Marker-weight departures instead use the bare
#'   family name `"student_t"`, whose degrees of freedom are fitted.
#'
#' @return A `joinme_prior_spec` object for use inside [jm_priors()].
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
#' @return A `joinme_prior_spec` object for use inside [jm_priors()].
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
#' whilst the local scales permit coefficient-specific escape from regularisation.
#'
#' @param df Positive fixed degrees of freedom for local regularisation scales.
#' @param global_df Positive fixed degrees of freedom for the global scale.
#' @param global_scale Positive global scale. Smaller values express stronger
#'   prior sparsity.
#' @param slab_df Positive fixed degrees of freedom for the regularising slab.
#' @param slab_scale Positive scale of the regularising slab.
#'
#' @return A `joinme_prior_spec` object for use inside [jm_priors()].
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
#' @return A `joinme_lkj_prior` object for use inside [jm_priors()].
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

#' Declare priors by scientific model component
#'
#' @description
#' `jm_priors()` names priors by the part of the statistical model they govern,
#' rather than by internal coefficient letters. `intercept` and `slope` are
#' global fallbacks. A component-specific declaration replaces only the named
#' role, leaving the other role to inherit its global fallback.
#'
#' Components with ordinary regression roles are `longitudinal`, `vcov`, and
#' `functional`. `survival` and `assoc` are slope-only because the baseline
#' hazard supplies the event-process intercept. `marker_weights` is special:
#' its `offset` gives the known marker-specific contribution, `shared` selects
#' one common set or term-specific sets, `intercept` governs each fitted common
#' weight location, and `family` selects only the centred, unit-scale
#' distribution of marker-specific departures.
#' Names appearing on the left-hand side of `formulaDist` (`sigma`, `nu`,
#' `phi`, `alpha`, `kappa`, and `tau`) may be supplied through `...`, each with
#' its own intercept and slope priors. Bracket selectors refine a declaration
#' to a family-scoped coefficient block, for example
#' `` `sigma[family='student']` = list(slope = prior_normal()) ``. A marker
#' selector such as `` `sigma[marker='y']` `` is accepted when that marker
#' uniquely identifies one family-scoped coefficient block. If several markers
#' share the block, use its family selector because the fitted coefficient is
#' shared and cannot receive different marker-specific priors.
#'
#' @param intercept,slope Global coefficient priors. These arguments must be
#'   created with a `prior_*()` function. Fixed numerical coefficients belong
#'   to [jm_truth()] and are rejected here because a constant is not a prior
#'   distribution.
#' @param longitudinal A bare prior declaration in complete longitudinal
#'   model-matrix order, or a named list containing `intercept` and/or `slope`
#'   declarations in their separate role-wise orders.
#' @param survival Slope prior for event-model covariates. It may be supplied
#'   directly or as `list(slope = ...)`.
#' @param baseline Intercept and slope priors for a formula-based baseline
#'   hazard.
#' @param vcov A shared covariance-regression declaration, an
#'   intercept/slope/latent list, or a list with `sd` and `corr` components.
#'   The last form permits distinct population declarations for marginal
#'   scales and off-diagonal partial correlations. `intercept` and `slope`
#'   govern their regression coefficients; `latent` governs the loading of
#'   unexplained subject variation in the corresponding linear predictor.
#' @param marker_weights A list with `offset`, `intercept`, `family`, and
#'   `shared`. The logical `shared` value determines whether active weighted
#'   association terms use one common marker-weight set (`TRUE`, the default)
#'   or distinct sets (`FALSE`).
#'   `offset` is the known contribution to each marker weight. It accepts one
#'   numeric vector whose entries are either all named by marker or all
#'   unnamed, or a named collection of such vectors for unshared association
#'   terms. Partly named vectors are rejected because their marker alignment is
#'   not statistically defined. `intercept` is the ordinary coefficient prior
#'   for each fitted common marker-weight location and may therefore declare
#'   its own location and scale; when omitted, it inherits the global
#'   `intercept` declaration. It must be created with a `prior_*()` function.
#'   `family`
#'   governs only the centred unit-scale
#'   departure law for marker-specific deviations.
#'   Accepted stochastic values are `"student_t"`, `"normal"`, `"laplace"`, and
#'   `"horseshoe"`. Values `"constant"` and `"none"` are synonymous: the
#'   offset is then used exactly, with no fitted location or departure. The
#'   family name `"student_t"` fits a single degrees-of-freedom
#'   value above two shared by every active weight set under a shifted `Gamma(2, 0.1)`
#'   prior. `prior_*()` objects are not accepted for `family` because
#'   marker-weight departures retain location zero and ordinary scale one.
#'   `marker_weights` itself must be a named list; bare
#'   prior declarations are rejected to avoid confusing the common-location
#'   prior with the departure family. Fixed common locations belong to
#'   [jm_truth()].
#' @param assoc Slope prior for longitudinal--event association coefficients.
#' @param functional Intercept and slope priors for fitted affine shifts inside
#'   functional association transformations.
#' @param marker A family-only prior declaration for standardised marker-level
#'   random effects. Location zero and scale one are enforced.
#' @param ... Named distributional-parameter prior components, such as
#'   `sigma = list(intercept = ..., slope = ...)`,
#'   `` `sigma[family='student']` = list(slope = ...) ``, or
#'   `` `sigma[marker='y']` = list(intercept = ...) ``.
#' @param lkj An object from [prior_lkj()] or a positive numeric concentration.
#' @param class A named list with `baseline_prob`, the positive Dirichlet
#'   concentration for baseline class probabilities; `slope`, the prior
#'   declaration for coefficients from `formulaClass`; and `family`, the
#'   centred unit-scale distribution of latent coordinates within each class.
#'   `family` must be created with [prior_student_t()], [prior_normal()], or
#'   [prior_laplace()]; its location and scale must be zero and one because
#'   class-specific locations and scales are fitted separately. A scalar
#'   `baseline_prob` is repeated over classes; a vector may provide one
#'   concentration per class. The slope remains one prior block because
#'   class-specific formula lists share a reference-class parameterisation.
#'   When omitted, `class$slope` inherits the global `slope` declaration.
#' @param .validate Logical; if `TRUE` (default), validate the resulting prior
#'   declarations immediately.
#'
#' @return An object of class `joinme_priors`.
#' @aliases joinme_prior jm_priors jm_prior
#' @export
#'
#' @examples
#' pri <- jm_priors(
#'   intercept = prior_student_t(df = 6, scale = 2),
#'   slope = prior_normal(scale = 1),
#'   longitudinal = list(slope = prior_normal(scale = 0.5)),
#'   marker_weights = list(
#'     offset = c(marker_a = 0, marker_b = 0.25),
#'     intercept = prior_normal(scale = 0.75),
#'     family = "student_t",
#'     shared = TRUE
#'   ),
#'   assoc = prior_student_t(df = 4, scale = 1),
#'   sigma = list(intercept = prior_normal(), slope = prior_normal(scale = 0.5)),
#'   `sigma[family='student']` = list(slope = prior_student_t(df = 4)),
#'   lkj = prior_lkj(2),
#'   class = list(
#'     baseline_prob = c(2, 2, 2),
#'     slope = prior_normal(scale = 1),
#'     family = prior_student_t(df = 6)
#'   )
#' )
#' print(pri)
joinme_priors <- function(
  ...,
  intercept = NULL, # global prior inherited by intercept coefficients
  slope = NULL, # global prior inherited by non-intercept regression coefficients
  longitudinal = NULL, # longitudinal population intercept and slope priors
  survival = NULL, # event-covariate slope prior; the baseline hazard is separate
  baseline = NULL, # formula-based baseline-hazard coefficient priors
  vcov = NULL, # SD and correlation regression intercept and slope priors
  marker_weights = NULL, # marker-weight offset, common-location intercept prior, departure family, and sharing rule
  assoc = NULL, # longitudinal-event association slope prior
  functional = NULL, # affine transformation intercept and slope priors
  marker = NULL, # unit-scale family for standardised marker random effects
  lkj = NULL, # LKJ concentration for correlation matrices
  class = NULL, # baseline probability concentration and formulaClass slope prior
  .validate = TRUE # whether to validate the assembled prior specification
) {
  distributional <- list(...) # formulaDist left-hand-side prior components
  if (length(distributional) > 0L &&
      (is.null(names(distributional)) || any(names(distributional) %in% c("", NA_character_)))) {
    cli::cli_abort(c(
      x = "Distributional priors supplied through {.arg ...} must be named.",
      i = "Use a left-hand-side name from {.arg formulaDist}, for example {.code sigma = list(intercept = prior_normal())}."
    ))
  }
  prior_request <- list(
      intercept = intercept,
      slope = slope,
      longitudinal = longitudinal,
      survival = survival,
      baseline = baseline,
      vcov = vcov,
      marker_weights = marker_weights,
      assoc = assoc,
      functional = functional,
      marker = marker,
      distributional = distributional,
      lkj = lkj,
      class = class
    ) # analyst's probability distributions grouped by scientific component

  # Every coefficient declaration in a fitting prior must be a probability
  # distribution. Numerical values are meaningful for a truth declaration but
  # would silently lose that meaning during fitting. Validate them before the
  # ordinary prior normaliser can interpret a number as an historical scale
  # shorthand.
  .assert_joinme_prior_distributions(prior_request)
  .joinme_priors_(prior_request, validate = .validate)
}

#' @rdname joinme_priors
#' @export
jm_priors <- joinme_priors

#' Declare population truths for JoiNMe simulation
#'
#' @description
#' `jm_truth()` describes the population quantities used once to generate a
#' complete data set. Its coefficient hierarchy follows [jm_priors()]: global
#' `intercept` and `slope` declarations provide fallbacks, whilst named model
#' components may replace either role.
#'
#' A finite numeric vector fixes the corresponding population coefficient. A
#' `prior_*()` declaration draws one coefficient vector once at the beginning
#' of a simulation. The realised vector is then held constant for every subject,
#' marker and observation and is recorded in the returned simulation's `truth`.
#' Thus, a distribution supplied here describes variation between simulated
#' data sets, not variation between observations within one data set.
#'
#' The declaration also contains the baseline hazard and the covariance
#' parameters for ordinary random-effect blocks. If a random-effect standard
#' deviation is omitted, one value per coefficient is drawn from
#' `Exponential(1)`. If its correlation matrix is omitted, one matrix is drawn
#' from the LKJ distribution selected by `lkj`. These are the same population
#' distributions used by the Stan fitting model.
#'
#' @param intercept,slope Global fixed values or `prior_*()` generating
#'   distributions for intercept and non-intercept population coefficients.
#' @param longitudinal Fixed values or generating distributions for the
#'   longitudinal population regression. A named list may contain `intercept`
#'   and `slope` separately.
#' @param survival Fixed values or a generating distribution for event-model
#'   covariate slopes.
#' @param baseline Fixed values or generating distributions for coefficients
#'   of a formula baseline hazard.
#' @param vcov Covariance-regression truth. It may be shared between the
#'   marginal-standard-deviation and correlation regressions, or supplied as
#'   `list(sd = ..., corr = ...)`. Each part accepts `intercept`, `slope`, and
#'   `latent` declarations.
#' @param marker_weights A named list containing `offset`, `intercept`, `family`,
#'   and `shared`. `offset` is a fixed marker-specific contribution;
#'   `intercept` is a fixed or once-drawn common marker-weight location;
#'   `family` governs centred unit-scale marker departures; and `shared`
#'   determines whether weighted association terms share one set.
#' @param assoc_coef Fixed values or a generating distribution for the active
#'   longitudinal--event association coefficients, in their fitted order.
#' @param functional Fixed values or generating distributions for affine
#'   transformation intercepts and slopes.
#' @param marker A `prior_*()` family declaration for standardised marker-level
#'   effects. Its location and ordinary scale must remain zero and one.
#' @param ... Named distributional-parameter truths corresponding to the
#'   left-hand sides of `formulaDist`, including family or marker selectors.
#' @param lkj An object created by [prior_lkj()]. Its concentration governs each
#'   random-effect correlation matrix drawn when `re_params` omits `corr`.
#' @param class A named list containing `baseline_prob`, `slope`, and `family`,
#'   following the class roles accepted by [jm_priors()]. The slope may be fixed
#'   or drawn once for each simulation. The family is a probability declaration
#'   governing every standardised latent coordinate drawn conditionally on its
#'   class; it is not itself a population coefficient draw.
#' @param basehaz Baseline-hazard truth. Supply a hazard function, a formula in
#'   `time`, a character family name, or a named list such as
#'   `list(type = "weibull", shape = 1.4, scale = 6)`.
#' @param re_params Named random-effect covariance declarations for `id`,
#'   `marker`, and distributional regressions under `dist`. Each ordinary block
#'   accepts `sd` and `corr`; either may be omitted and drawn once from the
#'   fitted model's corresponding population distribution.
#' @param .validate Logical; if `TRUE`, check the declaration immediately.
#'
#' @return An object of class `joinme_truth` for the `truth` argument of
#'   [simulate_joinme()] or [simulate_joinme_mix()].
#' @aliases joinme_truth
#' @export
#'
#' @examples
#' generating_truth <- jm_truth(
#'   longitudinal = list(
#'     intercept = 0,
#'     slope = prior_normal(mu = 0.5, scale = 0.2)
#'   ),
#'   survival = c(treatment = -0.4),
#'   assoc_coef = c(cv_mean = 0.3),
#'   basehaz = list(type = "weibull", shape = 1.2, scale = 7),
#'   re_params = list(id = list(sd = NULL, corr = NULL)),
#'   lkj = prior_lkj(2)
#' )
joinme_truth <- function(
  ...,
  intercept = NULL, # global intercept truth or its between-simulation generating distribution
  slope = NULL, # global slope truth or its between-simulation generating distribution
  longitudinal = NULL, # longitudinal population coefficient truths
  survival = NULL, # event-regression population coefficient truths
  baseline = NULL, # formula baseline-hazard coefficient truths
  vcov = NULL, # covariance-regression population truths
  marker_weights = NULL, # marker-weight offsets, locations, departure family, and sharing rule
  assoc_coef = NULL, # active longitudinal--event association coefficient truths
  functional = NULL, # affine transformation coefficient truths
  marker = NULL, # standardised marker-effect family
  lkj = prior_lkj(1), # correlation-generating LKJ distribution
  class = NULL, # class-membership probability and regression truths
  basehaz = list(type = "weibull", shape = 1.4, scale = 6.0), # baseline-hazard function or family declaration
  re_params = list(
    id = list(sd = NULL, corr = NULL),
    marker = list(sd = NULL, corr = NULL),
    dist = list()
  ), # ordinary random-effect covariance truths
  .validate = TRUE # whether to check the assembled truth declaration
) {
  distributional <- list(...) # distributional population truths named by formulaDist left-hand sides
  if (length(distributional) > 0L &&
      (is.null(names(distributional)) || any(names(distributional) %in% c("", NA_character_)))) {
    cli::cli_abort(c(
      x = "Distributional truths supplied through {.arg ...} must be named.",
      i = "Use a left-hand-side name from {.arg formulaDist}, for example {.code sigma = list(intercept = 0)}."
    ))
  }

  # Retain the analyst's exact fixed-or-random declarations. The simulator
  # resolves their dimensions only after the model matrices are known, which
  # permits named vectors to follow the fitted coefficient order exactly.
  truth_request <- list(
    intercept = intercept,
    slope = slope,
    longitudinal = longitudinal,
    survival = survival,
    baseline = baseline,
    vcov = vcov,
    marker_weights = marker_weights,
    assoc_coef = assoc_coef,
    functional = functional,
    marker = marker,
    distributional = distributional,
    lkj = lkj,
    class = class,
    basehaz = basehaz,
    re_params = re_params
  ) # complete data-generating declaration before formula dimensions are available

  if (!is.list(re_params) || is.null(names(re_params))) {
    cli::cli_abort("{.arg re_params} must be a named list.")
  }
  unknown_random_effect_blocks <- setdiff(names(re_params), c("id", "marker", "dist")) # unrecognised random-effect truth components
  if (length(unknown_random_effect_blocks) > 0L) {
    cli::cli_abort("Unknown {.arg re_params} component{?s}: {.field {unknown_random_effect_blocks}}.")
  }
  if (!(is.function(basehaz) || inherits(basehaz, "formula") || is.character(basehaz) || is.list(basehaz))) {
    cli::cli_abort("{.arg basehaz} must be a function, formula, family name, or named list.")
  }
  if (!inherits(lkj, "joinme_lkj_prior")) {
    cli::cli_abort("{.arg lkj} must be created with {.fn prior_lkj}.")
  }
  marker_family_declaration <- if (
    is.list(marker) && !inherits(marker, "joinme_prior_spec")
  ) marker$family else marker # unit-scale family declaration for standardised marker effects
  if (!is.null(marker_family_declaration) &&
      !inherits(marker_family_declaration, "joinme_prior_spec")) {
    cli::cli_abort(c(
      x = "{.arg marker} describes a centred unit-scale random-effect family and cannot be a fixed number.",
      i = "Use a {.fn prior_student_t}, {.fn prior_normal}, {.fn prior_laplace}, or {.fn prior_horseshoe} declaration."
    ))
  }

  # Construct the fitting-prior counterpart only from probability
  # distributions. Fixed truths are replaced by the ordinary component
  # defaults, whilst structural marker-weight settings remain identical so a
  # recovery fit represents the data-generating marker feature correctly.
  fitting_request <- truth_request[c(
    "intercept", "slope", "longitudinal", "survival", "baseline", "vcov",
    "marker_weights", "functional", "marker", "distributional", "lkj", "class"
  )]
  fitting_request$assoc <- truth_request$assoc_coef
  fitting_request <- .remove_joinme_truth_constants(fitting_request)
  fitting_priors <- .joinme_priors_(fitting_request, validate = .validate) # probability distributions retained for an optional recovery fit

  structure(
    truth_request,
    class = c("joinme_truth", "list"),
    fitting_priors = fitting_priors
  )
}

#' @rdname joinme_truth
#' @export
jm_truth <- joinme_truth

#' Print a JoiNMe simulation-truth declaration
#'
#' @param x A `joinme_truth` object.
#' @param ... Unused.
#'
#' @return `x`, invisibly.
#' @export
print.joinme_truth <- function(x, ...) {
  cat("Truth declaration for JoiNMe simulation\n")
  populated <- names(x)[vapply(x, function(value) !is.null(value) && length(value) > 0L, logical(1))] # scientific components explicitly or structurally represented
  cat("Components: ", paste(populated, collapse = ", "), "\n", sep = "")
  cat("Population coefficients are fixed once per simulated data set.\n")
  invisible(x)
}

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
  role_rows <- function(component, specification) {
    roles <- intersect(c("intercept", "slope", "latent", "family"), names(specification))
    data.frame(
      component = paste0(component, ".", roles),
      value = vapply(specification[roles], describe_prior, character(1)),
      stringsAsFactors = FALSE
    )
  }
  distributional_rows <- function(parameter, specification) {
    rows <- list(role_rows(parameter, specification)) # parameter-wide fallback roles
    if (length(specification$by_family) > 0L) {
      rows <- c(rows, lapply(names(specification$by_family), function(family_name) {
        role_rows(
          paste0(parameter, "[family=", family_name, "]"),
          specification$by_family[[family_name]]
        )
      }))
    }
    if (length(specification$by_marker) > 0L) {
      rows <- c(rows, lapply(names(specification$by_marker), function(marker_name) {
        role_rows(
          paste0(parameter, "[marker=", marker_name, "]"),
          specification$by_marker[[marker_name]]
        )
      }))
    }
    do.call(rbind, rows)
  }
  describe_marker_weight_offset <- function(offset) {
    if (is.null(offset)) return("default")
    describe_vector <- function(values) {
      value_names <- names(values)
      if (is.null(value_names) || !all(nzchar(value_names))) {
        return(paste(as.numeric(values), collapse = ", "))
      }
      paste0(value_names, " = ", as.numeric(values), collapse = ", ")
    }
    if (is.numeric(offset) && is.null(dim(offset))) return(describe_vector(offset))
    if (is.list(offset) && !is.data.frame(offset)) {
      return(paste(vapply(names(offset), function(term_name) {
        paste0(term_name, ": [", describe_vector(offset[[term_name]]), "]")
      }, character(1)), collapse = "; "))
    }
    paste(apply(as.matrix(offset), 1L, function(values) {
      paste0("[", paste(as.numeric(values), collapse = ", "), "]")
    }), collapse = "; ")
  }
  marker_weight_rows <- role_rows("marker_weights", priors$marker_weights) # common-location prior and standardised departure family
  marker_weight_rows <- rbind(
    data.frame(
      component = "marker_weights.offset",
      value = describe_marker_weight_offset(priors$marker_weights$offset),
      stringsAsFactors = FALSE
    ),
    data.frame(
      component = "marker_weights.shared",
      value = as.character(priors$marker_weights$shared),
      stringsAsFactors = FALSE
    ),
    marker_weight_rows
  )
  if (identical(priors$marker_weights$family$family, "student_t")) {
    family_row <- marker_weight_rows$component == "marker_weights.family" # Student-t departure row requiring fixed-versus-fitted df interpretation
    marker_weight_rows$value[family_row] <- "student_t(fitted df: 2 + Gamma(2, 0.1); location = 0; scale = 1)"
  } else if (identical(priors$marker_weights$family$family, "constant")) {
    marker_weight_rows$value[marker_weight_rows$component == "marker_weights.intercept"] <- "not fitted"
    marker_weight_rows$value[marker_weight_rows$component == "marker_weights.family"] <- "constant (offset used exactly)"
  }
  tbl <- do.call(rbind, c(
    list(role_rows("global", priors$global)),
    lapply(c("longitudinal", "survival", "baseline", "assoc", "functional"), function(component) {
      role_rows(component, priors[[component]])
    }),
    list(
      role_rows("vcov.sd", priors$vcov$sd),
      role_rows("vcov.corr", priors$vcov$corr),
      marker_weight_rows,
      role_rows("marker", priors$marker),
      do.call(rbind, lapply(names(priors$distributional), function(parameter) {
        distributional_rows(parameter, priors$distributional[[parameter]])
      })),
      data.frame(component = "lkj", value = describe_prior(priors$lkj), stringsAsFactors = FALSE),
      data.frame(component = "class.baseline_prob", value = paste(priors$class$baseline_prob, collapse = ", "), stringsAsFactors = FALSE),
      data.frame(component = "class.slope", value = describe_prior(priors$class$slope), stringsAsFactors = FALSE),
      data.frame(component = "class.family", value = describe_prior(priors$class$family), stringsAsFactors = FALSE)
    )
  ))
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
#' [conditional_effects.JoiNMeFit()] or [conditional_contrast()]. Conditioning
#' belongs to the conditional-effects or conditional-contrast estimand; neither
#' [predict.JoiNMeFit()] nor `plot.JoiNMeFit()` accepts a `condition` argument.
#'
#' @param x A data frame containing the baseline covariates used to define the
#'   conditioning rows.
#' @param ... Additional arguments passed to [brms::make_conditions()].
#'
#' @return A data frame with one row per requested condition and a `cond__`
#'   column containing the profile labels used in conditional-effects facets.
#'
#' @seealso [conditional_effects.JoiNMeFit()], [conditional_contrast()],
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
  if (inherits(priors, "joinme_priors")) {
    return(priors)
  }
  if (is.null(priors)) {
    priors <- list(
      intercept = NULL,
      slope = NULL,
      longitudinal = NULL,
      survival = NULL,
      baseline = NULL,
      vcov = NULL,
      marker_weights = NULL,
      assoc = NULL,
      functional = NULL,
      marker = NULL,
      distributional = list(),
      lkj = NULL,
      class = NULL
    )
  }
  if (!is.list(priors)) {
    cli::cli_abort(c(
      x = "{.arg priors} must be a named list or a {.fn joinme_priors} object.",
      i = "Example: jm_priors(intercept = prior_normal(scale = 2.5), slope = prior_student_t(df = 6), lkj = prior_lkj(2))."
    ))
  }
  if (length(priors) > 0 && (is.null(names(priors)) || any(names(priors) %in% c("", NA_character_)))) {
    cli::cli_abort(c(
      x = "{.arg priors} must be a named list.",
      i = "Prior components must be named by their scientific model role."
    ))
  }

  allowed <- c(
    "intercept",
    "slope",
    "longitudinal",
    "survival",
    "baseline",
    "vcov",
    "marker_weights",
    "assoc",
    "functional",
    "marker",
    "distributional",
    "lkj",
    "class"
  )
  bad <- setdiff(names(priors), allowed)
  if (length(bad) > 0) {
    cli::cli_abort(c(
      x = "Unknown prior component(s): {paste(bad, collapse = ', ')}.",
      i = "Allowed components are intercept, slope, longitudinal, survival, baseline, vcov, marker_weights, assoc, functional, marker, distributional, lkj, and class."
    ))
  }

  complete_prior_counter <- 0L # unique identifier for a bare prior spanning every role of one component
  normalise_roles <- function(x, component, roles = c("intercept", "slope"), defaults) {
    if (is.null(x)) x <- list()
    if (inherits(x, "joinme_prior_spec")) {
      complete_prior_counter <<- complete_prior_counter + 1L
      complete_prior_id <- paste0(component, "#", complete_prior_counter) # groups role-specific views of one complete coefficient prior
      x <- stats::setNames(lapply(roles, function(role) {
        declaration <- x # copy retaining the full coefficient-wise hyperparameter vectors
        attr(declaration, "complete_prior_id") <- complete_prior_id
        declaration
      }), roles)
    } else if (is.numeric(x)) {
      x <- stats::setNames(rep(list(x), length(roles)), roles)
    }
    if (!is.list(x) || (length(x) > 0L && (is.null(names(x)) || any(!names(x) %in% roles)))) {
      cli::cli_abort(c(
        x = "Prior component {.arg {component}} must be a prior declaration or a list using {.field {roles}}.",
        i = "Declare only the statistical roles that should override the global fallback."
      ))
    }
    stats::setNames(lapply(roles, function(role) {
      .normalise_joinme_prior_component(
        x[[role]],
        name = paste0(component, "$", role),
        default = defaults[[role]]
      )
    }), roles)
  }

  global_intercept <- .normalise_joinme_prior_component(
    priors$intercept, "intercept", prior_student_t(df = 6, scale = 2)
  )
  global_slope <- .normalise_joinme_prior_component(
    priors$slope, "slope", prior_student_t(df = 6, scale = 1)
  )
  global_defaults <- list(intercept = global_intercept, slope = global_slope)

  marker_weight_request <- priors$marker_weights %||% list()
  if (inherits(marker_weight_request, "joinme_prior_spec")) {
    cli::cli_abort(c(
      x = "{.arg marker_weights} must be a named list, not a bare prior declaration.",
      i = "Use a named declaration such as {.code marker_weights = list(offset = ..., intercept = ..., family = ..., shared = TRUE)}."
    ))
  }
  if (!is.list(marker_weight_request) ||
      (length(marker_weight_request) > 0L &&
        (is.null(names(marker_weight_request)) || any(!names(marker_weight_request) %in% c("offset", "intercept", "family", "shared"))))) {
    cli::cli_abort("{.arg marker_weights} must contain only {.field offset}, {.field intercept}, {.field family}, and {.field shared}.")
  }
  marker_weight_sets_shared <- marker_weight_request$shared %||% TRUE # whether active weighted association terms use one common marker-weight set
  if (!is.logical(marker_weight_sets_shared) || length(marker_weight_sets_shared) != 1L || is.na(marker_weight_sets_shared)) {
    cli::cli_abort("{.arg marker_weights$shared} must be TRUE or FALSE.")
  }
  marker_weight_offset <- .normalise_marker_weight_offset_request(
    marker_weight_request$offset,
    name = "marker_weights$offset"
  ) # fixed numerical contribution added to every fitted marker weight, or the complete constant weight when family is constant
  # Marker-weight departures are a standardised latent block with fixed
  # location zero and fixed ordinary scale one. Their family is selected by a
  # name because no coefficient-prior location or scale applies to this block.
  marker_weight_family <- .normalise_marker_weight_family(
    marker_weight_request$family,
    "marker_weights$family",
    default = "student_t"
  )

  vcov_request <- priors$vcov %||% list()
  vcov_roles <- c("intercept", "slope", "latent") # covariance-regression population roles used by fitting and simulation
  vcov_defaults <- c(global_defaults, list(latent = global_slope)) # latent loadings use the global slope prior unless declared separately
  if (inherits(vcov_request, "joinme_prior_spec") || is.numeric(vcov_request) ||
      any(names(vcov_request) %in% vcov_roles)) {
    shared_vcov <- normalise_roles(vcov_request, "vcov", roles = vcov_roles, defaults = vcov_defaults)
    vcov_request <- list(sd = shared_vcov, corr = shared_vcov)
  } else {
    if (!is.list(vcov_request) || any(!names(vcov_request) %in% c("sd", "corr"))) {
      cli::cli_abort("{.arg vcov} must use intercept/slope directly or contain {.field sd} and {.field corr} components.")
    }
    vcov_request <- list(
      sd = normalise_roles(vcov_request$sd, "vcov$sd", roles = vcov_roles, defaults = vcov_defaults),
      corr = normalise_roles(vcov_request$corr, "vcov$corr", roles = vcov_roles, defaults = vcov_defaults)
    )
  }

  distributional_names <- c("sigma", "nu", "phi", "alpha", "kappa", "tau") # supported distributional regression parameters
  distributional_request <- priors$distributional %||% list() # unparsed parameter-wide, family-scoped, or marker-scoped declarations
  if (length(distributional_request) > 0L &&
      (is.null(names(distributional_request)) || any(names(distributional_request) %in% c("", NA_character_)))) {
    cli::cli_abort("{.arg distributional} prior components must be named by a {.arg formulaDist} left-hand side.")
  }
  requested_parameter_names <- vapply(names(distributional_request), function(selector_name) {
    .canonical_dist_param(sub("\\[.*$", "", selector_name))
  }, character(1)) # parameter part before an optional bracket selector
  unknown_distributional <- unique(requested_parameter_names[!requested_parameter_names %in% distributional_names])
  if (length(unknown_distributional) > 0L) {
    cli::cli_abort(c(
      x = "Unknown distributional prior component{?s}: {.val {paste(unknown_distributional, collapse = ', ')}}.",
      i = "Use names that may occur on the left-hand side of formulaDist: sigma, nu, phi, alpha, kappa, or tau."
    ))
  }

  parsed_distributional_request <- stats::setNames(lapply(distributional_names, function(parameter) {
    list(default = NULL, by_family = list(), by_marker = list())
  }), distributional_names) # declarations separated by parameter and scientific scope
  for (selector_name in names(distributional_request)) {
    selector <- .parse_dist_selector_text(selector_name, allow_marker = TRUE) # checked parameter, family, or marker selector
    parameter_request <- parsed_distributional_request[[selector$param]] # accumulated declarations for this distributional parameter
    if (is.null(selector$scope_type)) {
      if (!is.null(parameter_request$default)) {
        cli::cli_abort("Distributional prior selector {.val {selector_name}} is declared more than once.")
      }
      parameter_request$default <- distributional_request[[selector_name]]
    } else {
      scope_collection <- if (identical(selector$scope_type, "family")) "by_family" else "by_marker"
      if (!is.null(parameter_request[[scope_collection]][[selector$scope_value]])) {
        cli::cli_abort("Distributional prior selector {.val {selector_name}} is declared more than once.")
      }
      parameter_request[[scope_collection]][[selector$scope_value]] <- distributional_request[[selector_name]]
    }
    parsed_distributional_request[[selector$param]] <- parameter_request
  }

  normalised_distributional <- stats::setNames(lapply(distributional_names, function(parameter) {
    request <- parsed_distributional_request[[parameter]] # all declarations for one distributional parameter
    parameter_default <- normalise_roles(
      request$default,
      parameter,
      defaults = global_defaults
    ) # parameter-wide intercept and slope fallback
    family_overrides <- lapply(names(request$by_family), function(family_name) {
      normalise_roles(
        request$by_family[[family_name]],
        paste0(parameter, "[family=", family_name, "]"),
        defaults = parameter_default
      )
    }) # role priors for explicitly family-scoped coefficients
    names(family_overrides) <- names(request$by_family)
    marker_overrides <- lapply(names(request$by_marker), function(marker_name) {
      normalise_roles(
        request$by_marker[[marker_name]],
        paste0(parameter, "[marker=", marker_name, "]"),
        defaults = parameter_default
      )
    }) # role priors resolved to a unique family-scoped block after marker levels are known
    names(marker_overrides) <- names(request$by_marker)
    c(parameter_default, list(by_family = family_overrides, by_marker = marker_overrides))
  }), distributional_names)

  class_request <- priors$class %||% list() # baseline allocation, formulaClass slope prior and latent-component family
  if (!is.list(class_request) ||
      (length(class_request) > 0L &&
        (is.null(names(class_request)) || any(!names(class_request) %in% c("baseline_prob", "slope", "family"))))) {
    cli::cli_abort(c(
      x = "{.arg class} must be a named list containing only {.field baseline_prob}, {.field slope}, and {.field family}.",
      i = "For example, use {.code class = list(baseline_prob = 1, slope = prior_normal(), family = prior_student_t(df = 6))}."
    ))
  }
  class_baseline_probability <- class_request$baseline_prob %||% 1 # Dirichlet concentration before expansion to n_classes
  if (
    !is.numeric(class_baseline_probability) ||
      length(class_baseline_probability) < 1L ||
      any(!is.finite(class_baseline_probability)) ||
      any(class_baseline_probability <= 0)
  ) {
    cli::cli_abort(
      "{.arg class$baseline_prob} must contain positive finite Dirichlet concentrations."
    )
  }
  class_slope_request <- class_request$slope # common family and hyperparameters for formulaClass slopes
  class_family_request <- class_request$family # location-zero, unit-scale component distribution for selected latent coordinates
  if (!is.null(class_family_request) && !inherits(class_family_request, "joinme_prior_spec")) {
    cli::cli_abort(c(
      x = "{.arg class$family} must be created with {.fn prior_student_t}, {.fn prior_normal}, or {.fn prior_laplace}.",
      i = "The component locations and scales are estimated separately, so fixed numbers are not meaningful here."
    ))
  }
  class_family <- .normalise_joinme_prior_component(
    class_family_request,
    name = "class$family",
    default = prior_student_t(df = 6),
    fixed_unit_scale = TRUE
  ) # validated standard component law before its fitted class-specific location and scale are applied
  if (identical(class_family$family, "horseshoe")) {
    cli::cli_abort(c(
      x = "{.arg class$family} does not accept {.fn prior_horseshoe}.",
      i = "Use {.fn prior_student_t}, {.fn prior_normal}, or {.fn prior_laplace} for a proper location--scale component density."
    ))
  }

  out <- list(
    global = global_defaults,
    longitudinal = normalise_roles(priors$longitudinal, "longitudinal", defaults = global_defaults),
    survival = normalise_roles(priors$survival, "survival", roles = "slope", defaults = list(slope = global_slope)),
    baseline = normalise_roles(priors$baseline, "baseline", defaults = global_defaults),
    vcov = vcov_request,
    marker_weights = list(
      offset = marker_weight_offset,
      shared = isTRUE(marker_weight_sets_shared),
      intercept = .normalise_joinme_prior_component(
        marker_weight_request$intercept, "marker_weights$intercept", global_intercept
      ),
      family = marker_weight_family
    ),
    assoc = normalise_roles(priors$assoc, "assoc", roles = "slope", defaults = list(slope = global_slope)),
    functional = normalise_roles(priors$functional, "functional", defaults = global_defaults),
    distributional = normalised_distributional,
    marker = list(family = .normalise_joinme_prior_component(
      if (is.list(priors$marker) && !inherits(priors$marker, "joinme_prior_spec")) {
        if (length(priors$marker) > 0L &&
            (is.null(names(priors$marker)) || !identical(names(priors$marker), "family"))) {
          cli::cli_abort("{.arg marker} may contain only the named {.field family} declaration.")
        }
        priors$marker$family
      } else priors$marker,
      "marker$family", prior_normal(), fixed_unit_scale = TRUE
    )),
    lkj = .normalise_joinme_lkj_prior(priors$lkj),
    class = list(
      baseline_prob = as.numeric(class_baseline_probability),
      slope = if (is.numeric(class_slope_request)) {
        global_slope
      } else {
        .normalise_joinme_prior_component(
          class_slope_request,
          name = "class$slope",
          default = global_slope
        )
      },
      family = class_family
    )
  )

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
          !identical(as.numeric(output$scale), 1)
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
.normalise_marker_weight_family <- function(x, name, default = "student_t") {
  family_request <- x %||% default # family label for the standardised marker-weight departure
  if (inherits(family_request, "joinme_prior_spec")) {
    cli::cli_abort(c(
      x = "{.arg {name}} accepts a family name rather than a prior declaration.",
      i = "Use {.code family = 'student_t'}, {.code 'normal'}, {.code 'laplace'}, {.code 'horseshoe'}, or {.code 'constant'}."
    ))
  }
  family_name <- family_request
  if (!is.character(family_name) || length(family_name) != 1L || !nzchar(family_name)) {
    cli::cli_abort(c(
      x = "{.arg {name}} must be one family name.",
      i = "Family names are {.val student_t}, {.val normal}, {.val laplace}, {.val horseshoe}, {.val constant}, and {.val none}."
    ))
  }
  family_name <- match.arg(
    tolower(family_name),
    c("student_t", "normal", "laplace", "horseshoe", "constant", "none")
  )
  if (identical(family_name, "none")) family_name <- "constant"

  # A constant family has no stochastic departure and no fitted common
  # location.  The ordinary prior fields are retained only so this declaration
  # has the same predictable shape as the stochastic family declarations.
  if (identical(family_name, "constant")) {
    return(structure(list(
      family = "constant",
      mu = 0,
      scale = 1,
      df = Inf,
      global_df = Inf,
      global_scale = 1,
      slab_df = Inf,
      slab_scale = 1,
      estimate_df = FALSE
    ), class = c("joinme_prior_spec", "list")))
  }

  # A family name selects the standardised departure law. Student-t deliberately
  # leaves df undeclared so fitting uses the shifted-Gamma distributional
  # parameter; other families retain their fixed-unit-scale records.
  output <- switch(
    family_name,
    student_t = prior_student_t(df = 6),
    normal = prior_normal(),
    laplace = prior_laplace(),
    horseshoe = prior_horseshoe()
  )
  output$estimate_df <- identical(family_name, "student_t") # only the family label requests the moving shifted-Gamma degrees of freedom
  output
}

#' Validate a marker-weight offset declaration before marker names are known
#'
#' @description
#' The prior declaration is parsed before the longitudinal data establish the
#' fitted marker order. This helper therefore checks the numerical and naming
#' rules without reordering values. Alignment to the observed marker levels is
#' performed later by `.normalise_marker_weight_vector()`.
#'
#' Every individual offset vector must use one of two complete conventions:
#' either all entries are named, or no entries are named. A partly named vector
#' is ambiguous because unnamed entries cannot be assigned to markers without
#' relying on their incidental position.
#'
#' @param offset `NULL`, a numeric vector, a numeric matrix or data frame, or a
#'   named list of numeric vectors indexed by weighted association term.
#' @param name Name used in diagnostic messages.
#'
#' @return The checked declaration without changing its numerical values or
#'   names.
#' @keywords internal
#' @noRd
.normalise_marker_weight_offset_request <- function(offset, name) {
  check_vector <- function(values, vector_name) {
    if (!is.numeric(values) || length(values) < 1L || any(!is.finite(values))) {
      cli::cli_abort("{.arg {vector_name}} must be a non-empty finite numeric vector.")
    }
    value_names <- names(values)
    if (!is.null(value_names)) {
      named <- !is.na(value_names) & nzchar(value_names)
      if (any(named) && !all(named)) {
        cli::cli_abort(c(
          x = "{.arg {vector_name}} must be wholly named or wholly unnamed.",
          i = "Name every marker entry, or remove every marker name and use fitted marker order."
        ))
      }
      if (all(named) && anyDuplicated(value_names)) {
        cli::cli_abort("{.arg {vector_name}} contains duplicated marker names.")
      }
    }
    values
  }

  if (is.null(offset)) return(NULL)
  if (is.numeric(offset) && is.null(dim(offset))) return(check_vector(offset, name))
  if (is.matrix(offset) || is.data.frame(offset)) {
    if (!all(vapply(offset, is.numeric, logical(1)))) {
      cli::cli_abort("{.arg {name}} must contain only numeric marker offsets.")
    }
    if (any(!is.finite(as.matrix(offset)))) {
      cli::cli_abort("{.arg {name}} must contain only finite marker offsets.")
    }
    column_names <- colnames(offset)
    if (!is.null(column_names)) {
      named_columns <- !is.na(column_names) & nzchar(column_names)
      if (any(named_columns) && !all(named_columns)) {
        cli::cli_abort("The marker columns of {.arg {name}} must be wholly named or wholly unnamed.")
      }
      if (all(named_columns) && anyDuplicated(column_names)) {
        cli::cli_abort("The marker columns of {.arg {name}} contain duplicated names.")
      }
    }
    row_names <- rownames(offset)
    if (is.null(row_names) || anyNA(row_names) || any(!nzchar(row_names)) || anyDuplicated(row_names)) {
      cli::cli_abort(c(
        x = "The rows of {.arg {name}} must be uniquely named.",
        i = "Name each row by its weighted association term."
      ))
    }
    unknown_terms <- setdiff(row_names, .weighted_assoc_term_keys())
    if (length(unknown_terms) > 0L) {
      cli::cli_abort("The rows of {.arg {name}} contain unknown weighted association terms: {.val {paste(unknown_terms, collapse = ', ')}}.")
    }
    return(offset)
  }
  if (is.list(offset)) {
    offset_names <- names(offset)
    if (is.null(offset_names) || anyNA(offset_names) || any(!nzchar(offset_names)) || anyDuplicated(offset_names)) {
      cli::cli_abort(c(
        x = "{.arg {name}} must be a completely named list.",
        i = "Name each vector by its weighted association term."
      ))
    }
    unknown_terms <- setdiff(offset_names, .weighted_assoc_term_keys())
    if (length(unknown_terms) > 0L) {
      cli::cli_abort("{.arg {name}} contains unknown weighted association terms: {.val {paste(unknown_terms, collapse = ', ')}}.")
    }
    return(stats::setNames(lapply(seq_along(offset), function(index) {
      check_vector(offset[[index]], paste0(name, "$", offset_names[[index]]))
    }), offset_names))
  }
  cli::cli_abort(c(
    x = "{.arg {name}} has an unsupported form.",
    i = "Use a numeric vector, or a named collection of numeric vectors for term-specific offsets."
  ))
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
    i = "A positive numeric concentration is also accepted."
  ))
}

#' @keywords internal
#' @noRd
.normalise_marker_weight_family <- function(x, name, default = "student_t") {
  family_request <- x %||% default # family label for the standardised marker-weight departure
  if (inherits(family_request, "joinme_prior_spec")) {
    cli::cli_abort(c(
      x = "{.arg {name}} accepts a family name rather than a prior declaration.",
      i = "Use {.code family = 'student_t'}, {.code 'normal'}, {.code 'laplace'}, {.code 'horseshoe'}, or {.code 'constant'}."
    ))
  }
  family_name <- family_request
  if (!is.character(family_name) || length(family_name) != 1L || !nzchar(family_name)) {
    cli::cli_abort(c(
      x = "{.arg {name}} must be one family name.",
      i = "Family names are {.val student_t}, {.val normal}, {.val laplace}, {.val horseshoe}, {.val constant}, and {.val none}."
    ))
  }
  family_name <- match.arg(
    tolower(family_name),
    c("student_t", "normal", "laplace", "horseshoe", "constant", "none")
  )
  if (identical(family_name, "none")) family_name <- "constant"

  # A constant family has no stochastic departure and no fitted common
  # location.  The ordinary prior fields are retained only so this declaration
  # has the same predictable shape as the stochastic family declarations.
  if (identical(family_name, "constant")) {
    return(structure(list(
      family = "constant",
      mu = 0,
      scale = 1,
      df = Inf,
      global_df = Inf,
      global_scale = 1,
      slab_df = Inf,
      slab_scale = 1,
      estimate_df = FALSE
    ), class = c("joinme_prior_spec", "list")))
  }

  # A family name selects the standardised departure law. Student-t deliberately
  # leaves df undeclared so fitting uses the shifted-Gamma distributional
  # parameter; other families retain their fixed-unit-scale records.
  output <- switch(
    family_name,
    student_t = prior_student_t(df = 6),
    normal = prior_normal(),
    laplace = prior_laplace(),
    horseshoe = prior_horseshoe()
  )
  output$estimate_df <- identical(family_name, "student_t") # only the family label requests the moving shifted-Gamma degrees of freedom
  output
}

#' Validate a marker-weight offset declaration before marker names are known
#'
#' @description
#' The prior declaration is parsed before the longitudinal data establish the
#' fitted marker order. This helper therefore checks the numerical and naming
#' rules without reordering values. Alignment to the observed marker levels is
#' performed later by `.normalise_marker_weight_vector()`.
#'
#' Every individual offset vector must use one of two complete conventions:
#' either all entries are named, or no entries are named. A partly named vector
#' is ambiguous because unnamed entries cannot be assigned to markers without
#' relying on their incidental position.
#'
#' @param offset `NULL`, a numeric vector, a numeric matrix or data frame, or a
#'   named list of numeric vectors indexed by weighted association term.
#' @param name Name used in diagnostic messages.
#'
#' @return The checked declaration without changing its numerical values or
#'   names.
#' @keywords internal
#' @noRd
.normalise_marker_weight_offset_request <- function(offset, name) {
  check_vector <- function(values, vector_name) {
    if (!is.numeric(values) || length(values) < 1L || any(!is.finite(values))) {
      cli::cli_abort("{.arg {vector_name}} must be a non-empty finite numeric vector.")
    }
    value_names <- names(values)
    if (!is.null(value_names)) {
      named <- !is.na(value_names) & nzchar(value_names)
      if (any(named) && !all(named)) {
        cli::cli_abort(c(
          x = "{.arg {vector_name}} must be wholly named or wholly unnamed.",
          i = "Name every marker entry, or remove every marker name and use fitted marker order."
        ))
      }
      if (all(named) && anyDuplicated(value_names)) {
        cli::cli_abort("{.arg {vector_name}} contains duplicated marker names.")
      }
    }
    values
  }

  if (is.null(offset)) return(NULL)
  if (is.numeric(offset) && is.null(dim(offset))) return(check_vector(offset, name))
  if (is.matrix(offset) || is.data.frame(offset)) {
    if (!all(vapply(offset, is.numeric, logical(1)))) {
      cli::cli_abort("{.arg {name}} must contain only numeric marker offsets.")
    }
    if (any(!is.finite(as.matrix(offset)))) {
      cli::cli_abort("{.arg {name}} must contain only finite marker offsets.")
    }
    column_names <- colnames(offset)
    if (!is.null(column_names)) {
      named_columns <- !is.na(column_names) & nzchar(column_names)
      if (any(named_columns) && !all(named_columns)) {
        cli::cli_abort("The marker columns of {.arg {name}} must be wholly named or wholly unnamed.")
      }
      if (all(named_columns) && anyDuplicated(column_names)) {
        cli::cli_abort("The marker columns of {.arg {name}} contain duplicated names.")
      }
    }
    row_names <- rownames(offset)
    if (is.null(row_names) || anyNA(row_names) || any(!nzchar(row_names)) || anyDuplicated(row_names)) {
      cli::cli_abort(c(
        x = "The rows of {.arg {name}} must be uniquely named.",
        i = "Name each row by its weighted association term."
      ))
    }
    unknown_terms <- setdiff(row_names, .weighted_assoc_term_keys())
    if (length(unknown_terms) > 0L) {
      cli::cli_abort("The rows of {.arg {name}} contain unknown weighted association terms: {.val {paste(unknown_terms, collapse = ', ')}}.")
    }
    return(offset)
  }
  if (is.list(offset)) {
    offset_names <- names(offset)
    if (is.null(offset_names) || anyNA(offset_names) || any(!nzchar(offset_names)) || anyDuplicated(offset_names)) {
      cli::cli_abort(c(
        x = "{.arg {name}} must be a completely named list.",
        i = "Name each vector by its weighted association term."
      ))
    }
    unknown_terms <- setdiff(offset_names, .weighted_assoc_term_keys())
    if (length(unknown_terms) > 0L) {
      cli::cli_abort("{.arg {name}} contains unknown weighted association terms: {.val {paste(unknown_terms, collapse = ', ')}}.")
    }
    return(stats::setNames(lapply(seq_along(offset), function(index) {
      check_vector(offset[[index]], paste0(name, "$", offset_names[[index]]))
    }), offset_names))
  }
  cli::cli_abort(c(
    x = "{.arg {name}} has an unsupported form.",
    i = "Use a numeric vector, or a named collection of numeric vectors for term-specific offsets."
  ))
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
    i = "A positive numeric concentration is also accepted."
  ))
}
