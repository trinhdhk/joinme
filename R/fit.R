#' Fit Joint Nested Mixed-effects (JoiNMe) model via cmdstanr or rstan
#'
#' @importFrom stats setNames
#' @importFrom utils modifyList
#'
#' @description
#' Fits the Joint Nested Mixed-effects (JoiNMe) Stan model using cmdstanr (default) or rstan.
#' The engine can be set globally via options(stan_preferred_engine = "cmdstanr"|"rstan") or via
#' control$engine.
#'
#' @details
#' Transformations are provided through a single `transforms` argument.
#' Each element (cv_total, cs_total, corr, vcov) is defined as a list specifying a type
#' and type-specific fields (see `build_standata_transforms()` for details). For
#' covariance-style associations, `corr` acts on the off-diagonal entries of the
#' subject-specific Cholesky-correlation factor `K`, while `vcov` acts on those
#' same off-diagonal `K` entries together with the subject-specific standard
#' deviations. If both `corr` and `vcov` are requested, `vcov` is kept and
#' `corr` is ignored with a warning.
#'
#' Distributional regression can be specified in two equivalent forms:
#' 1. A named list with RHS-only formulas, e.g. `list(sigma = ~ 1 + time)`.
#' 2. An unnamed list with explicit LHS parameter names, e.g. `list(sigma ~ 1 + time)`.
#'
#' The time axis is internally scaled for numerical stability; reported coefficients
#' are rescaled within Stan so fixed effects remain interpretable on the original scale.
#'
#' Marker weights (see `joinme_standata()`) are used to form marker-average summaries
#' for both current value (CV) and current slope (CS) association components. When
#' `shared_marker_weights = TRUE`, all weighted marker-based association terms share
#' one marker-weight structure. When `shared_marker_weights = FALSE`, each active
#' weighted marker-based association term (`cv_total`, `cs_total`, `cv_marker`,
#' `cs_marker`) gets its own marker-weight structure. When
#' `fixed_marker_weights = FALSE`, signed perturbations are estimated around the
#' supplied base weights. The effective marker intensities are used directly,
#' without additional normalisation, to scale the corresponding association
#' contribution.
#'
#' The returned `JoiNMeFit` object stores a compact association plotting bundle
#' containing only the posterior quantities needed to draw association curves
#' (`alpha_*`, marker weights, spline coefficients, and cached model-implied raw
#' support ranges). This keeps association plotting usable after serialization
#' without needing the full transient CmdStan CSV outputs.
#'
#' For CmdStanR fits, `joinme()` also eagerly imports the CSV-backed fit
#' contents in memory before returning. This mirrors the loading step used by
#' `cmdstanr::save_object()` so later `saveRDS()` calls do not rely on the
#' original CmdStan CSV files remaining on disk.
#'
#' `joinme()` can also fit the nested multivariate longitudinal model without
#' an event process. Leave both `formulaEvent` and `dataEvent` as `NULL`. The
#' longitudinal likelihood and every random-effect block are retained, whereas
#' the survival likelihood and longitudinal--survival association are removed.
#' Internally, a likelihood-neutral event scaffold is used only to reuse the
#' common design and Stan infrastructure. Its event parameters are not reported
#' as fitted scientific results.
#'
#' Formula-scoped subject weighting is supported in random-effect grouping terms via
#' `weighted(group, weights = <column>)`. For example:
#' - `(1 + time | weighted(id, weights = id_w))`
#' - `(0 + x1 + (1 + time | weighted(id, weights = id_w)) | weighted(marker, weights = marker_w))`
#' - `sigma ~ 1 + (1 | weighted(id, weights = id_w))`
#'
#' Weighted grouping affects latent random-effect priors for the referenced grouping
#' factor and applies subject-level likelihood weighting for the longitudinal and
#' survival contributions. Weight columns must be finite and strictly positive.
#'
#' Transformation for supported association terms
#' can be fed to `joinme()` in a named list. For example
#' `list(cv_total = list(type = "functional", expr = ~ log1p(x)),
#'      cs_total = list(type = "identity"),
#'      corr = list(type = "pwlin", knots = c(-2, 0, 2), direction = "increasing"),
#'      vcov = list(type = "identity"))`
#'
#' Transformation parameterisation:
#' - Omitted term, `NULL`, or `list(type = "identity")`:
#'   identity transform (default).
#' - `list(type = "functional", expr = ~ log1p(x))`:
#'   functional transform; `expr` is required and has no additional defaults.
#' - `list(type = "ispline", knots = c(-1, 0, 1), coeff = c(0, 0.3, 0.8, 1.1, 1.3), degree = 3)`:
#'   monotone I-spline; `knots` and `coeff` are required, `degree` defaults to `3`.
#' - `list(type = "ispline_penalised", x = seq(-2, 2, length.out = 50), y = exp(seq(-2, 2, length.out = 50)), n_knots = 6, degree = 3, lambda = 1)`:
#'   penalised monotone I-spline. This is intended for simulation but it works here too, off-labelly; defaults are
#'   `n_knots = 6`, `degree = 3`, `lambda = 1` when those values are omitted.
#' - `list(type = "ispline_penalised", x = seq(-2, 2, length.out = 50), n_knots = 6, degree = 3, lambda = 1)`:
#'   penalised monotone I-spline with Stan-estimated coefficients; if `knots`
#'   are omitted they are derived from `x` quantiles using `n_knots = 6` by default.
#' - `list(type = "ispline_expit", knots = c(0.05, 0.5, 0.95), coeff = c(0, 0.25, 0.8, 1.0, 1.1), degree = 3)`:
#'   monotone I-spline evaluated on `plogis(x)`; explicit `knots` are specified
#'   on the expit scale because the spline basis itself lives on that bounded
#'   `expit(x)` domain. This maybe useful for numerical stability
#'   but highly experimental and may be useful
#'   to put for knots at the each (near 0 and 1) if users is concerned
#'   about the curve at the tails.
#' - `list(type = "ispline_expit_penalised", x = seq(0.02, 0.98, length.out = 50), y = seq(0.02, 0.98, length.out = 50)^0.8, n_knots = 6, degree = 3, lambda = 1)`:
#'   penalised monotone I-spline on `plogis(x)`. This is intended for simulation but it works here too, off-labelly.
#' - `list(type = "ispline_expit_penalised", x = seq(0.02, 0.98, length.out = 50), n_knots = 6, degree = 3, lambda = 1)`:
#'   penalised monotone I-spline on `plogis(x)` with Stan-estimated coefficients.
#' - `list(type = "pwlin", knots = c(-2, -1, 0, 1, 2), direction = "increasing")`:
#'   ordered piecewise-linear association estimated jointly in Stan.
#'
#' For fitted piecewise-linear associations, the knots divide the raw
#' association feature into linear intervals. With $K$ knots, Stan estimates
#' a $K-1$ simplex. Its cumulative sums give ordered relative association
#' ordinates running from zero to one for an increasing curve, or from zero to
#' minus one for a decreasing curve. The signed association coefficient
#' estimates the total log-hazard span. The first ordinate is anchored at zero
#' because a free common ordinate would be confounded with the baseline hazard.
#'
#'
#' For monotone spline transforms:
#' - `type = "ispline"`: provide `knots` and `coeff` directly (plus optional
#'   `degree`); `x`/`y`/`lambda` are not used.
#' - `type = "ispline_penalised"` (alias: `"ispline_penalized"`): provide
#'   spline structure through `knots` (or `n_knots`) and smoothness penalty `lambda`.
#'   - If `y` is supplied, `JoiNMe` first fits the monotone spline to training
#'     pairs `(x, y)` in R and passes fixed coefficients to Stan. These plug-in
#'     coefficients use the same anchored convention as the Stan-estimated path:
#'     increasing splines run from `0` to `1`, while decreasing splines run from
#'     `1` to `0`.
#'   - If `y` is omitted, Stan estimates the monotone spline coefficients
#'     directly. Use `direction = "decreasing"` when the monotone transform
#'     should fall as the raw association feature increases.
#' - `x`: raw association-feature values used either to define training pairs or
#'   to help derive knot locations.
#' - `y`: optional target transformed values at those `x` points. Supplying `y`
#'   activates the legacy plug-in fit; omitting it activates Stan estimation.
#' - `lambda`: smoothness control (larger = smoother transform).
#' - `direction`: monotone orientation for penalised spline families. Accepted
#'   values are `"increasing"` and `"decreasing"`. When `y` is supplied, the
#'   plug-in fit infers the direction from the training pairs if you omit it.
#'   When `y` is omitted and Stan estimates the spline directly, the default
#'   remains `"increasing"`.
#' - `type = "ispline_expit"` / `"ispline_expit_penalised"`:
#'   same semantics as the I-spline variants above, except the spline basis is
#'   built on `plogis(x)`. This keeps the spline input on the bounded interval
#'   $(0, 1)$ and is often more numerically stable when the raw association
#'   feature spans a wide range. My experiments showed that this may help with
#'   a workaround for the Boundary knots but may cause considerable suppress at
#'   the tails if mis-specificied. Again, this may not affect the predictive value
#'   at much because the tail of the association term is often scarced. However,
#'   if you think the distribution of the latent association term is heavy-tailed,
#'   you may want to put more knots at the tails (e.g knots = c(0.001, 0.01, 0.5, 0.99, 0.999))
#'   to fight the impact of the `plogis()` transformation. In the future, a
#'   normalisation of the latent association term may be added to the model; may be
#'   dividing by sqrt of the second moment.
#'
#' Additional arguments are forwarded to `joinme_standata()` (e.g., `assoc`,
#' `basehaz`, `basehaz_degree`, `n_knots`, `time_var`, and `shrinkage`). Use
#' lme4-style `||` in `formulaLong` to request diagonal random-effect covariance.
#'
#' @param formulaLong Longitudinal formula defining fixed effects, id-level effects,
#'   and the marker block. The marker block may optionally include an inner
#'   `( ... | id )` term for marker-by-id random effects. When omitted,
#'   marker-by-id effects are disabled (Q_idm = 0). Grouping terms may use
#'   `||` at either level: outer `( ... || marker )` keeps marker-only and
#'   marker-by-id blocks independent, while inner `( ... || id )` keeps the
#'   marker-by-id covariance diagonal. Grouping terms may also use
#'   `weighted(group, weights = <column>)` to declare
#'   formula-scoped subject/group weights.
#' @param dataLong Long-format longitudinal data with columns for id, marker, time,
#'   outcome, and covariates referenced in `formulaLong`.
#' @param formulaEvent Optional survival formula for baseline covariates and
#'   the event model. Supply it together with `dataEvent`, or leave both `NULL`
#'   to fit only the nested longitudinal mixed model.
#'   Supported LHS forms are `survival::Surv(time, status)`,
#'   `survival::Surv(start, stop, status)`,
#'   `survival::Surv(time, status, type = "left")`, and
#'   `survival::Surv(time1, time2, type = "interval2")`.
#'   The legacy `type = "interval"` representation is not supported.
#' @param dataEvent Optional event-process data with either one row per id
#'   (`Surv(time, status)`) or multiple interval rows per id
#'   (`Surv(start, stop, status)`). Covariates in `formulaEvent` may vary by
#'   interval. Leave this and `formulaEvent` as `NULL` for longitudinal-only
#'   fitting.
#' @param formulaVCov Covariance regression formula for id-specific marker-by-id effects.
#'   If the marker block omits the inner `( ... | id )`, marker-by-id effects are
#'   absent and covariance-style associations (`corr`, `vcov`) are not allowed.
#'   When present, `corr` associations use the off-diagonal entries of the
#'   subject-specific Cholesky-correlation factor `K`, whereas `vcov`
#'   associations use those same off-diagonal `K` entries together with the
#'   subject-specific standard deviations. If both `corr` and `vcov` are
#'   requested, `vcov` is kept and `corr` is ignored with a warning. The
#'   default `~ 1` remains a supported intercept-only covariance regression.
#' @param formulaDist Optional list of formulas for distributional regression.
#'   Supported LHS parameters are:
#'   - `sigma` (scale parameter)
#'   - `nu`    (degrees of freedom for Student-t)
#'   - `phi`   (precision for Negative Binomial 2)
#'   - `alpha` (skewness parameter)
#'   - `kappa` (positive sample-size parameter for the Beta distribution, with
#'     shapes \eqn{\mu\kappa} and \eqn{(1-\mu)\kappa})
#'   - `tau` (quantile/asymmetry parameter in \eqn{(0,1)} for the skew double
#'     exponential distribution)
#'
#'   Three input styles are supported:
#'   1. Named list with RHS-only formulas, e.g. `list(sigma = ~ 1 + time)`.
#'   2. Unnamed list with explicit LHS, e.g. `list(sigma ~ 1 + time)`.
#'   3. Family-scoped LHS using brackets, e.g.
#'      `list(sigma[family=student_t] ~ 1 + time, sigma[family=gaussian] ~ 1 + x1)`.
#'
#'   Semantics of scoping:
#'   - Without `[family=...]`, the formula applies to all marker rows where that
#'     parameter is defined.
#'   - With `[family=...]`, the formula is applied only to rows of that family;
#'     other rows receive zero contribution from that scoped block.
#'   - Multiple scoped formulas for the same parameter are allowed and are
#'     estimated jointly with separate coefficients.
#'
#'   Random-effects terms with `|` are supported, but nested random-effects formulas
#'   are not. Grouping factors in random-effects terms may use
#'   `weighted(group, weights = <column>)`.
#' @param control Named list containing sampling configuration.
#'   - cmdstanr::model$sample() arguments (e.g., chains, parallel_chains,
#'     iter_warmup, iter_sampling, seed, refresh, adapt_delta, max_treedepth).
#'   - engine: "cmdstanr" or "rstan". Defaults to options(stan_preferred_engine).
#'   - threads_per_chain: integer; the threaded Stan program is always used.
#'     `threads_per_chain = 1` keeps execution serial while preserving the
#'     thread-capable kernel. For engine = "rstan", threading uses
#'     options(stan.thread = threads_per_chain).
#'   - grainsize: integer; reduce_sum grainsize for the threaded kernel.
#'     Defaults to the full subject count when `threads_per_chain = 1`, and to
#'     `max(1, min(n_cores, ceiling(n_id/(4*threads_per_chain*chains))))`
#'     otherwise.
#'   - force_recompile: logical; recompile the Stan model if needed.
#'   - quadrature_nodes: optional positive integer total node target for survival
#'     integration. Allowed values are exactly 7/15/31/41/51/61. Only the node
#'     count is passed to Stan; GK nodes/weights are fixed in the Stan code.
#'   - vcov_diag_link: "softplus" or "exp" for covariance regression diagonals.
#' @param draws Optional number of posterior draws used for summaries (not sampling).
#' @param families Marker-specific family specification (optional).
#'   Can be a character vector of family names aligned to marker order, or a
#'   list of `jm_family(...)` entries with per-marker links.
#'   A skew-Laplace marker may fix its quantile/asymmetry parameter through
#'   `jm_family("skew_laplace", tau = 0.8)`; otherwise `tau` is estimated.
#'   Supported named forward links are `identity`, `log`, `logit`, `probit`,
#'   and `exp`; `jm_family()` also accepts an invertible formula link or a
#'   directly specified inverse-link formula. In formula syntax,
#'   `inv_Phi`/`qnorm`/`probit` denote the standard normal quantile and
#'   `Phi`/`pnorm` denote the standard normal CDF. Thus a probit forward link
#'   is inverted to `Phi` before fitting.
#' @param transforms Transformation specifications for association terms.
#'   Prefer declaring them with `joinme_tf(...)`; raw named lists remain
#'   supported. Fit-time functional transforms may also request a free affine
#'   shift through `intercept = TRUE` and/or `slope = TRUE`, for example
#'   `joinme_tf(cv_total = ~ expit(x, intercept = TRUE, slope = TRUE))`, which
#'   is fitted as `expit(iota_1 + iota_2 * x)`. See details.
#'
#' @param priors Prior declaration. Prefer `joinme_priors(...)`; raw named lists
#'   with components `beta`, `alpha`, `iota`, and `lkj` remain supported.
#' @param fixed_marker_weights Logical; if TRUE, marker weights are fixed at the
#'   supplied base values. If FALSE, marker-weight perturbations are estimated
#'   using the family selected by `shrinkage` (0 = Student-t(6), 1 = Laplace,
#'   2 = Normal).
#' @param shared_marker_weights Logical; if TRUE, all weighted marker-based
#'   association terms share one marker-weight structure. If FALSE, each active
#'   weighted marker-based association term gets its own marker-weight structure.
#'   When `marker_weights` is a named list, use names `cv_total`, `cs_total`,
#'   `cv_marker`, and `cs_marker`.
#' @param basehaz An object of class `joinme_basehaz` created by `joinme_basehaz()`.
#' This controls the baseline hazard parameterisation and spline basis. See `?joinme_basehaz` for details.
#' @param fit logical; if TRUE, the model is fitted and a `JoiNMeFit` object is returned. If FALSE, only the Stan data list is returned.
#' @param seed Optional random seed for reproducibility.
#' @param ... Additional args passed to joinme_standata().
#'
#' @examples
#' \dontrun{
#' longitudinal_fit <- joinme(
#'   formulaLong = y ~ time + (1 + time | id) +
#'     (1 + time + (1 + time | id) | marker),
#'   dataLong = dataLong,
#'   families = rep("gaussian", 3)
#' )
#'
#' formulaDist <- list(
#'   sigma[family=student_t] ~ 1 + time + (1 | id),
#'   sigma[family=gaussian] ~ 1 + x1,
#'   nu[family=student_t] ~ 1,
#'   alpha[family=skew_normal] ~ 1 + x1,
#'   phi[family=negbin2] ~ 1,
#'   kappa[family=beta] ~ 1,
#'   tau[family=skew_double_exponential] ~ 1
#' )
#' }
#'
#' @export
# File overview:
# - Validate inputs and build Stan data.
# - Resolve threading/engine settings and select the Stan program.
# - Fit with cmdstanr/rstan and wrap results in a JoiNMeFit object.

joinme <- function(
  formulaLong,
  dataLong,
  formulaEvent = NULL,
  dataEvent = NULL,
  formulaVCov = ~1,
  formulaDist = NULL,
  control = list(),
  draws = NULL,
  families = NULL,
  transforms = NULL,
  priors = joinme_priors(),
  fixed_marker_weights = FALSE,
  shared_marker_weights = TRUE,
  basehaz = joinme_basehaz(),
  fit = TRUE,
  seed = NULL,
  ...
) {
  # Step 1: identify whether the user supplied a complete event process before
  # constructing any internal scaffold. A formula without event data, or event
  # data without its formula, has no unambiguous statistical interpretation and
  # is therefore rejected rather than silently treated as longitudinal-only.
  has_survival_process <- !is.null(formulaEvent) || !is.null(dataEvent)
  if (xor(is.null(formulaEvent), is.null(dataEvent))) {
    cli::cli_abort(c(
      x = "{.arg formulaEvent} and {.arg dataEvent} must be supplied together.",
      i = "Leave both NULL to fit only the nested longitudinal mixed model."
    ))
  }

  # Step 2: resolve association terms before calling `joinme_standata()`,
  # whose historical default is `cv_mean`. That default remains appropriate
  # for a joint longitudinal--survival model, but it must not accidentally
  # introduce a prior-only association parameter in a longitudinal-only fit.
  dot_arguments <- list(...)
  association_terms <- dot_arguments$assoc
  if (is.null(association_terms)) {
    association_terms <- if (has_survival_process) {
      "cv_mean"
    } else {
      character(0)
    }
  }
  if (!has_survival_process && length(association_terms) > 0L) {
    cli::cli_abort(c(
      x = "Association terms require a survival process.",
      i = "Remove {.arg assoc}, or supply both {.arg formulaEvent} and {.arg dataEvent}."
    ))
  }
  dot_arguments$assoc <- association_terms

  # Step 3: construct one administrative, event-free row per subject when the
  # analysis contains no survival outcome. Positive follow-up keeps the common
  # time-design code well-defined; `include_survival = FALSE` below then sets
  # all event likelihood contributions to exactly zero before sampling.
  if (!has_survival_process) {
    dataEvent <- .longitudinal_only_event_scaffold(
      data_long = dataLong,
      id_variable = dot_arguments$id_var %||% "id",
      time_variable = dot_arguments$time_var %||% "time"
    )
    formulaEvent <- survival::Surv(
      .joinme_follow_up,
      .joinme_event
    ) ~ 1
  }

  if (!is.list(control)) {
    cli::cli_abort(c(
      x = "{.arg control} must be a named list.",
      i = "Provide list(threads_per_chain=..., grainsize=..., adapt_delta=..., max_treedepth=...)."
    ))
  }
  transforms <- unclass(.normalise_joinme_tf_input(
    transforms,
    validate = FALSE
  ))
  priors <- unclass(.joinme_priors_(priors, validate = TRUE))

  if (length(control) > 0 && is.null(names(control))) {
    cli::cli_abort(c(
      x = "{.arg control} must be a named list.",
      i = "Provide list(threads_per_chain=..., grainsize=..., adapt_delta=..., max_treedepth=...)."
    ))
  }

  # Workflow: build standata -> resolve threading -> choose engine -> fit -> wrap
  # - sd: prepared Stan data list with all dimensions and transforms
  vcov_diag_link <- control$vcov_diag_link %||% "softplus"
  quadrature_nodes <- control$quadrature_nodes %||% NULL
  formulaVCov <- .get_vcov_formula(
    formulaVCov = formulaVCov,
    default = ~1,
    context = "joinme()"
  )
  assertthat::assert_that(
    inherits(basehaz, "joinme_basehaz"),
    msg = "{.arg basehaz} must be a {.cls joinme_basehaz} object, created by {.fn joinme_basehaz()} or {.fn jm_basehaz()}."
  )

  arg_list <- dot_arguments
  arg_list <- arg_list[names(arg_list) %in% names(formals(joinme_standata))]
  arg_list <- modifyList(
    arg_list,
    list(
      formulaLong = formulaLong,
      dataLong = dataLong,
      formulaEvent = formulaEvent,
      dataEvent = dataEvent,
      include_survival = has_survival_process,
      formulaVCov = formulaVCov,
      formulaDist = formulaDist,
      families = families,
      transforms = transforms,
      beta_prior = priors$beta,
      alpha_prior = priors$alpha,
      iota_prior = priors$iota,
      lkj_prior = priors$lkj,
      fixed_marker_weights = fixed_marker_weights,
      shared_marker_weights = shared_marker_weights,
      quadrature_nodes = quadrature_nodes,
      vcov_diag_link = vcov_diag_link,
      basehaz = basehaz$type,
      basehaz_n_knots = basehaz$n_knots,
      basehaz_knots = basehaz$knots,
      basehaz_degree = basehaz$degree,
      basehaz_formula = basehaz$formula,
      seed = seed
    )
  )
  sd <- do.call(joinme_standata, arg_list)

  # Threading: honor explicit control overrides, fall back to mc.cores or 1
  threads_per_chain <- control$threads_per_chain %||%
    control$threads %||%
    control$mc.cores %||%
    1L
  if (!is.numeric(threads_per_chain) || length(threads_per_chain) != 1) {
    cli::cli_abort(c(
      x = "{.arg control$threads_per_chain} must be a single numeric value.",
      i = "Example: control = list(threads_per_chain = 2)."
    ))
  }
  threads_per_chain <- as.integer(threads_per_chain)
  if (threads_per_chain < 1) {
    cli::cli_abort(c(
      x = "{.arg control$threads_per_chain} must be >= 1.",
      i = "Use 1 to disable threading."
    ))
  }
  # Core cap: never exceed available physical cores or subject count
  n_cores <- parallel::detectCores(logical = FALSE) %||% 1L
  max_threads <- min(sd$n_id %||% 1L, n_cores)
  if (threads_per_chain > max_threads) {
    cli::cli_warn(c(
      x = "Requested {threads_per_chain} threads exceeds max {max_threads}.",
      i = "Capping threads_per_chain to {max_threads}."
    ))
    threads_per_chain <- max_threads
  }

  # reduce_sum grainsize: default depends on id count, chains, and threads
  grainsize <- control$grainsize
  if (is.null(grainsize)) {
    n_id <- sd$n_id %||% 1L
    if (threads_per_chain <= 1L) {
      grainsize <- as.integer(n_id)
    } else {
      n_chains <- control$parallel_chains %||% control$chains %||% 4L
      denom <- 4L * as.integer(threads_per_chain) * as.integer(n_chains)
      denom <- max(1L, denom)
      grainsize <- max(1L, as.integer(ceiling(n_id / denom)))
      grainsize <- min(as.integer(n_cores), grainsize)
    }
  }
  if (!is.numeric(grainsize) || length(grainsize) != 1) {
    cli::cli_abort(c(
      x = "{.arg control$grainsize} must be a single numeric value.",
      i = "Example: control = list(grainsize = 5)."
    ))
  }
  grainsize <- as.integer(grainsize)
  if (grainsize < 1) {
    cli::cli_abort(c(
      x = "{.arg control$grainsize} must be >= 1.",
      i = "Use a positive integer for reduce_sum grainsize."
    ))
  }

  stan_file <- .get_stan_file(
    program = .stan_fit_program(sd),
    threaded = TRUE
  )

  engine <- .get_stan_engine(control$engine)
  if (engine == "rstan") {
    old_stan_thread <- getOption("stan.thread")
    options(stan.thread = threads_per_chain)
    on.exit(options(stan.thread = old_stan_thread), add = TRUE)
  }

  cpp_opts <- list(stan_threads = TRUE)
  if (engine == "cmdstanr") {
    mod <- .get_cmdstan_model(
      stan_file,
      cpp_options = cpp_opts,
      force_recompile = control$force_recompile %||% FALSE
    )
  } else {
    mod <- .get_rstan_model(
      stan_file
    )
  }
  if (engine == "cmdstanr" && !.cmdstan_threads_enabled(mod)) {
    cli::cli_abort(c(
      x = "The CmdStan model is not compiled with {.code stan_threads = TRUE}.",
      i = "Retry with {.code control = list(force_recompile = TRUE)} or call {.fn precompile_cmdstanr_models}."
    ))
  }

  # Remove non-Stan fields only
  sd_stan <- sd
  sd_stan$x_cols <- NULL
  sd_stan$w_cols <- NULL
  sd_stan$zid_cols <- NULL
  sd_stan$zmk_cols <- NULL
  sd_stan$zidm_cols <- NULL
  sd_stan$time_var <- NULL
  sd_stan$marker_levels <- NULL
  sd_stan$family_codes <- NULL
  sd_stan$family_names <- NULL
  sd_stan$tf_compositions <- NULL
  sd_stan$basehaz <- NULL
  sd_stan$n_knots <- NULL
  sd_stan$basehaz_n_knots <- NULL
  sd_stan$basehaz_knots <- NULL
  sd_stan$basehaz_degree <- NULL
  sd_stan$basehaz_formula <- NULL
  sd_stan$basehaz_col_means <- NULL
  sd_stan$basehaz_cols <- NULL
  sd_stan$Bs_obj <- NULL
  sd_stan$dist_cols <- NULL
  sd_stan$dist_re_terms <- NULL
  sd_stan$dist_formulas <- NULL

  if (is.null(sd_stan$id) || is.null(sd_stan$n_id)) {
    cli::cli_abort(c(
      x = "Threaded Stan execution requires {.arg id} and {.arg n_id} in Stan data.",
      i = "Check the standata builder output."
    ))
  }
  id_vec <- sd_stan$id
  idx <- split(seq_along(id_vec), id_vec)
  id_start <- vapply(idx, min, integer(1))
  id_end <- vapply(idx, max, integer(1))
  if (length(id_start) != sd_stan$n_id) {
    cli::cli_abort(c(
      x = "Threaded Stan execution requires contiguous ids from 1..n_id.",
      i = "Check the id mapping in standata."
    ))
  }
  sd_stan$id_start <- as.integer(id_start)
  sd_stan$id_end <- as.integer(id_end)
  sd_stan$grainsize <- grainsize

  # Ensure time index arrays are preserved for cmdstanr JSON (avoid auto-unbox)
  sd_stan <- .coerce_rstan_time_indices(sd_stan)
  sd_stan <- .coerce_rstan_mixture_data(sd_stan)

  # Coerce arrays/vectors consistently for both cmdstanr and rstan
  sd_stan <- .coerce_rstan_dist_arrays(sd_stan)
  sd_stan <- .coerce_rstan_vectors(
    sd_stan,
    c(
      "beta_scale",
      "const_data_cv",
      "const_data_cs",
      "const_data_corr",
      "const_data_vcov",
      "const_data_cv_mean",
      "const_data_cv_marker",
      "const_data_cs_mean",
      "const_data_cs_marker"
    )
  )

  has_nonstan_metadata <- function(x) {
    if (is.null(x)) {
      return(FALSE)
    }
    if (
      is.character(x) ||
        is.factor(x) ||
        is.language(x) ||
        inherits(x, c("formula", "call"))
    ) {
      return(TRUE)
    }
    if (is.list(x) && !is.data.frame(x)) {
      return(any(vapply(x, has_nonstan_metadata, logical(1))))
    }
    FALSE
  }
  keep_idx <- !vapply(sd_stan, has_nonstan_metadata, logical(1))
  sd_stan <- sd_stan[keep_idx]

  allowed_data <- if (engine == "cmdstanr") {
    vars <- tryCatch(mod$variables(), error = function(e) NULL)
    if (!is.null(vars$data)) {
      names(vars$data) %||% character(0)
    } else {
      character(0)
    }
  } else {
    .stan_data_names(stan_file)
  }
  required_data <- .stan_data_names(stan_file)
  if (length(allowed_data) == 0) {
    allowed_data <- required_data
  }
  if (length(required_data) > 0) {
    missing <- setdiff(required_data, names(sd_stan))
    if (length(missing) > 0) {
      cli::cli_abort(c(
        x = "Stan data missing required fields: {paste(missing, collapse = ', ')}.",
        i = "Check joinme_standata() and Stan data block consistency."
      ))
    }
  }
  if (length(allowed_data) > 0) {
    keep_data <- union(allowed_data, required_data)
    sd_stan <- sd_stan[names(sd_stan) %in% keep_data]
  }

  defaults <- list(
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    seed = 1,
    init = 1,
    refresh = 100
  )
  sample_control <- control[setdiff(
    names(control),
    c(
      "engine",
      "threads_per_chain",
      "threads",
      "mc.cores",
      "grainsize",
      "quadrature_nodes",
      "force_recompile"
    )
  )]
  args <- modifyList(defaults, sample_control)
  chains_val <- args$parallel_chains %||%
    args$chains %||%
    defaults$parallel_chains
  args$parallel_chains <- chains_val
  args$chains <- chains_val
  args$threads_per_chain <- threads_per_chain
  args$data <- sd_stan
  args <- args[!vapply(args, is.null, logical(1))]

  if (engine == "cmdstanr") {
    allowed <- names(formals(mod$sample))
    args <- args[names(args) %in% allowed]
  } else {
    control_list <- list()
    if (!is.null(args$adapt_delta)) {
      control_list$adapt_delta <- args$adapt_delta
    }
    if (!is.null(args$max_treedepth)) {
      control_list$max_treedepth <- args$max_treedepth
    }
    iter_warmup <- args$iter_warmup %||% defaults$iter_warmup
    iter_sampling <- args$iter_sampling %||% defaults$iter_sampling
    iter_total <- iter_warmup + iter_sampling
    rstan_args <- list(
      object = mod,
      data = sd_stan,
      chains = args$chains %||% 4,
      iter = iter_total,
      warmup = iter_warmup,
      seed = args$seed %||% defaults$seed,
      refresh = args$refresh %||% defaults$refresh,
      open_progress = FALSE,
      cores = min(
        args$chains %||% 1,
        parallel::detectCores(logical = FALSE) %||% 1
      )
    )
    if (length(control_list) > 0) rstan_args$control <- control_list
  }

  call <- match.call()
  sd_recipe <- list(
    formulaLong = formulaLong,
      formulaEvent = formulaEvent,
      formulaVCov = formulaVCov,
      parent_call = call,
      dataLong = dataLong,
      dataEvent = dataEvent,
      transforms = transforms,
      draws = draws,
      threads_per_chain = threads_per_chain,
      # Retain the complete model description for posterior reporting.  The
      # sampling arguments already contain `sd_stan`, which is the numeric
      # subset accepted by Stan.  Keeping the unabridged `sd` here preserves
      # design-column labels, marker labels, family labels, and baseline-hazard
      # metadata needed to translate fitted parameters into statistical terms.
      stan_data = sd,
      stan_mod = mod,
      stan_engine = engine,
      stan_args = if (engine == "cmdstanr") args else rstan_args
  )
  jm_sd <- JoiNMeStanData$new(sd_recipe)
  # If no fitting is requested, return the prepared Stan data list, program sample args
  if (!fit) {
    return(jm_sd)
  }

  jm_sd$sample()
}

#' Construct a likelihood-neutral event scaffold
#'
#' @param data_long Long-format longitudinal data.
#' @param id_variable Name of the subject identifier.
#' @param time_variable Name of the longitudinal time variable.
#'
#' @return One data-frame row per subject, with positive administrative
#'   follow-up and a zero event indicator.
#' @keywords internal
#' @noRd
.longitudinal_only_event_scaffold <- function(
  data_long,
  id_variable,
  time_variable
) {
  if (!is.data.frame(data_long)) {
    cli::cli_abort("{.arg dataLong} must be a data frame.")
  }
  missing_variables <- setdiff(
    c(id_variable, time_variable),
    names(data_long)
  ) # longitudinal columns required to identify subjects and follow-up
  if (length(missing_variables) > 0L) {
    cli::cli_abort(
      "Longitudinal-only fitting requires columns: {paste(missing_variables, collapse = ', ')}."
    )
  }

  valid_rows <- !is.na(data_long[[id_variable]]) &
    is.finite(
      as.numeric(data_long[[time_variable]])
    ) # observations carrying both an identifier and a finite measurement time
  observed_data <- data_long[
    valid_rows,
    ,
    drop = FALSE
  ] # eligible longitudinal rows from which subject records can be copied
  if (nrow(observed_data) == 0L) {
    cli::cli_abort("No finite longitudinal observation time is available.")
  }

  ordered_rows <- order(
    as.character(observed_data[[id_variable]]),
    as.numeric(observed_data[[time_variable]])
  ) # subject-major ordering with the latest observation last
  observed_data <- observed_data[
    ordered_rows,
    ,
    drop = FALSE
  ]
  scaffold <- observed_data[
    !duplicated(
      as.character(observed_data[[id_variable]]),
      fromLast = TRUE
    ),
    ,
    drop = FALSE
  ] # one covariate-complete row copied from the last observation of each subject
  overall_follow_up <- max(
    as.numeric(observed_data[[time_variable]]),
    na.rm = TRUE
  ) # common positive administrative time needed only by the design builder
  if (!is.finite(overall_follow_up) || overall_follow_up <= 0) {
    overall_follow_up <- 1
  }
  scaffold$.joinme_follow_up <-
    overall_follow_up # neutral administrative endpoint for every subject
  scaffold$.joinme_event <- 0L # zero event indicator removing endpoint hazards
  rownames(scaffold) <- NULL
  scaffold
}

#' Determine whether a fitted model contains an observed event process
#'
#' @param object A fitted JoiNMe object carrying configuration or Stan data.
#'
#' @return A single logical value. Objects created before longitudinal-only
#'   support are conservatively treated as joint models.
#' @keywords internal
#' @noRd
.fit_includes_survival <- function(object) {
  if (inherits(object, "JoiNMeMixFit")) {
    mixture_specification <- object$mixture %||%
      object$config$mixture # checked latent-class model description
    if (!is.null(mixture_specification$include_survival)) {
      return(isTRUE(mixture_specification$include_survival))
    }
  }
  isTRUE(
    object$config$include_survival %||%
      (as.integer(object$stan_data$include_survival %||% 1L) == 1L)
  )
}

#' Require an observed event process for a survival-specific operation
#'
#' @param object A fitted JoiNMe object.
#' @param operation Plain-language name of the requested operation.
#'
#' @return The input object, invisibly, when the requirement is met.
#' @keywords internal
#' @noRd
.require_fitted_survival_process <- function(object, operation) {
  if (!.fit_includes_survival(object)) {
    cli::cli_abort(c(
      x = "{operation} is unavailable for a longitudinal-only fit.",
      i = "Fit with both {.arg formulaEvent} and {.arg dataEvent} to model an event process."
    ))
  }
  invisible(object)
}
