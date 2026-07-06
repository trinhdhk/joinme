#' Fit JoiNMe model via cmdstanr or rstan
#'
#' @importFrom stats setNames
#' @importFrom utils modifyList
#'
#' @description
#' Fits the JoiNMe Stan model using cmdstanr (default) or rstan. The engine can be
#' set globally via options(stan_preferred_engine = "cmdstanr"|"rstan") or via
#' control$engine.
#' following the notebook, and supports association transformations and time-internal
#' scaling metadata.
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
#' For CmdStanR fits, `joinme()` also eagerly materializes the CSV-backed fit
#' contents in memory before returning. This mirrors the loading step used by
#' `cmdstanr::save_object()` so later `saveRDS()` calls do not rely on the
#' original CmdStan CSV files remaining on disk.
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
#' Association coefficients are denoted with the `alpha_` prefix to match joint-model
#' conventions and to avoid confusion with the linear predictor eta used throughout
#' the longitudinal and survival submodels. Total, marker, and mean associations use
#' a positive non-centred parameterisation: `alpha = z_alpha * sd_alpha`, which
#' stabilizes sampling while preserving the intended sign structure via marker weights.
#'
#' The typical workflow is:
#' 1. Prepare `dataLong` and `dataEvent` with aligned ids and time scales.
#' 2. Specify `formulaLong` and `formulaEvent` to define the longitudinal and survival submodels.
#' 3. Choose `assoc` components and optional `transforms` for association terms.
#' 4. Fit with `control` to manage sampling, threading, and reproducibility.
#' 5. Validate convergence (R-hat, ESS, divergences), then predict.
#'
#' Example:
#' `list(cv_total = list(type = "functional", expr = ~ log1p(x)),
#'      cs_total = list(type = "identity"),
#'      corr = list(type = "pwlin", x = c(-2, 0, 2), y = c(0.2, 1, 0.2)),
#'      vcov = list(type = "identity"))`
#'
#' Transformation parameterisation cheat sheet:
#' - Omitted term, `NULL`, or `list(type = "identity")`:
#'   identity transform (default).
#' - `list(type = "functional", expr = ~ log1p(x))`:
#'   functional transform; `expr` is required and has no additional defaults.
#' - `list(type = "ispline", knots = c(-1, 0, 1), coeff = c(0, 0.3, 0.8, 1.1, 1.3), degree = 3)`:
#'   monotone I-spline; `knots` and `coeff` are required, `degree` defaults to `3`.
#' - `list(type = "ispline_penalised", x = seq(-2, 2, length.out = 50), y = exp(seq(-2, 2, length.out = 50)), n_knots = 6, degree = 3, lambda = 1)`:
#'   penalised monotone I-spline in legacy plug-in mode; defaults are
#'   `n_knots = 6`, `degree = 3`, `lambda = 1` when those values are omitted.
#' - `list(type = "ispline_penalised", x = seq(-2, 2, length.out = 50), n_knots = 6, degree = 3, lambda = 1)`:
#'   penalised monotone I-spline with Stan-estimated coefficients; if `knots`
#'   are omitted they are derived from `x` quantiles using `n_knots = 6` by default.
#' - `list(type = "ispline_expit", knots = c(0.05, 0.5, 0.95), coeff = c(0, 0.25, 0.8, 1.0, 1.1), degree = 3)`:
#'   monotone I-spline evaluated on `plogis(x)`; explicit `knots` are specified
#'   on the expit scale because the spline basis itself lives on that bounded
#'   `expit(x)` domain.
#' - `list(type = "ispline_expit_penalised", x = seq(0.02, 0.98, length.out = 50), y = seq(0.02, 0.98, length.out = 50)^0.8, n_knots = 6, degree = 3, lambda = 1)`:
#'   penalised monotone I-spline on `plogis(x)` in legacy plug-in mode.
#' - `list(type = "ispline_expit_penalised", x = seq(0.02, 0.98, length.out = 50), n_knots = 6, degree = 3, lambda = 1)`:
#'   penalised monotone I-spline on `plogis(x)` with Stan-estimated coefficients.
#' - `list(type = "pwlin", x = c(-2, -1, 0, 1, 2), y = c(0.2, 0.5, 1, 0.5, 0.2))`:
#'   piecewise-linear transform; both `x` and `y` are required and there are no
#'   additional defaults.
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
#'   remains `"increasing"` for backward compatibility.
#' - `type = "ispline_expit"` / `"ispline_expit_penalised"`:
#'   same semantics as the I-spline variants above, except the spline basis is
#'   built on `plogis(x)`. This keeps the spline input on the bounded interval
#'   $(0, 1)$ and is often more numerically stable when the raw association
#'   feature spans a wide range. Explicit `knots` and training `x` values must
#'   be specified on that expit scale.
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
#' @param formulaEvent Survival formula for baseline covariates and event model.
#' @param dataEvent One row per id event data with event time, event indicator, and
#'   covariates referenced in `formulaEvent`.
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
#'   - `sigma`
#'   - `nu`
#'   - `phi`
#'   - `alpha` (aliases: `alpha_skew`, `skew`)
#'   - `phi_beta`
#'   - `tau_sde`
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
#'   - tau_sde_fixed: optional fixed tau in (0,1) for skew-double-exponential.
#' @param draws Optional number of posterior draws used for summaries (not sampling).
#' @param families Marker-specific family specification (optional).
#'   Can be a character vector of family names aligned to marker order, or a
#'   list of `jm_family(...)` entries with per-marker links.
#'   Supported links/inverse-links: `identity`, `log`, `logit`, `probit`, `exp`.
#' @param transforms Transformation specifications for association terms.
#'   Prefer declaring them with `joinme_tf(...)`; raw named lists remain
#'   supported. Fit-time functional transforms may also request a free affine
#'   shift through `intercept = TRUE` and/or `slope = TRUE`, for example
#'   `joinme_tf(cv_total = ~ expit(x, intercept = TRUE, slope = TRUE))`, which
#'   is fitted as `expit(iota_1 + iota_2 * x)`. This fit-only affine shift is
#'   not mirrored by `simulate_joinme()`.
#' @param priors Prior declaration. Prefer `joinme_priors(...)`; raw named lists
#'   with components `beta`, `alpha`, `iota`, and `lkj` remain supported.
#' @param fixed_marker_weights Logical; if TRUE, marker weights are fixed at the
#'   supplied base values. If FALSE, marker-weight perturbations are estimated.
#' @param shared_marker_weights Logical; if TRUE, all weighted marker-based
#'   association terms share one marker-weight structure. If FALSE, each active
#'   weighted marker-based association term gets its own marker-weight structure.
#'   When `marker_weights` is a named list, use names `cv_total`, `cs_total`,
#'   `cv_marker`, and `cs_marker`.
#' @param ... Additional args passed to joinme_standata().
#'
#' @examples
#' \dontrun{
#' formulaDist <- list(
#'   sigma[family=student_t] ~ 1 + time + (1 | id),
#'   sigma[family=gaussian] ~ 1 + x1,
#'   nu[family=student_t] ~ 1,
#'   alpha[family=skew_normal] ~ 1 + x1,
#'   phi[family=negbin2] ~ 1,
#'   phi_beta[family=beta] ~ 1,
#'   tau_sde[family=skew_double_exponential] ~ 1
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
  formulaEvent,
  dataEvent,
  formulaVCov = ~1,
  formulaDist = NULL,
  control = list(),
  draws = NULL,
  families = NULL,
  transforms = NULL,
  priors = joinme_priors(),
  fixed_marker_weights = FALSE,
  shared_marker_weights = TRUE,
  ...
) {
  if (!is.list(control)) {
    cli::cli_abort(c(
      x = "{.arg control} must be a named list.",
      i = "Provide list(threads_per_chain=..., grainsize=..., adapt_delta=..., max_treedepth=...)."
    ))
  }
  transforms <- unclass(.normalise_joinme_tf_input(transforms, validate = FALSE))
  priors <- unclass(.normalise_joinme_priors_input(priors, validate = TRUE))

  if (length(control) > 0 && is.null(names(control))) {
    cli::cli_abort(c(
      x = "{.arg control} must be a named list.",
      i = "Provide list(threads_per_chain=..., grainsize=..., adapt_delta=..., max_treedepth=...)."
    ))
  }

  # Workflow: build standata -> resolve threading -> choose engine -> fit -> wrap
  # - sd: prepared Stan data list with all dimensions and transforms
  vcov_diag_link <- control$vcov_diag_link %||% "softplus"
  tau_sde_fixed <- control$tau_sde_fixed %||% NULL
  quadrature_nodes <- control$quadrature_nodes %||% NULL
  formulaVCov <- .resolve_vcov_formula(
    formulaVCov = formulaVCov,
    default = ~ 1,
    context = "joinme()"
  )
  sd <- joinme_standata(
    formulaLong = formulaLong,
    dataLong = dataLong,
    formulaEvent = formulaEvent,
    dataEvent = dataEvent,
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
    tau_sde_fixed = tau_sde_fixed,
    ...
  )

  # Threading: honor explicit control overrides, fall back to mc.cores or 1
  threads_per_chain <- control$threads_per_chain %||% control$threads %||% control$mc.cores %||% 1L
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
    program = "JoiNMe_fit",
    threaded = TRUE
  )

  engine <- .resolve_stan_engine(control$engine)
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
  sd_stan$basehaz_degree <- NULL
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

  # Coerce arrays/vectors consistently for both cmdstanr and rstan
  sd_stan <- .coerce_rstan_dist_arrays(sd_stan)
  sd_stan <- .coerce_rstan_vectors(sd_stan, c(
    "beta_scale",
    "const_data_cv",
    "const_data_cs",
    "const_data_corr",
    "const_data_vcov",
    "const_data_cv_mean",
    "const_data_cv_marker",
    "const_data_cs_mean",
    "const_data_cs_marker"
  ))

  has_nonstan_metadata <- function(x) {
    if (is.null(x)) return(FALSE)
    if (is.character(x) || is.factor(x) || is.language(x) || inherits(x, c("formula", "call"))) return(TRUE)
    if (is.list(x) && !is.data.frame(x)) {
      return(any(vapply(x, has_nonstan_metadata, logical(1))))
    }
    FALSE
  }
  keep_idx <- !vapply(sd_stan, has_nonstan_metadata, logical(1))
  sd_stan <- sd_stan[keep_idx]

  allowed_data <- if (engine == "cmdstanr") {
    vars <- tryCatch(mod$variables(), error = function(e) NULL)
    if (!is.null(vars$data)) names(vars$data) %||% character(0) else character(0)
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
  sample_control <- control[setdiff(names(control), c(
    "engine",
    "threads_per_chain",
    "threads",
    "mc.cores",
    "grainsize",
    "quadrature_nodes",
    "force_recompile"
  ))]
  args <- modifyList(defaults, sample_control)
  chains_val <- args$parallel_chains %||% args$chains %||% defaults$parallel_chains
  args$parallel_chains <- chains_val
  args$chains <- chains_val
  args$threads_per_chain <- threads_per_chain
  args$data <- sd_stan
  args <- args[!vapply(args, is.null, logical(1))]

  if (engine == "cmdstanr") {
    allowed <- names(formals(mod$sample))
    args <- args[names(args) %in% allowed]
    fit <- do.call(mod$sample, args)
    fit <- .materialize_cmdstanr_fit(fit)
  } else {
    control_list <- list()
    if (!is.null(args$adapt_delta)) control_list$adapt_delta <- args$adapt_delta
    if (!is.null(args$max_treedepth)) control_list$max_treedepth <- args$max_treedepth
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
      cores = min(args$chains %||% 1, parallel::detectCores(logical = FALSE) %||% 1)
    )
    if (length(control_list) > 0) rstan_args$control <- control_list
    fit <- do.call(rstan::sampling, rstan_args)
  }

  cfg <- list(
    family_long = sd$family_long,
    family_link = sd$link_long,
    family_link_names = sd$link_names,
    family_inv_link_n_ops = sd$inv_link_n_ops,
    family_inv_link_ops = sd$inv_link_ops,
    family_inv_link_n_const = sd$inv_link_n_const,
    family_inv_link_const = sd$inv_link_const,
    assoc = c(
      cv_total = sd$assoc_cv_total,
      cv_mean = sd$assoc_cv_mean,
      cv_marker = sd$assoc_cv_marker,
      cs_total = sd$assoc_cs_total,
      cs_mean = sd$assoc_cs_mean,
      cs_marker = sd$assoc_cs_marker,
      corr = sd$assoc_corr,
      vcov = sd$assoc_vcov
    ),
    transforms = list(
      tf_mode_cv_tot = sd$tf_mode_cv_tot,
      tf_mode_cs_tot = sd$tf_mode_cs_tot,
      tf_mode_cv_mean = sd$tf_mode_cv_mean,
      tf_mode_cs_mean = sd$tf_mode_cs_mean,
      tf_mode_cv_marker = sd$tf_mode_cv_marker,
      tf_mode_cs_marker = sd$tf_mode_cs_marker,
      tf_mode_corr = sd$tf_mode_corr,
      tf_mode_vcov = sd$tf_mode_vcov
    ),
    transforms_spec = transforms,
    dist = list(
      dist_cols = sd$dist_cols,
      dist_re_terms = sd$dist_re_terms,
      dist_formulas = sd$dist_formulas
    ),
    indep = c(
      id = sd$indep_id_re,
      marker = sd$indep_marker_re,
      idmarker_cov = sd$indep_idmarker_cov
    ),
    allow_marker_crosscorr = sd$allow_marker_crosscorr,
    shrinkage = sd$shrinkage,
    dims = c(n_id = sd$n_id, N = sd$N, D = sd$D, P = sd$P, R_id = sd$R_id, R_mk = sd$R_mk, Q_idm = sd$Q_idm),
    draws_default = draws,
    threads_per_chain = threads_per_chain,
    tmax_internal = sd$tmax_internal,
    time_indices = list(
      idx_time_beta = sd$idx_time_beta,
      idx_time_uid = sd$idx_time_uid,
      idx_time_vmk = sd$idx_time_vmk,
      idx_time_idm = sd$idx_time_idm
    ),
    marker_weights = sd$marker_weights,
    fixed_marker_weights = sd$fixed_marker_weights,
    basehaz = sd$basehaz,
    n_knots = sd$n_knots,
    basehaz_degree = sd$basehaz_degree,
    K_event = sd$K_event,
    vcov_diag_link = sd$vcov_diag_link,
    use_tau_sde_fixed = sd$use_tau_sde_fixed,
    tau_sde_fixed = sd$tau_sde_fixed
  )
  cfg$engine <- engine

  fit_obj <- JoiNMeFit$new(
    fit = fit,
    stan_data = sd,
    formulaLong = formulaLong,
    formulaEvent = formulaEvent,
    formulaVCov = formulaVCov,
    config = cfg,
    call = match.call(),
    tmax = sd$tmax_internal,
    dataLong = dataLong,
    dataEvent = dataEvent
  )

  # Store a compact, self-contained plotting bundle so association plots remain
  # usable even when cmdstanr CSV outputs are no longer available.
  fit_obj$config$association_plot_payload <- tryCatch(
    .build_JoiNMefit_association_plot_payload(
      fit = fit,
      stan_data = sd,
      config = fit_obj$config,
      dataLong = dataLong,
      seed = args$seed %||% defaults$seed
    ),
    error = function(e) NULL
  )

  fit_obj
}
