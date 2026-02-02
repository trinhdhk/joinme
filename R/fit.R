#' Fit joinme model via cmdstanr or rstan
#'
#' @importFrom stats setNames
#' @importFrom utils modifyList
#'
#' @description
#' Fits the joinme Stan model using cmdstanr (default) or rstan. The engine can be
#' set globally via options(stan_preferred_engine = "cmdstanr"|"rstan") or via
#' control$engine.
#' following the notebook, and supports association transformations and time-internal
#' scaling metadata.
#'
#' @details
#' Transformations are provided through a single `transforms` argument.
#' Each element (cv_total, cs_total, vcov) is defined as a list specifying a type
#' and type-specific fields (see `build_standata_transforms()` for details).
#'
#' Distributional regression can be specified in two equivalent forms:
#' 1. A named list with RHS-only formulas, e.g. `list(sigma = ~ 1 + time)`.
#' 2. An unnamed list with explicit LHS parameter names, e.g. `list(sigma ~ 1 + time)`.
#'
#' The time axis is internally scaled for numerical stability; reported coefficients
#' are rescaled within Stan so fixed effects remain interpretable on the original scale.
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
#'      vcov = list(type = "pwlin", x = c(-2, 0, 2), y = c(0.2, 1, 0.2)))`
#'
#' Additional arguments are forwarded to `joinme_standata()` (e.g., `assoc`,
#' `basehaz`, `basehaz_degree`, `n_knots`, `time_var`, `shrinkage`, and `indep_*`).
#'
#' @param formulaLong Longitudinal formula defining fixed effects, id-level effects,
#'   and the marker block. The marker block may optionally include an inner
#'   `( ... | id )` term for marker-by-id random effects. When omitted,
#'   marker-by-id effects are disabled (Q_idm = 0).
#' @param dataLong Long-format longitudinal data with columns for id, marker, time,
#'   outcome, and covariates referenced in `formulaLong`.
#' @param formulaEvent Survival formula for baseline covariates and event model.
#' @param dataEvent One row per id event data with event time, event indicator, and
#'   covariates referenced in `formulaEvent`.
#' @param formulaVcov Covariance regression formula for id-specific marker-by-id effects.
#'   If the marker block omits the inner `( ... | id )`, marker-by-id effects are
#'   absent and `vcov` associations are not allowed.
#' @param formulaDist Optional list of formulas for distributional regression.
#'   Two forms are supported:
#'   1. Named list with RHS-only formulas, e.g. `list(sigma = ~ 1 + time)`.
#'   2. Unnamed list with LHS parameter names, e.g. `list(sigma ~ 1 + time)`.
#'   Random-effects terms with `|` are supported, but nested random-effects formulas
#'   are not.
#' @param control Named list containing sampling configuration.
#'   - cmdstanr::model$sample() arguments (e.g., chains, parallel_chains,
#'     iter_warmup, iter_sampling, seed, refresh, adapt_delta, max_treedepth).
#'   - engine: "cmdstanr" or "rstan". Defaults to options(stan_preferred_engine).
#'   - threads_per_chain: integer; if > 1 uses `joinme_fit_threading.stan`.
#'     For engine = "rstan", threading uses options(stan.thread = threads_per_chain).
#'   - grainsize: integer; reduce_sum grainsize for threading (default max(1, min(n_cores, ceiling(n_id/(4*threads_per_chain*chains))))).
#'   - force_recompile: logical; recompile the Stan model if needed.
#' @param draws Optional number of posterior draws used for summaries (not sampling).
#' @param families Marker-specific family specification (optional). Can be a vector
#'   of family names aligned to marker order.
#' @param transforms Transformation specifications for association terms
#'   (cv_total, cs_total, vcov). Each entry is a list with a `type` and fields
#'   required by that type.
#' @param priors Named list for priors: list(beta = ..., alpha = ..., lkj = ...).
#' @param ... Additional args passed to joinme_standata().
#'
#' @export
joinme <- function(
  formulaLong,
  dataLong,
  formulaEvent,
  dataEvent,
  formulaVcov = ~1,
  formulaDist = NULL,
  control = list(),
  draws = NULL,
  families = NULL,
  transforms = NULL,
  priors = list(beta = NULL, alpha = NULL, lkj = NULL),
  ...
) {
  if (!is.list(control)) {
    cli::cli_abort(c(
      x = "{.arg control} must be a named list.",
      i = "Provide list(threads_per_chain=..., grainsize=..., adapt_delta=..., max_treedepth=...)."
    ))
  }
  if (!is.list(priors)) {
    cli::cli_abort(c(
      x = "{.arg priors} must be a named list.",
      i = "Example: list(beta = NULL, alpha = NULL, lkj = NULL)."
    ))
  }

  if (length(control) > 0 && is.null(names(control))) {
    cli::cli_abort(c(
      x = "{.arg control} must be a named list.",
      i = "Provide list(threads_per_chain=..., grainsize=..., adapt_delta=..., max_treedepth=...)."
    ))
  }

  sd <- joinme_standata(
    formulaLong = formulaLong,
    dataLong = dataLong,
    formulaEvent = formulaEvent,
    dataEvent = dataEvent,
    formulaVcov = formulaVcov,
    formulaDist = formulaDist,
    families = families,
    transforms = transforms,
    beta_prior = priors$beta,
    alpha_prior = priors$alpha,
    lkj_prior = priors$lkj,
    ...
  )

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
  n_cores <- parallel::detectCores(logical = FALSE) %||% 1L
  max_threads <- min(sd$n_id %||% 1L, n_cores)
  if (threads_per_chain > max_threads) {
    cli::cli_warn(c(
      x = "Requested {threads_per_chain} threads exceeds max {max_threads}.",
      i = "Capping threads_per_chain to {max_threads}."
    ))
    threads_per_chain <- max_threads
  }

  grainsize <- control$grainsize
  if (is.null(grainsize)) {
    n_id <- sd$n_id %||% 1L
    n_chains <- control$parallel_chains %||% control$chains %||% 4L
    denom <- 4L * as.integer(threads_per_chain) * as.integer(n_chains)
    denom <- max(1L, denom)
    grainsize <- max(1L, as.integer(ceiling(n_id / denom)))
    grainsize <- min(as.integer(n_cores), grainsize)
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

  stan_candidates <- if (threads_per_chain > 1) {
    c(
      system.file("stan/joinme_fit_threading.stan", package = "joinme"),
      file.path("inst", "stan", "joinme_fit_threading.stan"),
      file.path("..", "inst", "stan", "joinme_fit_threading.stan"),
      file.path("..", "..", "inst", "stan", "joinme_fit_threading.stan")
    )
  } else {
    c(
      system.file("stan/joinme_fit.stan", package = "joinme"),
      file.path("inst", "stan", "joinme_fit.stan"),
      file.path("..", "inst", "stan", "joinme_fit.stan"),
      file.path("..", "..", "inst", "stan", "joinme_fit.stan")
    )
  }
  stan_file <- stan_candidates[file.exists(stan_candidates)][1]
  if (is.na(stan_file) || !nzchar(stan_file)) {
    cli::cli_abort(c(
      x = "Stan model file not found.",
      i = "Reinstall the package or restore inst/stan files."
    ))
  }
  use_threading <- threads_per_chain > 1

  engine <- .resolve_stan_engine(control$engine)
  if (engine == "rstan") {
    old_stan_thread <- getOption("stan.thread")
    options(stan.thread = threads_per_chain)
    on.exit(options(stan.thread = old_stan_thread), add = TRUE)
  }

  cpp_opts <- if (use_threading) list(stan_threads = TRUE) else NULL
  if (engine == "cmdstanr") {
    mod <- .get_cmdstan_model(
      stan_file,
      cpp_options = cpp_opts,
      force_recompile = control$force_recompile %||% FALSE
    )
  } else {
    mod <- .get_rstan_model(
      stan_file,
      force_recompile = control$force_recompile %||% FALSE
    )
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

  if (use_threading) {
    if (is.null(sd_stan$id) || is.null(sd_stan$n_id)) {
      cli::cli_abort(c(
        x = "Threading requires {.arg id} and {.arg n_id} in Stan data.",
        i = "Check the standata builder output."
      ))
    }
    id_vec <- sd_stan$id
    idx <- split(seq_along(id_vec), id_vec)
    id_start <- vapply(idx, min, integer(1))
    id_end <- vapply(idx, max, integer(1))
    if (length(id_start) != sd_stan$n_id) {
      cli::cli_abort(c(
        x = "Threading requires contiguous ids from 1..n_id.",
        i = "Check the id mapping in standata."
      ))
    }
    sd_stan$id_start <- as.integer(id_start)
    sd_stan$id_end <- as.integer(id_end)
    sd_stan$grainsize <- grainsize
  }

  if (engine == "rstan") {
    sd_stan <- .coerce_rstan_dist_arrays(sd_stan)
    sd_stan <- .coerce_rstan_time_indices(sd_stan)
    sd_stan <- .coerce_rstan_vectors(sd_stan, c(
      "beta_scale",
      "const_data_cv",
      "const_data_cs",
      "const_data_vcov"
    ))
    allowed_data <- .stan_data_names(stan_file)
    if (length(allowed_data) > 0) {
      missing <- setdiff(allowed_data, names(sd_stan))
      if (length(missing) > 0) {
        cli::cli_abort(c(
          x = "Stan data missing required fields for rstan: {paste(missing, collapse = ', ')}.",
          i = "Check joinme_standata() and Stan data block consistency."
        ))
      }
      sd_stan <- sd_stan[names(sd_stan) %in% allowed_data]
    }
  }

  defaults <- list(
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
    "force_recompile"
  ))]
  args <- modifyList(defaults, sample_control)
  chains_val <- args$parallel_chains %||% args$chains %||% defaults$parallel_chains
  args$parallel_chains <- chains_val
  args$chains <- chains_val
  if (use_threading) args$threads_per_chain <- threads_per_chain
  args$data <- sd_stan
  args <- args[!vapply(args, is.null, logical(1))]

  if (engine == "cmdstanr") {
    allowed <- names(formals(mod$sample))
    args <- args[names(args) %in% allowed]
    fit <- do.call(mod$sample, args)
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
      chains = args$chains %||% 1,
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
    assoc = c(
      cv_total = sd$assoc_cv_total,
      cv_mean = sd$assoc_cv_mean,
      cv_marker = sd$assoc_cv_marker,
      cs_total = sd$assoc_cs_total,
      cs_mean = sd$assoc_cs_mean,
      cs_marker = sd$assoc_cs_marker,
      vcov = sd$assoc_vcov
    ),
    transforms = list(
      tf_mode_cv_tot = sd$tf_mode_cv_tot,
      tf_mode_cs_tot = sd$tf_mode_cs_tot,
      tf_mode_cv_mean = sd$tf_mode_cv_mean,
      tf_mode_cs_mean = sd$tf_mode_cs_mean,
      tf_mode_cv_marker = sd$tf_mode_cv_marker,
      tf_mode_cs_marker = sd$tf_mode_cs_marker,
      tf_mode_vcov = sd$tf_mode_vcov
    ),
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
    tmax = sd$tmax,
    time_indices = list(
      idx_time_beta = sd$idx_time_beta,
      idx_time_uid = sd$idx_time_uid,
      idx_time_vmk = sd$idx_time_vmk,
      idx_time_widm = sd$idx_time_widm
    ),
    basehaz = sd$basehaz,
    n_knots = sd$n_knots,
    basehaz_degree = sd$basehaz_degree,
    K_event = sd$K_event
  )
  cfg$engine <- engine

  JoinMeFit$new(
    fit = fit,
    stan_data = sd,
    formulaLong = formulaLong,
    formulaEvent = formulaEvent,
    formulaVcov = formulaVcov,
    config = cfg,
    call = match.call(),
    tmax = sd$tmax
  )
}
