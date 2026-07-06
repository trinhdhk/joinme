#' Dynamic Prediction for Joint Models
#'
#' @importFrom stats predict median sd quantile na.omit optim terms
#' @importFrom utils head modifyList
#' @importFrom stats setNames
#' @keywords internal
#' @name predict.JoiNMeFit
NULL

# File overview:
# - Build subject-specific Stan data for dynamic prediction.
# - Run the dynpred model and summarise longitudinal/survival draws.

#' Dynamic prediction for a fitted `JoiNMeFit` model
#'
#' @rdname predict.JoiNMeFit
#' @description
#' Performs dynamic prediction for a new subject using a fitted `JoiNMeFit` joint model.
#' It calculates the posterior predictive distribution of longitudinal trajectories
#' and survival probabilities conditional on the subject's observed history up to a specific time point.
#'
#' @param object A fitted object of class `JoiNMeFit`.
#' @param newdataLong Data frame containing longitudinal histories for one or more subjects.
#' Must contain columns for id, time, marker, and response variables as specified in the original model formula.
#' @param newdataEvent Data frame containing one row per subject with event time, event status,
#' and baseline covariates used in the survival model.
#' @param process Character vector specifying which predictions to compute.
#' Options: "longitudinal" (future trajectory), "event" (conditional survival probability). Default: both.
#' @param pred_type Character. Type of longitudinal predictions:
#'   - "per_marker_id" (default): Subject-specific predictions for each marker.
#'   - "marginal_marker": Average across subjects for each marker (population-level by marker).
#'   - "marginal_id": Average across markers for each subject (subject-specific average).
#'   - "marker_subject": Alias for "marginal_marker".
#'   - "subject_marker": Alias for "marginal_id".
#' @param scale Character vector. Longitudinal prediction scale(s):
#'   - "epred": expected response (inverse-link),
#'   - "linpred": linear predictor,
#'   - "predict": predictive draw (includes noise).
#'   Defaults to all scales when not explicitly set.
#'   For `epred`/`predict`, family-specific inverse-links configured in
#'   `families` (via `jm_family(...)`) are respected.
#' @param times Numeric vector of times at which to predict the longitudinal and survival trajectories.
#' Can also be a named list of numeric vectors (one per subject id). If NULL, a grid
#' from `time_start` to `time_start + time_horizon` is generated (and truncated to
#' training support when needed). If fewer than 50 points are supplied for a subject,
#' a warning is emitted when `control$n_times = 50`.
#' @param time_start Numeric scalar or character. The conditioning time (last observation time).
#' If numeric, a single value is reused for all subjects, or a named vector supplies
#' subject-specific values. If character, it is interpreted as a column name in
#' `newdataEvent` holding subject-specific conditioning times. If NULL, defaults to
#' the maximum observed time in `newdataLong` for each subject.
#' @param time_horizon Numeric scalar. Prediction horizon (in time units) used when
#' `times` is NULL. The default is `tmax`.
#' @param ci_levels Numeric vector of credible interval levels for plotting.
#' Must be strictly between 0 and 1.
#' @param tmax Numeric scalar. The maximum time used for scaling during model fitting.
#' **CRITICAL**: If this is not provided and cannot be inferred from the object, predictions will be on the wrong time scale.
#' @param control Named list for prediction configuration.
#'   - cmdstanr::model$sample() arguments (e.g., `chains`, `parallel_chains`,
#'     `iter_warmup`, `iter_sampling`, `seed`, `refresh`, `adapt_delta`).
#'     Values override defaults except `data`.
#'   - engine: "cmdstanr" or "rstan". Defaults to options(stan_preferred_engine).
#'   - n_samples: integer; number of posterior parameter draws extracted from the
#'     fitted object before dynamic prediction. Default 200.
#'   - n_times: integer; number of points in the prediction time grid when
#'     `times` is `NULL`. Default 50.
#'   - n_pred_draws: integer; number of prediction draws produced by the dynpred
#'     object. This is independent of `n_samples` (posterior parameter draw
#'     extraction count). If omitted, defaults to the extracted posterior draw
#'     count.
#'   - threads_per_chain: integer; the threaded dynpred Stan program is always
#'     used. `threads_per_chain = 1` keeps execution serial while preserving the
#'     thread-capable kernel. For engine = "rstan", threading uses
#'     options(stan.thread = threads_per_chain).
#'   - grainsize: integer; reduce_sum grainsize for threaded prediction.
#'     Defaults to the full prediction draw count when `threads_per_chain = 1`,
#'     and to `max(1, ceiling(n_pred_draws/(4*threads_per_chain*chains)))`
#'     otherwise.
#'   - quadrature_nodes: optional positive integer target for total quadrature
#'     points during dynamic prediction. Allowed values are exactly `7`, `15`,
#'     `31`, `41`, `51`, and `61`. Only the node count is passed to Stan; GK
#'     nodes/weights are fixed in the Stan code.
#'   - progress: logical; show sampling progress bar (default TRUE).
#'   - initialisation fallback (cmdstanr): if a subject-level dynamic prediction
#'     chain fails to initialise, prediction automatically retries that subject
#'     with a narrower random init range (`init = 0.1`).
#' @param seed Integer. Random seed for reproducibility of random effect sampling.
#' @importFrom stats predict median sd quantile na.omit optim terms
#' @param ... Additional arguments (unused).
#' 
#' @return A list with two components:
#' \item{longitudinal}{A data.frame containing longitudinal predictions (Mean, Median, SD, 95% CrI) for each time point in `times`. Columns: id, time, marker, Estimate, Median, Est.Error, L95, U95.}
#' \item{longitudinal_fitted}{A data.frame containing fitted values for observed history on the linpred and epred scales. Includes a `scale` column.}
#' \item{survival}{A data.frame containing survival probabilities (S(t|time_start)) for each time point in `times`. Columns: id, time, Survival, Median, Est.Error, L95, U95.}
#' \item{cumhaz}{A data.frame containing conditional cumulative hazards H(t|time_start). Columns: id, time, Cumhaz, Median, Est.Error, L95, U95.}
#' \item{draws}{A list containing raw posterior draws (`longitudinal`, `longitudinal_fitted`, `survival`, `cumhaz`) and reconstructed subject-level random effects (`random_effects_id`, `random_effects_marker_id`, including marker-by-id covariance draws when id-dependent covariance is active).}
#'
#' @details
#' The function uses a Bayesian approach in two draw layers:
#' 1.  Extract posterior parameter draws from the fitted object (`n_samples`).
#' 2.  Re-index those draws to the dynpred draw count (`control$n_pred_draws`),
#'     allowing independent control of prediction Monte Carlo size.
#' 3.  For each dynpred draw, sample subject-specific random effects
#'     *conditional* on observed history in `newdataLong`.
#' 4.  Evaluate longitudinal and survival quantities on requested future grids.
#' 5.  Pool draws to form marginal predictive summaries and intervals.
#'
#' @export
predict.JoiNMeFit <- function(object,
                           newdataLong,
                           newdataEvent,
                           process = c("longitudinal", "event"),
                           pred_type = c("per_marker_id", "marginal_marker", "marginal_id", "marker_subject", "subject_marker"),
                           scale = c("epred", "linpred", "predict"),
                           times = NULL,
                           time_start = NULL,
                           time_horizon = NULL,
                           tmax = NULL,
                           ci_levels = c(0.5, 0.95),
                           control = list(),
                           seed = 123,
                           ...) {
    if (!inherits(object, "JoiNMeFit")) {
        cli::cli_abort(c(
            x = "Object must be a {.cls JoiNMeFit} fit.",
            i = "Fit the model with joinme() before predicting."
        ))
    }
    process <- match.arg(process, several.ok = TRUE)
    pred_type <- match.arg(pred_type)
    pred_type <- switch(pred_type,
                        marker_subject = "marginal_marker",
                        subject_marker = "marginal_id",
                        pred_type)
    scale <- .normalize_prediction_scales(scale)

    ci_levels <- .validate_ci_levels(ci_levels)
    quantile_probs <- .quantile_probs_from_ci(ci_levels)

    # Workflow: recover metadata -> build prediction data -> run dynpred -> summarise
    # 1. Recover Metadata
    meta <- .recover_metadata(object, tmax)
    tmax_val <- meta$tmax
    knots <- meta$knots
    col_means <- meta$col_means

    time_horizon_val <- if (is.null(time_horizon)) tmax_val else time_horizon
    if (!is.numeric(time_horizon_val) || length(time_horizon_val) != 1 ||
        !is.finite(time_horizon_val) || time_horizon_val <= 0) {
        cli::cli_abort(c(
            x = "{.arg time_horizon} must be a single positive finite numeric value.",
            i = "Set {.arg time_horizon} explicitly, or ensure {.arg tmax} can be recovered from the fitted object."
        ))
    }

    # 2. Parse Formulas
    forms <- .parse_formulas(object)
    sd <- object$stan_data

    # 3. Identify Subjects
    id_var <- eval(object$call$id_var) %||% "id"
    ids <- unique(newdataEvent[[id_var]])
    if (length(ids) == 0) {
        cli::cli_abort("No subjects found in {.code newdataEvent}.")
    }

    # Conditioning time inputs
    # - support scalar numeric, per-subject numeric, or column name in newdataEvent
    time_start_col <- NULL
    time_start_map <- NULL
    if (!is.null(time_start)) {
        if (is.character(time_start)) {
            if (length(time_start) != 1) {
                cli::cli_abort(c(
                    x = "{.arg time_start} must be a single column name when character.",
                    i = "Example: time_start = 'time_start'"
                ))
            }
            time_start_col <- time_start
            if (!(time_start_col %in% names(newdataEvent))) {
                cli::cli_abort(c(
                    x = "{.arg time_start} column {time_start_col} not found in {.code newdataEvent}.",
                    i = "Provide a valid column name or a numeric value."
                ))
            }
        } else if (is.numeric(time_start)) {
            if (length(time_start) > 1) {
                if (!is.null(names(time_start))) {
                    time_start_map <- time_start
                } else if (length(time_start) == length(ids)) {
                    time_start_map <- stats::setNames(as.numeric(time_start), ids)
                } else {
                    cli::cli_abort(c(
                        x = "{.arg time_start} must be a single numeric value or match the number of subjects.",
                        i = "Use a named vector to specify subject-specific values."
                    ))
                }
            }
        } else {
            cli::cli_abort(c(
                x = "{.arg time_start} must be numeric or a column name string.",
                i = "Example: time_start = 2.5 or time_start = 'time_start'."
            ))
        }
    }

    if (!is.list(control)) {
        cli::cli_abort(c(
            x = "{.arg control} must be a named list.",
            i = "Provide list(threads_per_chain=..., grainsize=..., iter_warmup=..., iter_sampling=...)."
        ))
    }

    if (length(control) > 0 && is.null(names(control))) {
        cli::cli_abort(c(
            x = "{.arg control} must be a named list.",
            i = "Provide list(threads_per_chain=..., grainsize=..., iter_warmup=..., iter_sampling=...)."
        ))
    }
    if ("data" %in% names(control)) {
        cli::cli_warn(c(
            x = "{.arg control$data} is ignored for prediction.",
            i = "Prediction always supplies subject-specific Stan data."
        ))
        control$data <- NULL
    }

    dots <- list(...)
    if ("n_samples" %in% names(dots)) {
        if (is.null(control$n_samples)) {
            control$n_samples <- dots$n_samples
        }
    }
    if ("n_times" %in% names(dots)) {
        if (is.null(control$n_times)) {
            control$n_times <- dots$n_times
        }
    }

    n_samples <- control$n_samples %||% 200L
    if (!is.numeric(n_samples) || length(n_samples) != 1 || !is.finite(n_samples)) {
        cli::cli_abort(c(
            "x" = "{.arg control$n_samples} must be a single finite numeric value.",
            "i" = "Use a positive integer, e.g. {.code control = list(n_samples = 200)}."
        ))
    }
    n_samples <- as.integer(n_samples)
    if (n_samples < 1L) {
        cli::cli_abort(c(
            "x" = "{.arg control$n_samples} must be >= 1.",
            "i" = "Use a positive integer for posterior draw extraction."
        ))
    }

    n_times <- control$n_times %||% 50L
    if (!is.numeric(n_times) || length(n_times) != 1 || !is.finite(n_times)) {
        cli::cli_abort(c(
            "x" = "{.arg control$n_times} must be a single finite numeric value.",
            "i" = "Use a positive integer, e.g. {.code control = list(n_times = 50)}."
        ))
    }
    n_times <- as.integer(n_times)
    if (n_times < 1L) {
        cli::cli_abort(c(
            "x" = "{.arg control$n_times} must be >= 1.",
            "i" = "Use a positive integer for prediction grid size."
        ))
    }

    threads_per_chain <- control$threads_per_chain %||%
        control$threads %||% control$mc.cores %||% 1L
    if (!is.numeric(threads_per_chain) || length(threads_per_chain) != 1) {
        cli::cli_abort(c(
            x = "{.arg threads_per_chain} must be a single numeric value.",
            i = "Use 1 to disable threading."
        ))
    }
    threads_per_chain <- as.integer(threads_per_chain)
    if (threads_per_chain < 1) {
        cli::cli_abort(c(
            x = "{.arg threads_per_chain} must be >= 1.",
            i = "Use 1 to disable threading."
        ))
    }

    # 4) Draw-layer setup
    #
    # Algorithm:
    # a) Extract posterior parameter draws from the fitted object (n_samples).
    # b) Apply time-scaling corrections to draw-dependent coefficients.
    # c) Resolve dynpred draw count from control$n_pred_draws (independent knob).
    # d) Re-index draw arrays so Stan receives exactly n_pred_draws rows.
    draws_list_raw <- .extract_draws_for_pred(object, n_samples, seed)
    draws_list_raw <- .scale_draw_dependent_time_terms(draws_list_raw, object$stan_data, tmax_val)
    n_samples_extracted <- .n_draws_in_prediction_list(draws_list_raw)
    n_pred_draws <- .resolve_n_pred_draws(control$n_pred_draws, n_samples_extracted)
    pred_draw_index <- .prediction_draw_index(n_samples_extracted, n_pred_draws, seed)
    draws_list <- .subset_draws_for_prediction(draws_list_raw, pred_draw_index)

    # Finalize threads_per_chain from control only, capped by available work.
    n_cores <- parallel::detectCores(logical = FALSE) %||% 1L
    max_threads <- min(n_pred_draws, n_cores)
    if (threads_per_chain > max_threads) {
        cli::cli_warn(c(
            x = "Requested {threads_per_chain} threads exceeds max {max_threads}.",
            i = "Capping threads_per_chain to {max_threads}."
        ))
        threads_per_chain <- max_threads
    }

    stan_file <- .get_stan_file(
        program = "JoiNMe_dynpred",
        threaded = TRUE
    )

    engine_default <- getOption("stan_preferred_engine", object$config$engine %||% "cmdstanr")
    engine <- .resolve_stan_engine(control$engine %||% engine_default)
    if (!is.null(object$config$engine) && !identical(engine, object$config$engine)) {
        cli::cli_warn(c(
            x = "Prediction engine {engine} differs from fitted engine {object$config$engine}.",
            i = "Set options(stan_preferred_engine=...) to align engines."
        ))
    }
    if (engine == "rstan") {
        old_stan_thread <- getOption("stan.thread")
        options(stan.thread = threads_per_chain)
        on.exit(options(stan.thread = old_stan_thread), add = TRUE)
    }

    cpp_opts <- list(stan_threads = TRUE)
    force_recompile <- isTRUE(control$force_recompile %||% getOption("JoiNMe.force_recompile", FALSE))
    if (engine == "cmdstanr") {
        mod <- .get_cmdstan_model(
            stan_file,
            cpp_options = cpp_opts,
            force_recompile = force_recompile
        )
    } else {
        mod <- .get_rstan_model(
            stan_file
        )
    }

    if (engine == "cmdstanr" && !.cmdstan_threads_enabled(mod)) {
        cli::cli_warn(c(
            x = "The CmdStan prediction model is not compiled with stan_threads.",
            i = "Recompiling the prediction model with stan_threads enabled."
        ))
        mod <- .get_cmdstan_model(
            stan_file,
            cpp_options = list(stan_threads = TRUE),
            force_recompile = TRUE
        )
        if (!.cmdstan_threads_enabled(mod)) {
            cli::cli_abort(c(
                x = "Thread-enabled compilation failed for the dynpred Stan program.",
                i = "Retry with {.code control = list(force_recompile = TRUE)} or precompile the threaded model ahead of time."
            ))
        }
    }

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

    # reduce_sum grainsize for prediction (draw-level parallelism)
    grainsize <- control$grainsize
    if (is.null(grainsize)) {
        n_pred_draws_local <- n_pred_draws
        if (threads_per_chain <= 1L) {
            grainsize <- as.integer(n_pred_draws_local)
        } else {
            # reduce_sum in dynpred parallelizes over posterior draws, not subjects.
            # Tune default grainsize using number of draws to avoid overly tiny slices
            # and to keep work balanced across threads/chains.
            n_chains <- control$chains %||% 1L
            denom <- 4L * as.integer(threads_per_chain) * as.integer(n_chains)
            denom <- max(1L, denom)
            grainsize <- max(1L, as.integer(ceiling(n_pred_draws_local / denom)))
            grainsize <- min(as.integer(n_cores), grainsize)
        }
    }
    if (!is.numeric(grainsize) || length(grainsize) != 1) {
        cli::cli_abort(c(
            x = "{.arg grainsize} must be a single numeric value.",
            i = "Use an integer >= 1."
        ))
    }
    grainsize <- as.integer(grainsize)
    if (grainsize < 1) {
        cli::cli_abort(c(
            x = "{.arg grainsize} must be >= 1.",
            i = "Use a positive integer for reduce_sum grainsize."
        ))
    }

    # 6. Loop
    # - iterate over subjects, build per-subject standata, run dynpred
    sample_control <- control[setdiff(names(control), c(
        "threads_per_chain",
        "threads",
        "mc.cores",
        "grainsize",
        "n_pred_draws",
        "n_samples",
        "n_times"
    ))]
    results_long <- list()
    results_surv <- list()
    results_cumhaz <- list()
    results_long_fit <- list()
    draws_long_list <- list()
    draws_surv_list <- list()
    draws_cumhaz_list <- list()
    draws_long_fit <- list()
    draws_re_id_list <- list()
    draws_re_marker_id_list <- list()
    quantiles_long_list <- list()
    quantiles_surv_list <- list()
    quantiles_cumhaz_list <- list()
    quantiles_long_fit <- list()
    pred_sampler_diag_list <- list()

    time_var <- eval(object$call$time_var) %||% "time"
    marker_var <- eval(object$call$marker_var) %||% "marker"

    time_start_by_id <- numeric(0)
    show_progress <- control$progress %||% TRUE
    if (show_progress) {
        has_progressr <- requireNamespace("progressr", quietly = TRUE)
    
        if (has_progressr) {
            progressr::handlers(control$progress_handler %||% progressr::handler_progress())
            pb <- progressr::progressor(along = ids)
        } else {
            pb <- utils::txtProgressBar(min = 0, max = length(ids), style = 3, width = 60)
        }
    }
    marker_corr_depends_on_id <- .marker_corr_depends_on_id(object, newdataEvent)

    .safe_progress(show_progress = show_progress, { 
        for (id in ids) {
            if (show_progress && has_progressr) {
                this_id <- as.character(id)
                pb(amount = 0, message = sprintf("Subject %s", this_id))
            }
            dE <- newdataEvent[newdataEvent[[id_var]] == id, , drop = FALSE]
            dL <- newdataLong[newdataLong[[id_var]] == id, , drop = FALSE]

            if (nrow(dL) == 0) {
                cli::cli_warn("No longitudinal data for ID {id}.")
                next
            }

            # T_cond
            # - conditioning time is last observed time by default
            if (!is.null(time_start_col)) {
                t_cond <- dE[[time_start_col]][1]
            } else if (!is.null(time_start_map)) {
                t_cond <- time_start_map[[as.character(id)]]
            } else if (!is.null(time_start)) {
                t_cond <- time_start
            } else {
                t_cond <- max(dL[[time_var]], na.rm = TRUE)
            }
            if (!is.finite(t_cond)) {
                cli::cli_abort(c(
                    x = "Conditioning time is not finite for subject {id}.",
                    i = "Check {.arg time_start} or the time column in {.code newdataLong}."
                ))
            }
            time_start_by_id[as.character(id)] <- t_cond

            # Grids
            # - longitudinal and survival grids can differ in density
            t_grid <- numeric(0)
            if ("longitudinal" %in% process) {
                t_grid <- .resolve_time_grid(
                    times,
                    id,
                    t_cond,
                    tmax_val,
                    time_horizon_val,
                    default_n = n_times,
                    min_points = 50,
                    kind = "longitudinal"
                )
            }

            t_surv_grid <- numeric(0)
            if ("event" %in% process) {
                t_surv_grid <- .resolve_time_grid(
                    times,
                    id,
                    t_cond,
                    tmax_val,
                    time_horizon_val,
                    default_n = n_times,
                    min_points = 50,
                    kind = "survival"
                )
            }

            # SD Prep
            # - assemble subject-specific data list for Stan
            grainsize_data <- grainsize
            sd_pred <- .prepare_subject_standata(
                dE, dL, object, tmax_val, knots, col_means,
                t_cond, t_grid, t_surv_grid,
                forms, draws_list, grainsize_data,
                control = control,
                degree = meta$degree
            )
            # Ensure time index arrays are preserved for cmdstanr JSON (avoid auto-unbox)
            sd_pred <- .coerce_rstan_time_indices(sd_pred)
            sd_pred <- .coerce_rstan_dist_arrays(sd_pred)
            sd_pred <- .coerce_rstan_vectors(sd_pred, c(
                "vec_cov_vcov",
                "const_data_cv",
                "const_data_cs",
                "const_data_corr",
                "const_data_vcov",
                "const_data_cv_mean",
                "const_data_cv_marker",
                "const_data_cs_mean",
                "const_data_cs_marker"
            ))
            if (is.null(sd_pred$n_draws) || length(sd_pred$n_draws) != 1) {
                sd_pred$n_draws <- as.integer(n_pred_draws)
            }
            if (length(required_data) > 0) {
                missing <- setdiff(required_data, names(sd_pred))
                if (length(missing) > 0) {
                    cli::cli_abort(c(
                        x = "Prediction Stan data missing required fields: {paste(missing, collapse = ', ')}.",
                        i = "Check prediction standata construction and Stan data block."
                    ))
                }
            }
            if (!is.null(allowed_data) && length(allowed_data) > 0) {
                keep_data <- union(allowed_data, required_data)
                sd_pred <- sd_pred[names(sd_pred) %in% keep_data]
            }

            # Run Stan
            # - 1 chain, 1 post-warmup iteration, draws taken from input arrays
            sample_args <- list(
                chains = 1,
                iter_warmup = 100,
                iter_sampling = 1, # Single pass over n_samples array
                fixed_param = FALSE,
                refresh = 0,
                show_messages = FALSE,
                seed = seed,
                adapt_delta = 0.8
            )
            sample_args <- utils::modifyList(sample_args, sample_control)
            sample_args$threads_per_chain <- threads_per_chain
            sample_args$data <- sd_pred
            sample_args <- sample_args[!vapply(sample_args, is.null, logical(1))]
            if (engine == "cmdstanr") {
                allowed <- names(formals(mod$sample))
                sample_args <- sample_args[names(sample_args) %in% allowed]
                fit_pred <- do.call(mod$sample, sample_args)
                pred_diag <- .JoiNMe_sampler_diagnostics(fit_pred)
                if (is.na(pred_diag$draws) || pred_diag$draws < 1) {
                    retry_args <- sample_args
                    retry_args$init <- 0.1
                    cli::cli_warn(c(
                        x = "Prediction sampling failed to initialise for subject {id}.",
                        i = "Retrying with narrower random init range ({.code init = 0.1})."
                    ))
                    fit_pred <- do.call(mod$sample, retry_args)
                    pred_diag <- .JoiNMe_sampler_diagnostics(fit_pred)
                }
            } else {
                control_list <- list()
                if (!is.null(sample_args$adapt_delta)) control_list$adapt_delta <- sample_args$adapt_delta
                if (!is.null(sample_args$max_treedepth)) control_list$max_treedepth <- sample_args$max_treedepth
                iter_warmup <- sample_args$iter_warmup
                iter_sampling <- sample_args$iter_sampling
                iter_total <- iter_warmup + iter_sampling
                rstan_args <- list(
                    object = mod,
                    data = sample_args$data,
                    chains = sample_args$chains %||% 1,
                    iter = iter_total,
                    warmup = iter_warmup,
                    seed = sample_args$seed,
                    refresh = sample_args$refresh,
                    cores = min(sample_args$chains %||% 1, parallel::detectCores(logical = FALSE) %||% 1)
                )
                if (isTRUE(sample_args$fixed_param)) {
                    rstan_args$algorithm <- "Fixed_param"
                }
                if (length(control_list) > 0) rstan_args$control <- control_list
                fit_pred <- do.call(rstan::sampling, rstan_args)
                pred_diag <- .JoiNMe_sampler_diagnostics(fit_pred)
            }
            pred_sampler_diag_list[[as.character(id)]] <- pred_diag

            draw_variables <- .prediction_draw_variables(scale, sd_pred)
            draws_mat <- suppressMessages(
                .get_draws_matrix(
                    fit_pred,
                    variables = draw_variables
                )
            )

            u_id_draws <- .reconstruct_subject_u_id_draws(
                draws_matrix = draws_mat,
                standata_subject = sd_pred,
                n_draws_target = n_pred_draws
            )
            if (!is.null(u_id_draws)) {
                draws_re_id_list[[as.character(id)]] <- list(
                    matrix = u_id_draws,
                    terms = object$stan_data$zid_cols %||% paste0("u_id[", seq_len(ncol(u_id_draws)), "]")
                )
            }

            marker_id_draws <- .reconstruct_subject_marker_id_draws(
                draws_matrix = draws_mat,
                standata_subject = sd_pred,
                n_draws_target = n_pred_draws,
                marker_levels = object$stan_data$marker_levels
            )
            if (!is.null(marker_id_draws)) {
                draws_re_marker_id_list[[as.character(id)]] <- marker_id_draws
            }

            if (nrow(dL) > 0) {
                # Align marker levels to fitted model ordering
                marker_levels <- object$stan_data$marker_levels %||% levels(dL[[marker_var]]) %||% sort(unique(dL[[marker_var]]))
                marker_values_obs <- as.character(dL[[marker_var]])
                unknown_markers <- setdiff(unique(marker_values_obs), marker_levels)
                if (length(unknown_markers) > 0) {
                    cli::cli_abort(c(
                        x = "Found marker level(s) in {.arg newdataLong} not seen during fitting: {.val {paste(unknown_markers, collapse = ', ')}}.",
                        i = "Use only training marker levels: {.val {paste(marker_levels, collapse = ', ')}}."
                    ))
                }
                marker_int_obs <- match(marker_values_obs, marker_levels)
                if (any(is.na(marker_int_obs))) {
                    cli::cli_abort(c(
                        x = "Markers in newdataLong don't match fitted model levels.",
                        i = "Ensure marker levels match the fitted model."
                    ))
                }

                fitted_scales <- intersect(c("linpred", "epred"), scale)
                fit_draws_for_id <- list(time = dL[[time_var]], marker_idx = marker_int_obs)
                fit_sum_rows <- list()
                fit_quant_rows <- list()

                if ("linpred" %in% fitted_scales) {
                    fit_lin_mat <- .extract_matrix_from_stan(draws_mat, "y_fit_linpred", nrow(dL), n_pred_draws)
                    fit_draws_for_id$linpred <- fit_lin_mat

                    fit_lin_sum <- .summarize_pred_long(fit_lin_mat, dL[[time_var]], marker_int_obs, marker_levels, id)
                    fit_lin_sum$scale <- "linpred"
                    fit_sum_rows[[length(fit_sum_rows) + 1L]] <- fit_lin_sum

                    fit_lin_quant <- .compute_quantiles_long(
                        fit_lin_mat, dL[[time_var]], marker_int_obs, marker_levels, id, probs = quantile_probs
                    )
                    fit_lin_quant$scale <- "linpred"
                    fit_quant_rows[[length(fit_quant_rows) + 1L]] <- fit_lin_quant
                }

                if ("epred" %in% fitted_scales) {
                    fit_epred_mat <- .extract_matrix_from_stan(draws_mat, "y_fit_epred", nrow(dL), n_pred_draws)
                    fit_draws_for_id$epred <- fit_epred_mat

                    fit_epred_sum <- .summarize_pred_long(fit_epred_mat, dL[[time_var]], marker_int_obs, marker_levels, id)
                    fit_epred_sum$scale <- "epred"
                    fit_sum_rows[[length(fit_sum_rows) + 1L]] <- fit_epred_sum

                    fit_epred_quant <- .compute_quantiles_long(
                        fit_epred_mat, dL[[time_var]], marker_int_obs, marker_levels, id, probs = quantile_probs
                    )
                    fit_epred_quant$scale <- "epred"
                    fit_quant_rows[[length(fit_quant_rows) + 1L]] <- fit_epred_quant
                }

                if (length(fit_draws_for_id) > 2L) {
                    draws_long_fit[[as.character(id)]] <- fit_draws_for_id
                }
                if (length(fit_sum_rows) > 0) {
                    results_long_fit[[as.character(id)]] <- do.call(rbind, fit_sum_rows)
                }
                if (length(fit_quant_rows) > 0) {
                    quantiles_long_fit[[as.character(id)]] <- do.call(rbind, fit_quant_rows)
                }
            }

            if (length(t_grid) > 0) {
                # Extract predictions for all (time, marker) combinations
                n_obs_pred_total <- nrow(sd_pred$mat_fixed_pred)  # This is n_time_points * n_markers
                n_markers_actual <- length(unique(dL[[marker_var]]))
                marker_levels <- object$stan_data$marker_levels %||% levels(dL[[marker_var]]) %||% sort(unique(dL[[marker_var]]))
                time_rep <- rep(t_grid, n_markers_actual)

                scale_var_map <- c(
                    epred = "y_pred_epred",
                    linpred = "y_pred_linpred",
                    predict = "y_pred"
                )

                per_id_draws <- list()
                per_id_results <- list()
                per_id_quantiles <- list()
                for (scale_i in scale) {
                    long_var <- unname(scale_var_map[scale_i])
                    long_mat <- .extract_matrix_from_stan(draws_mat, long_var, n_obs_pred_total, n_pred_draws)

                    per_id_draws[[scale_i]] <- list(
                        matrix = long_mat,
                        marker_idx = sd_pred$idx_marker_pred,
                        time = time_rep,
                        scale = scale_i
                    )

                    long_sum <- .summarize_pred_long(long_mat, time_rep, sd_pred$idx_marker_pred, marker_levels, id)
                    long_sum$scale <- scale_i
                    per_id_results[[length(per_id_results) + 1L]] <- long_sum

                    long_quant <- .compute_quantiles_long(
                        long_mat, time_rep, sd_pred$idx_marker_pred,
                        marker_levels, id, probs = quantile_probs
                    )
                    long_quant$scale <- scale_i
                    per_id_quantiles[[length(per_id_quantiles) + 1L]] <- long_quant
                }

                if (length(scale) == 1L) {
                    draws_long_list[[as.character(id)]] <- per_id_draws[[1L]]
                } else {
                    draws_long_list[[as.character(id)]] <- per_id_draws
                }
                results_long[[as.character(id)]] <- do.call(rbind, per_id_results)
                quantiles_long_list[[as.character(id)]] <- do.call(rbind, per_id_quantiles)
            }
            if (length(t_surv_grid) > 0) {
                surv_mat <- .extract_matrix_from_stan(draws_mat, "surv_prob", length(t_surv_grid), n_pred_draws)
                cumhaz_mat <- .extract_matrix_from_stan(draws_mat, "cumhaz_cond", length(t_surv_grid), n_pred_draws)

                draws_surv_list[[as.character(id)]] <- list(
                    matrix = surv_mat,
                    time = t_surv_grid
                )
                draws_cumhaz_list[[as.character(id)]] <- list(
                    matrix = cumhaz_mat,
                    time = t_surv_grid
                )

                results_surv[[as.character(id)]] <- .summarize_pred_surv(surv_mat, t_surv_grid, id)
                results_cumhaz[[as.character(id)]] <- .summarize_pred_cumhaz(cumhaz_mat, t_surv_grid, id)
                
                # Also compute rich quantiles for plotting
                quantiles_surv_list[[as.character(id)]] <- .compute_quantiles_surv(
                    surv_mat, t_surv_grid, id, probs = quantile_probs
                )
                quantiles_cumhaz_list[[as.character(id)]] <- .compute_quantiles_surv(
                    cumhaz_mat, t_surv_grid, id, probs = quantile_probs
                )
            }

            if (control$progress %||% TRUE) {
                if (has_progressr) {
                    pb()
                } else {
                    utils::setTxtProgressBar(pb, which(ids == id))
                }
            }
        }
    })

    # Compute population-level trajectories (average over subjects, keeping marker)
    quantiles_long_marker_pop <- NULL
    quantiles_long_overall_pop <- NULL
    results_long_aggregated <- NULL
    results_surv_aggregated <- NULL
    
    if (length(quantiles_long_list) > 0) {
        # Aggregate across subjects to get marker-specific population trajectories
        all_marker_preds <- do.call(rbind, quantiles_long_list)
        q_cols <- intersect(.quantile_colnames(quantile_probs), names(all_marker_preds))
        
        # Apply pred_type aggregation to results_long
        if (length(results_long) > 0) {
            results_long_combined <- do.call(rbind, results_long)
            
            if (pred_type == "marginal_marker") {
                # Average posterior summaries over subjects within each
                # marker-time cell to obtain a marker-specific marginal mean.
                results_long_aggregated <- tidytable::as_tidytable(results_long_combined) |>
                    tidytable::summarize(
                        Estimate = mean(Estimate, na.rm = TRUE),
                        Median = mean(Median, na.rm = TRUE),
                        Est.Error = mean(Est.Error, na.rm = TRUE),
                        L95 = mean(L95, na.rm = TRUE),
                        U95 = mean(U95, na.rm = TRUE),
                        .by = c(scale, time, marker)
                    ) |>
                    tidytable::mutate(id = "Marginal_Marker") |>
                    tidytable::relocate(id, scale, time, marker, Estimate, Median, Est.Error, L95, U95)
            } else if (pred_type == "marginal_id") {
                # Average posterior summaries over markers within each
                # subject-time cell to obtain a subject-specific marginal mean.
                results_long_aggregated <- tidytable::as_tidytable(results_long_combined) |>
                    tidytable::summarize(
                        Estimate = mean(Estimate, na.rm = TRUE),
                        Median = mean(Median, na.rm = TRUE),
                        Est.Error = mean(Est.Error, na.rm = TRUE),
                        L95 = mean(L95, na.rm = TRUE),
                        U95 = mean(U95, na.rm = TRUE),
                        .by = c(scale, id, time)
                    )
            }
            # else per_marker_id: keep as is (results_long will be used)
        }
        
        # Group by time and marker, average quantiles across subjects
        quantiles_long_marker_pop <- tidytable::as_tidytable(all_marker_preds) |>
            tidytable::summarize(
                tidytable::across(tidytable::any_of(q_cols), ~ mean(.x, na.rm = TRUE)),
                mean = mean(mean, na.rm = TRUE),
                sd = mean(sd, na.rm = TRUE),
                .by = c(scale, time, marker)
            )
        quantiles_long_marker_pop$id <- "Population_Marker"
        quantiles_long_marker_pop$marker_idx <- NA_integer_
        quantiles_long_marker_pop <- quantiles_long_marker_pop[
            , c("id", "scale", "time", "marker", q_cols, "mean", "sd", "marker_idx"), drop = FALSE
        ]
        
        # Overall population trajectory (average across subjects and markers)
        quantiles_long_overall_pop <- tidytable::as_tidytable(all_marker_preds) |>
            tidytable::summarize(
                tidytable::across(tidytable::any_of(q_cols), ~ mean(.x, na.rm = TRUE)),
                mean = mean(mean, na.rm = TRUE),
                sd = mean(sd, na.rm = TRUE),
                .by = c(scale, time)
            )
        quantiles_long_overall_pop$id <- "Population_Overall"
        quantiles_long_overall_pop$marker <- "All"
        quantiles_long_overall_pop$marker_idx <- NA_integer_
        quantiles_long_overall_pop <- quantiles_long_overall_pop[
            , c("id", "scale", "time", "marker", q_cols, "mean", "sd", "marker_idx"), drop = FALSE
        ]
    }

    # Prepare metadata
    # Extract variable names from forms for use in plotting
    formula_long <- object$formulaLong %||% forms$formulaLong
    resp_var <- tryCatch(
        all.vars(formula_long)[1],
        error = function(e) NA_character_
    )
    
    metadata <- list(
        conditioning_time = if (length(time_start_by_id) > 0) max(time_start_by_id, na.rm = TRUE) else NA_real_,
        conditioning_time_by_id = time_start_by_id,
        n_samples = n_pred_draws,
        n_posterior_draws = n_samples_extracted,
        n_pred_draws = n_pred_draws,
        n_subjects = length(ids),
        marker_corr_depends_on_id = marker_corr_depends_on_id,
        indep_id_re = sd$indep_id_re,
        indep_marker_re = sd$indep_marker_re,
        indep_idmarker_cov = sd$indep_idmarker_cov,
        family_links = if (!is.null(sd$link_names)) {
            ifelse(is.na(sd$link_names), "custom", sd$link_names)
        } else {
            vapply(sd$link_long %||% integer(0), .link_name_from_code, character(1))
        },
        sampler_diagnostics = .aggregate_sampler_diagnostics(pred_sampler_diag_list),
        pred_type = pred_type,
        scale = if (length(scale) == 1L) scale else scale[1],
        scales = scale,
        fitted_scales = intersect(c("linpred", "epred"), scale),
        ci_levels = ci_levels,
        quantile_levels = quantile_probs,
        time_max = tmax_val,
        time_horizon = time_horizon_val,
        formula_long = formula_long,
        response_var = resp_var,
        id_var = eval(object$call$id_var) %||% "id",
        time_var = eval(object$call$time_var) %||% "time",
        marker_var = eval(object$call$marker_var) %||% "marker"
    )

    JoiNMeDynPred$new(
        predictions = list(
            longitudinal = if (!is.null(results_long_aggregated)) results_long_aggregated else if (length(results_long) > 0) do.call(rbind, results_long) else NULL,
            longitudinal_fitted = if (length(results_long_fit) > 0) do.call(rbind, results_long_fit) else NULL,
            survival = if (!is.null(results_surv_aggregated)) results_surv_aggregated else if (length(results_surv) > 0) do.call(rbind, results_surv) else NULL,
            cumhaz = if (length(results_cumhaz) > 0) do.call(rbind, results_cumhaz) else NULL
        ),
        quantiles = list(
            longitudinal = if (length(quantiles_long_list) > 0) do.call(rbind, quantiles_long_list) else NULL,
            longitudinal_fitted = if (length(quantiles_long_fit) > 0) do.call(rbind, quantiles_long_fit) else NULL,
            longitudinal_marker_pop = quantiles_long_marker_pop,
            longitudinal_overall_pop = quantiles_long_overall_pop,
            survival = if (length(quantiles_surv_list) > 0) do.call(rbind, quantiles_surv_list) else NULL,
            cumhaz = if (length(quantiles_cumhaz_list) > 0) do.call(rbind, quantiles_cumhaz_list) else NULL
        ),
        draws = list(
            longitudinal = draws_long_list,
            longitudinal_fitted = draws_long_fit,
            random_effects_id = draws_re_id_list,
            random_effects_marker_id = draws_re_marker_id_list,
            survival = draws_surv_list,
            cumhaz = draws_cumhaz_list
        ),
        data = list(
            longitudinal = newdataLong,
            event = newdataEvent
        ),
        metadata = metadata,
        call = match.call(),
        tmax = tmax_val,
        n_samples = n_pred_draws
    )
}

# Return number of posterior rows in the prediction draw list.
.n_draws_in_prediction_list <- function(draws_list) {
    if (!is.null(draws_list$beta_fixed) && is.matrix(draws_list$beta_fixed)) {
        return(as.integer(nrow(draws_list$beta_fixed)))
    }
    cli::cli_abort("Unable to determine number of posterior draws for prediction.")
}

# Resolve prediction draw count independent of posterior extraction count.
.resolve_n_pred_draws <- function(n_pred_draws, n_available) {
    if (is.null(n_pred_draws)) {
        return(as.integer(n_available))
    }
    if (!is.numeric(n_pred_draws) || length(n_pred_draws) != 1 || !is.finite(n_pred_draws)) {
        cli::cli_abort(c(
            x = "{.arg control$n_pred_draws} must be a single finite numeric value.",
            i = "Use a positive integer, e.g. {.code control = list(n_pred_draws = 400)}."
        ))
    }
    n_pred_draws <- as.integer(n_pred_draws)
    if (n_pred_draws < 1L) {
        cli::cli_abort(c(
            x = "{.arg control$n_pred_draws} must be >= 1.",
            i = "Use a positive integer for dynpred output draw count."
        ))
    }
    n_pred_draws
}

.normalize_prediction_scales <- function(scale) {
    allowed_scales <- c("epred", "linpred", "predict")
    if (is.null(scale)) {
        return(allowed_scales)
    }
    unique(match.arg(scale, choices = allowed_scales, several.ok = TRUE))
}

.prediction_draw_variables <- function(scale, standata_subject) {
    scale_var_map <- c(
        epred = "y_pred_epred",
        linpred = "y_pred_linpred",
        predict = "y_pred"
    )
    fit_var_map <- c(
        epred = "y_fit_epred",
        linpred = "y_fit_linpred"
    )
    long_vars <- unname(scale_var_map[scale])
    fitted_vars <- unname(fit_var_map[intersect(scale, names(fit_var_map))])

    draw_variables <- unique(c(long_vars, "surv_prob", "cumhaz_cond", fitted_vars, "z_u"))
    if (as.integer(standata_subject$n_random_marker %||% 0L) > 0L) {
        draw_variables <- c(draw_variables, "z_v")
    }
    if (as.integer(standata_subject$n_random_marker_id %||% 0L) > 0L) {
        draw_variables <- c(draw_variables, "z_w_lat", "z_L")
    }
    unique(draw_variables)
}

# Build draw re-index map from extracted posterior draws to prediction draws.
.prediction_draw_index <- function(n_available, n_target, seed) {
    n_available <- as.integer(n_available)
    n_target <- as.integer(n_target)
    if (n_available < 1L || n_target < 1L) {
        cli::cli_abort("Both available and target draw counts must be >= 1.")
    }
    if (n_available == n_target) {
        return(seq_len(n_available))
    }

    # Use deterministic pseudo-random indexing for reproducibility.
    has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (has_seed) {
        old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
    } else {
        on.exit(rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
    }
    set.seed(as.integer(seed) + 7001L)
    sample.int(n_available, size = n_target, replace = n_target > n_available)
}

# Apply draw re-indexing to all draw-dependent containers.
.subset_draws_for_prediction <- function(draws_list, draw_index) {
    reindex_first_dim <- function(x, idx) {
        if (is.null(x)) return(NULL)
        if (is.matrix(x)) {
            if (nrow(x) == 0L) return(x)
            return(x[idx, , drop = FALSE])
        }
        if (is.array(x)) {
            dims <- dim(x)
            if (length(dims) < 1L) return(x)
            if (is.na(dims[1]) || dims[1] == 0L) return(x)
            if (length(dims) > 1L && any(dims[-1L] == 0L)) return(x)
            subs <- c(list(idx), rep(list(TRUE), length(dims) - 1L), list(drop = FALSE))
            return(do.call(`[`, c(list(x), subs)))
        }
        if (is.atomic(x) && is.vector(x) && length(x) >= max(idx)) {
            return(x[idx])
        }
        x
    }

    out <- draws_list
    for (nm in names(out)) {
        out[[nm]] <- reindex_first_dim(out[[nm]], draw_index)
    }
    out
}

#' Posterior Linear Predictor
#'
#' @rdname predict.JoiNMeFit
#' @description
#' Convenience wrapper around the [predict] method returning the posterior
#' linear predictor (linpred scale).
#'
#' @param object A fitted object of class `JoiNMeFit`.
#' @param ... Additional arguments passed to [predict].
#' 
#' @return A `JoiNMeDynPred` object with `metadata$scale = "linpred"` and
#'   `metadata$scales = "linpred"`.
#'
#' @export
posterior_linpred.JoiNMeFit <- function(object, ...) {
    call_ <- match.call()
    call_[[1]] <- quote(predict)
    call_$scale <- "linpred"
    eval(call_, parent.frame())
}

#' Posterior Expected Predictor
#'
#' @rdname predict.JoiNMeFit
#' @description
#' Convenience wrapper around the [predict] method returning the posterior
#' expected predictor (epred scale).
#'
#' @param object A fitted object of class `JoiNMeFit`.
#' @param ... Additional arguments passed to [predict].
#'
#' @return A `JoiNMeDynPred` object with `metadata$scale = "epred"` and
#'   `metadata$scales = "epred"`.
#' @export
posterior_epred.JoiNMeFit <- function(object, ...) {
    call_ <- match.call()
    call_[[1]] <- quote(predict)
    call_$scale <- "epred"
    eval(call_, parent.frame())
}

#' Posterior Predictive Draws
#' 
#' @rdname predict.JoiNMeFit
#' @description
#' Convenience wrapper around the [predict] method returning posterior
#' predictive draws (includes observation noise).
#'
#' @param object A fitted object of class `JoiNMeFit`.
#' @param ... Additional arguments passed to [predict].
#'
#' @return A `JoiNMeDynPred` object with `metadata$scale = "predict"` and
#'   `metadata$scales = "predict"`.
#' @export
posterior_predict.JoiNMeFit <- function(object, ...) {
    call_ <- match.call()
    call_[[1]] <- quote(predict)
    call_$scale <- "predict"
    eval(call_, parent.frame())
}

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

.extract_matrix_from_stan <- function(draws_mat, var_name, N_cols, N_rows) {
    # Stan output flattened across posterior rows.
    # We recover a matrix (N_rows x N_cols), where N_rows is the draw index
    # embedded in Stan variable indices and N_cols is the observation/time index.
    # To preserve conditional uncertainty (and avoid over-shrunk intervals),
    # pair each draw index with a posterior row instead of averaging rows.
    mat <- matrix(0, N_rows, N_cols)
    n_post_rows <- nrow(draws_mat)
    row_idx <- ((seq_len(N_rows) - 1L) %% max(1L, n_post_rows)) + 1L
    for (n in 1:N_cols) {
        nms <- paste0(var_name, "[", 1:N_rows, ",", n, "]")
        alt_nms <- paste0(var_name, "[", n, ",", 1:N_rows, "]")

        if (all(nms %in% colnames(draws_mat))) {
            cur <- draws_mat[row_idx, nms, drop = FALSE]
            mat[, n] <- as.numeric(diag(cur))
        } else if (all(alt_nms %in% colnames(draws_mat))) {
            cur <- draws_mat[row_idx, alt_nms, drop = FALSE]
            mat[, n] <- as.numeric(diag(cur))
        }
    }
    mat
}

.scale_draw_dependent_time_terms <- function(draws_list, stan_data, tmax) {
    if (is.null(draws_list) || is.null(stan_data)) {
        return(draws_list)
    }
    tmax_num <- suppressWarnings(as.numeric(tmax))
    if (!is.finite(tmax_num) || length(tmax_num) != 1L || abs(tmax_num - 1) < 1e-12) {
        return(draws_list)
    }

    .scale_matrix_cols <- function(mat, idx, scale_factor) {
        if (is.null(mat) || !is.matrix(mat) || ncol(mat) == 0) {
            return(mat)
        }
        idx <- as.integer(idx %||% integer(0))
        idx <- idx[is.finite(idx) & idx >= 1L & idx <= ncol(mat)]
        if (length(idx) == 0) {
            return(mat)
        }
        mat[, idx] <- mat[, idx, drop = FALSE] * scale_factor
        mat
    }

    draws_list$beta_fixed <- .scale_matrix_cols(
        draws_list$beta_fixed,
        stan_data$idx_time_beta,
        tmax_num
    )
    draws_list$tau_id <- .scale_matrix_cols(
        draws_list$tau_id,
        stan_data$idx_time_uid,
        tmax_num
    )
    draws_list$tau_marker <- .scale_matrix_cols(
        draws_list$tau_marker,
        stan_data$idx_time_vmk,
        tmax_num
    )
    draws_list
}

.recover_metadata <- function(object, tmax_arg) {
    tmax <- tmax_arg
    if (is.null(tmax)) {
        tmax <- object$config$tmax_internal %||% object$stan_data$tmax_internal
    }
    if (is.null(tmax)) {
        cli::cli_warn(c(
            x = "{.arg tmax} not provided or found; assuming 1.0.",
            i = "If training used time scaling, predictions may be incorrect."
        ))
        tmax <- 1.0
    }

    S_event <- object$stan_data$S_event
    n_knots <- object$config$n_knots %||% eval(object$call$n_knots) %||% 5
    degree <- object$config$basehaz_degree %||% eval(object$call$basehaz_degree) %||% 3
    basehaz_type <- object$config$basehaz %||% eval(object$call$basehaz) %||% "bs"

    probs <- (seq_len(n_knots)) / (n_knots + 1)
    knots <- as.numeric(quantile(S_event, probs = probs, names = FALSE, type = 7))
    knots <- unique(pmin(pmax(knots, 1e-6), 1 - 1e-6))

    # Compute col_means using the stored basis object if available (best practice)
    # This avoids ill-conditioned basis warnings from reconstructing on training data
    if (!is.null(object$config$Bs_obj)) {
        # Use the stored spline basis object from training
        B_raw <- as.matrix(predict(object$config$Bs_obj, newx = S_event))
    } else {
        # Fallback: reconstruct basis (suppress warnings about boundary knots)
        if (basehaz_type == "bs") {
            B_raw <- suppressWarnings(
                splines::bs(S_event, knots = knots, Boundary.knots = c(0, 1), degree = degree, intercept = TRUE)
            )
        } else if (basehaz_type == "ns") {
            B_raw <- suppressWarnings(
                splines::ns(S_event, knots = knots, Boundary.knots = c(0, 1), intercept = TRUE)
            )
        } else {
            # formula or other
            B_raw <- object$stan_data$Bs_event_c
        }
    }

    if (ncol(B_raw) != object$stan_data$Kbs) {
        stop(sprintf(
            "Could not recover baseline hazard basis. Training Kbs=%d, but prediction basis has %d columns. Ensure n_knots and degree are consistent.",
            object$stan_data$Kbs, ncol(B_raw)
        ))
    }

    col_means <- colMeans(B_raw - object$stan_data$Bs_event_c)

    list(
        tmax = tmax,
        knots = knots,
        col_means = col_means,
        degree = degree
    )
}

.parse_formulas <- function(object) {
    call_dist <- object$call$formulaDist
    if (is.null(call_dist) || identical(call_dist, quote(NULL))) {
        call_dist <- NULL
    }
    if (!is.null(call_dist) &&
        !is.list(call_dist) &&
        !inherits(call_dist, "formula") &&
        !is.character(call_dist)) {
        call_dist <- NULL
    }
    list(
        formulaLong = object$formulaLong,
        formulaEvent = object$formulaEvent,
        formulaVCov = object$formulaVCov,
        formulaDist = call_dist %||% object$config$dist$dist_formulas
    )
}

.extract_draws_for_pred <- function(object, n_samples, seed) {
    fit <- object$fit
    d <- suppressMessages(.get_draws_obj(fit))
    nd <- posterior::ndraws(d)
    if (n_samples < nd) {
        set.seed(seed)
        d <- suppressMessages(posterior::subset_draws(d, draw = sample.int(nd, n_samples)))
    }
    dmat <- suppressMessages(posterior::as_draws_matrix(d))

    sd <- object$stan_data
    if ((isTRUE(sd$assoc_corr) || isTRUE(sd$assoc_vcov)) && sd$Q_idm < 1) {
        cli::cli_abort(c(
            x = "Covariance-style association requires marker-by-id random effects (Q_idm > 0).",
            i = "Refit with an inner ( ... | id ) term inside the marker block, or drop {.arg corr}/{.arg vcov} from {.arg assoc}."
        ))
    }
    n <- nrow(dmat)

    get_mat <- function(nms) {
        nms <- nms[nms %in% colnames(dmat)]
        if (length(nms) == 0) return(matrix(0, n, 0))
        as.matrix(dmat[, nms, drop = FALSE])
    }

    get_col <- function(nm, default = 0) {
        if (nm %in% colnames(dmat)) return(as.numeric(dmat[, nm]))
        rep(default, n)
    }

    get_transform_coeff_draws <- function(eff_prefix, base_coeff, n_coeff, n_components = 1L) {
        n_coeff <- as.integer(n_coeff %||% 0L)
        n_components <- as.integer(n_components %||% 1L)
        if (n_coeff < 1L) {
            if (n_components > 1L) {
                return(array(0, dim = c(n, n_components, 0L)))
            }
            return(matrix(0, n, 0))
        }
        if (n_components > 1L) {
            eff_names <- as.vector(outer(seq_len(n_components), seq_len(n_coeff), function(m, j) paste0(eff_prefix, "[", m, ",", j, "]")))
            eff_names <- eff_names[eff_names %in% colnames(dmat)]
            if (length(eff_names) > 0L) {
                out <- array(0, dim = c(n, n_components, n_coeff))
                for (m in seq_len(n_components)) {
                    for (j in seq_len(n_coeff)) {
                        nm <- paste0(eff_prefix, "[", m, ",", j, "]")
                        if (nm %in% colnames(dmat)) {
                            out[, m, j] <- dmat[, nm]
                        }
                    }
                }
                return(out)
            }

            base_mat <- as.matrix(base_coeff %||% matrix(0, nrow = n_components, ncol = n_coeff))
            if (nrow(base_mat) == 0L) {
                base_mat <- matrix(0, nrow = n_components, ncol = n_coeff)
            }
            if (nrow(base_mat) < n_components) {
                base_mat <- rbind(base_mat, matrix(0, nrow = n_components - nrow(base_mat), ncol = ncol(base_mat)))
            }
            if (ncol(base_mat) < n_coeff) {
                base_mat <- cbind(base_mat, matrix(0, nrow = nrow(base_mat), ncol = n_coeff - ncol(base_mat)))
            }
            base_mat <- base_mat[seq_len(n_components), seq_len(n_coeff), drop = FALSE]
            out <- array(0, dim = c(n, n_components, n_coeff))
            for (m in seq_len(n_components)) {
                out[, m, ] <- matrix(rep(base_mat[m, ], times = n), nrow = n, byrow = TRUE)
            }
            return(out)
        }

        eff_names <- paste0(eff_prefix, "[", seq_len(n_coeff), "]")
        if (all(eff_names %in% colnames(dmat))) {
            return(get_mat(eff_names))
        }
        base_coeff <- as.numeric(base_coeff %||% rep(0, n_coeff))
        if (length(base_coeff) < n_coeff) {
            base_coeff <- c(base_coeff, rep(0, n_coeff - length(base_coeff)))
        }
        matrix(rep(base_coeff[seq_len(n_coeff)], times = n), nrow = n, byrow = TRUE)
    }

    get_transform_iota_draws <- function(eff_prefix, default = 0, n_components = 1L, n_iota = 1L) {
        n_components <- as.integer(n_components %||% 1L)
        n_iota <- as.integer(n_iota %||% 1L)
        if (n_iota < 1L) {
            return(matrix(default, nrow = n, ncol = 0L))
        }
        if (n_components > 1L || n_iota > 1L) {
            eff_names <- paste0(eff_prefix, "[", seq_len(max(0L, n_components * n_iota)), "]")
            if (all(eff_names %in% colnames(dmat))) {
                return(get_mat(eff_names))
            }
            if (n_components > 1L && n_iota <= 1L) {
                legacy_names <- paste0(eff_prefix, "[", seq_len(n_components), "]")
                if (all(legacy_names %in% colnames(dmat))) {
                    return(get_mat(legacy_names))
                }
            }
            return(matrix(default, nrow = n, ncol = max(0L, n_components * n_iota)))
        }

        if (eff_prefix %in% colnames(dmat)) {
            return(matrix(as.numeric(dmat[, eff_prefix]), ncol = 1L))
        }
        matrix(default, nrow = n, ncol = 1L)
    }

    # Prep Arrays (Renamed map)
    # beta -> beta_fixed
    beta_fixed <- get_mat(paste0("beta[", 1:sd$P, "]"))

    # tau_u -> tau_id
    tau_id <- get_mat(paste0("tau_u[", 1:sd$R_id, "]"))

    Lcorr_id <- array(0, dim = c(n, sd$R_id, sd$R_id))
    for (r in 1:sd$R_id) {
        for (c in 1:sd$R_id) {
            nm <- paste0("Lcorr_u[", r, ",", c, "]")
            if (nm %in% colnames(dmat)) Lcorr_id[, r, c] <- dmat[, nm]
        }
    }

    if (sd$R_mk > 0) {
        tau_marker <- get_mat(paste0("tau_v[", 1:sd$R_mk, "]"))
        Lcorr_marker <- array(0, dim = c(n, sd$R_mk, sd$R_mk))
        for (r in 1:sd$R_mk) {
            for (c in 1:sd$R_mk) {
                nm <- paste0("Lcorr_v[", r, ",", c, "]")
                if (nm %in% colnames(dmat)) Lcorr_marker[, r, c] <- dmat[, nm]
            }
        }
        B_cross <- array(0, dim = c(n, sd$Q_idm, sd$R_mk))
        if (sd$Q_idm > 0) {
            for (q in 1:sd$Q_idm) {
                for (r in 1:sd$R_mk) {
                    nm <- paste0("B_cross[", q, ",", r, "]")
                    if (nm %in% colnames(dmat)) B_cross[, q, r] <- dmat[, nm]
                }
            }
        }
    } else {
        tau_marker <- matrix(0, n, 0)
        Lcorr_marker <- array(0, dim = c(n, 0, 0))
        B_cross <- array(0, dim = c(n, sd$Q_idm, 0))
    }

    # Covariance Regression
    M_cov <- if (sd$Q_idm > 0) sd$Q_idm * (sd$Q_idm + 1) / 2 else 0
    if (sd$indep_idmarker_cov == 1) M_cov <- if (sd$Q_idm > 0) sd$Q_idm else 0

    beta_vcov_reg_flat <- array(0, dim = c(n, M_cov * sd$K_cov))
    if (M_cov > 0) {
        for (m in 1:M_cov) {
            for (k in 1:sd$K_cov) {
                nm <- paste0("beta_L[", m, ",", k, "]")
                idx <- (m - 1) * sd$K_cov + k
                if (nm %in% colnames(dmat)) beta_vcov_reg_flat[, idx] <- dmat[, nm]
            }
        }
    }

    bs_gamma_c <- array(0, dim = c(n, sd$K_event %||% 1L, sd$Kbs))
    for (k_ev in 1:(sd$K_event %||% 1L)) {
        for (j in 1:sd$Kbs) {
            nm <- paste0("bs_gamma_c[", k_ev, ",", j, "]")
            if (nm %in% colnames(dmat)) bs_gamma_c[, k_ev, j] <- dmat[, nm]
        }
    }

    gamma_hazard <- array(0, dim = c(n, sd$K_event %||% 1L, sd$p_w))
    for (k_ev in 1:(sd$K_event %||% 1L)) {
        for (j in 1:sd$p_w) {
            nm <- paste0("gamma_w[", k_ev, ",", j, "]")
            if (nm %in% colnames(dmat)) gamma_hazard[, k_ev, j] <- dmat[, nm]
        }
    }

    marker_weight_draws_by_term <- setNames(vector("list", length(.weighted_assoc_term_keys())), .weighted_assoc_term_keys())
    for (term_key in .weighted_assoc_term_keys()) {
        eff_names <- paste0(.marker_weight_var_prefix(term_key, effective = TRUE), "[", 1:sd$D, "]")
        base_names <- paste0(.marker_weight_var_prefix(term_key, effective = FALSE), "[", 1:sd$D, "]")

        if (all(eff_names %in% colnames(dmat))) {
            marker_weight_draws_by_term[[term_key]] <- get_mat(eff_names)
        } else if (all(base_names %in% colnames(dmat))) {
            marker_weight_draws_by_term[[term_key]] <- get_mat(base_names)
        } else if (isTRUE(as.integer(sd$shared_marker_weights %||% 1L) == 1L) && all(paste0("marker_weights_eff[", 1:sd$D, "]") %in% colnames(dmat))) {
            marker_weight_draws_by_term[[term_key]] <- get_mat(paste0("marker_weights_eff[", 1:sd$D, "]"))
        } else if (isTRUE(as.integer(sd$shared_marker_weights %||% 1L) == 1L) && all(paste0("marker_weights[", 1:sd$D, "]") %in% colnames(dmat))) {
            marker_weight_draws_by_term[[term_key]] <- get_mat(paste0("marker_weights[", 1:sd$D, "]"))
        } else {
            base_by_term <- sd$marker_weights_by_term %||% list()
            base_weights <- as.numeric(base_by_term[[term_key]] %||% sd$marker_weights %||% rep(1, sd$D))
            if (length(base_weights) == sd$D) {
                marker_weight_draws_by_term[[term_key]] <- matrix(rep(base_weights, each = n), nrow = n, byrow = TRUE)
            }
        }
    }

    n_family_sigma <- as.integer(sd$n_family_sigma %||% if ((sd$flag_resid_dim %||% 1L) == 0L) 1L else sd$D)
    n_family_nu <- as.integer(sd$n_family_nu %||% sd$D)
    n_family_phi <- as.integer(sd$n_family_phi %||% sd$D)
    n_family_alpha <- as.integer(sd$n_family_alpha %||% sd$D)
    n_family_phi_beta <- as.integer(sd$n_family_phi_beta %||% sd$D)
    n_family_tau_sde <- as.integer(sd$n_family_tau_sde %||% sd$D)

    sigma_family <- if (n_family_sigma > 0 && paste0("sigma_family[1]") %in% colnames(dmat)) {
        get_mat(paste0("sigma_family[", 1:n_family_sigma, "]"))
    } else if (sd$D > 0 && n_family_sigma == sd$D && "sigma_marker[1]" %in% colnames(dmat)) {
        get_mat(paste0("sigma_marker[", 1:sd$D, "]"))
    } else if (n_family_sigma == 1 && "sigma_y" %in% colnames(dmat)) {
        matrix(get_col("sigma_y", default = 0), ncol = 1)
    } else {
        matrix(0, n, n_family_sigma)
    }
    nu_family <- if (n_family_nu > 0 && paste0("nu_family[1]") %in% colnames(dmat)) {
        get_mat(paste0("nu_family[", 1:n_family_nu, "]"))
    } else if (sd$D > 0 && n_family_nu == sd$D && "nu_marker[1]" %in% colnames(dmat)) {
        get_mat(paste0("nu_marker[", 1:sd$D, "]"))
    } else {
        matrix(0, n, n_family_nu)
    }
    phi_family <- if (n_family_phi > 0 && paste0("phi_family[1]") %in% colnames(dmat)) {
        get_mat(paste0("phi_family[", 1:n_family_phi, "]"))
    } else if (sd$D > 0 && n_family_phi == sd$D && "phi_nb_marker[1]" %in% colnames(dmat)) {
        get_mat(paste0("phi_nb_marker[", 1:sd$D, "]"))
    } else {
        matrix(0, n, n_family_phi)
    }
    alpha_family <- if (n_family_alpha > 0 && paste0("alpha_family[1]") %in% colnames(dmat)) {
        get_mat(paste0("alpha_family[", 1:n_family_alpha, "]"))
    } else if (sd$D > 0 && n_family_alpha == sd$D && "alpha_skew_marker[1]" %in% colnames(dmat)) {
        get_mat(paste0("alpha_skew_marker[", 1:sd$D, "]"))
    } else {
        matrix(0, n, n_family_alpha)
    }
    phi_beta_family <- if (n_family_phi_beta > 0 && paste0("phi_beta_family[1]") %in% colnames(dmat)) {
        get_mat(paste0("phi_beta_family[", 1:n_family_phi_beta, "]"))
    } else if (sd$D > 0 && n_family_phi_beta == sd$D && "phi_beta_marker[1]" %in% colnames(dmat)) {
        get_mat(paste0("phi_beta_marker[", 1:sd$D, "]"))
    } else {
        matrix(0, n, n_family_phi_beta)
    }
    tau_sde_family <- if (n_family_tau_sde > 0 && paste0("tau_sde_family[1]") %in% colnames(dmat)) {
        get_mat(paste0("tau_sde_family[", 1:n_family_tau_sde, "]"))
    } else if (sd$D > 0 && n_family_tau_sde == sd$D && "tau_sde_marker[1]" %in% colnames(dmat)) {
        get_mat(paste0("tau_sde_marker[", 1:sd$D, "]"))
    } else {
        matrix(0, n, n_family_tau_sde)
    }

    corr_coef_vars <- grep("^alpha_corr_eff\\[", colnames(dmat), value = TRUE)
    if (length(corr_coef_vars) == 0) {
        corr_coef_vars <- grep("^alpha_corr\\[", colnames(dmat), value = TRUE)
    }
    if (length(corr_coef_vars) > 0) {
        corr_idx <- suppressWarnings(as.integer(sub("^alpha_corr(?:_eff)?\\[(\\d+)\\]$", "\\1", corr_coef_vars)))
        corr_coef_vars <- corr_coef_vars[order(corr_idx)]
    }
    corr_coef_raw <- if (length(corr_coef_vars) > 0) get_mat(corr_coef_vars) else matrix(0, n, 0)
    M_corr_assoc <- if (sd$Q_idm > 1) sd$Q_idm * (sd$Q_idm - 1) / 2 else 0
    corr_coef_padded <- if (M_corr_assoc > 0) {
        out <- matrix(0, nrow = n, ncol = M_corr_assoc)
        n_copy <- min(ncol(corr_coef_raw), M_corr_assoc)
        if (n_copy > 0) out[, seq_len(n_copy)] <- corr_coef_raw[, seq_len(n_copy), drop = FALSE]
        out
    } else {
        matrix(0, n, 0)
    }

    vcov_coef_vars <- grep("^alpha_vcov_eff\\[", colnames(dmat), value = TRUE)
    if (length(vcov_coef_vars) == 0) {
        vcov_coef_vars <- grep("^alpha_vcov\\[", colnames(dmat), value = TRUE)
    }
    if (length(vcov_coef_vars) > 0) {
        vcov_idx <- suppressWarnings(as.integer(sub("^alpha_vcov(?:_eff)?\\[(\\d+)\\]$", "\\1", vcov_coef_vars)))
        vcov_coef_vars <- vcov_coef_vars[order(vcov_idx)]
    }
    vcov_coef_raw <- if (length(vcov_coef_vars) > 0) get_mat(vcov_coef_vars) else matrix(0, n, 0)
    M_vcov_assoc <- if (sd$Q_idm > 0) {
        if (as.integer(sd$indep_idmarker_cov %||% 0L) == 1L) sd$Q_idm else (sd$Q_idm * (sd$Q_idm + 1)) / 2
    } else {
        0
    }
    vcov_coef_padded <- if (M_vcov_assoc > 0) {
        out <- matrix(0, nrow = n, ncol = M_vcov_assoc)
        n_copy <- min(ncol(vcov_coef_raw), M_vcov_assoc)
        if (n_copy > 0) out[, seq_len(n_copy)] <- vcov_coef_raw[, seq_len(n_copy), drop = FALSE]
        out
    } else {
        matrix(0, n, 0)
    }

    list(
        n_samples = n,
        beta_fixed = beta_fixed,
        tau_id = tau_id, Lcorr_id = Lcorr_id,
        tau_marker = tau_marker, Lcorr_marker = Lcorr_marker,
        B_cross = B_cross,
        beta_sigma = if (!is.null(sd$P_sigma) && sd$P_sigma > 0 && "beta_sigma[1]" %in% colnames(dmat)) {
            get_mat(paste0("beta_sigma[", 1:sd$P_sigma, "]"))
        } else {
            matrix(0, n, 0)
        },
        beta_nu = if (!is.null(sd$P_nu) && sd$P_nu > 0 && "beta_nu[1]" %in% colnames(dmat)) {
            get_mat(paste0("beta_nu[", 1:sd$P_nu, "]"))
        } else {
            matrix(0, n, 0)
        },
        beta_phi = if (!is.null(sd$P_phi) && sd$P_phi > 0 && "beta_phi[1]" %in% colnames(dmat)) {
            get_mat(paste0("beta_phi[", 1:sd$P_phi, "]"))
        } else {
            matrix(0, n, 0)
        },
        beta_alpha = if (!is.null(sd$P_alpha) && sd$P_alpha > 0 && "beta_alpha[1]" %in% colnames(dmat)) {
            get_mat(paste0("beta_alpha[", 1:sd$P_alpha, "]"))
        } else {
            matrix(0, n, 0)
        },
        beta_phi_beta = if (!is.null(sd$P_phi_beta) && sd$P_phi_beta > 0 && "beta_phi_beta[1]" %in% colnames(dmat)) {
            get_mat(paste0("beta_phi_beta[", 1:sd$P_phi_beta, "]"))
        } else {
            matrix(0, n, 0)
        },
        beta_tau_sde = if (!is.null(sd$P_tau_sde) && sd$P_tau_sde > 0 && "beta_tau_sde[1]" %in% colnames(dmat)) {
            get_mat(paste0("beta_tau_sde[", 1:sd$P_tau_sde, "]"))
        } else {
            matrix(0, n, 0)
        },
        alpha_vcov_reg = if (M_cov > 0) get_mat(paste0("alpha_L[", 1:M_cov, "]")) else matrix(0, n, 0),
        beta_vcov_reg_flat = beta_vcov_reg_flat,
        lambda_vcov_reg = if (M_cov > 0) get_mat(paste0("lambda_L[", 1:M_cov, "]")) else matrix(0, n, 0),
        bs_gamma_c = bs_gamma_c,
        gamma_hazard = gamma_hazard,
        marker_weights_draws = marker_weight_draws_by_term,
        sigma_family = sigma_family,
        nu_family = nu_family,
        phi_family = phi_family,
        alpha_family = alpha_family,
        phi_beta_family = phi_beta_family,
        tau_sde_family = tau_sde_family,
        coeff_assoc_cv_total = get_col("alpha_cv_total"),
        coeff_assoc_cs_total = get_col("alpha_cs_total"),
        coeff_assoc_cv_mean = get_col("alpha_cv_mean"),
        coeff_assoc_cs_mean = get_col("alpha_cs_mean"),
        coeff_assoc_cv_marker = get_col("alpha_cv_marker"),
        coeff_assoc_cs_marker = get_col("alpha_cs_marker"),
        coeff_assoc_corr = corr_coef_padded,
        coeff_assoc_vcov = vcov_coef_padded,
        iota_intercept_cv = get_transform_iota_draws("iota_intercept_cv_eff", default = 0, n_iota = sd$estimate_iota_intercept_cv %||% 1L),
        iota_slope_cv = get_transform_iota_draws("iota_slope_cv_eff", default = 1, n_iota = sd$estimate_iota_slope_cv %||% 1L),
        iota_intercept_cs = get_transform_iota_draws("iota_intercept_cs_eff", default = 0, n_iota = sd$estimate_iota_intercept_cs %||% 1L),
        iota_slope_cs = get_transform_iota_draws("iota_slope_cs_eff", default = 1, n_iota = sd$estimate_iota_slope_cs %||% 1L),
        iota_intercept_corr = get_transform_iota_draws("iota_intercept_corr_eff", default = 0, n_components = sd$M_corr_tf %||% ncol(corr_coef_padded), n_iota = sd$estimate_iota_intercept_corr %||% 1L),
        iota_slope_corr = get_transform_iota_draws("iota_slope_corr_eff", default = 1, n_components = sd$M_corr_tf %||% ncol(corr_coef_padded), n_iota = sd$estimate_iota_slope_corr %||% 1L),
        iota_intercept_vcov = get_transform_iota_draws("iota_intercept_vcov_eff", default = 0, n_components = sd$M_vcov_tf %||% ncol(vcov_coef_padded), n_iota = sd$estimate_iota_intercept_vcov %||% 1L),
        iota_slope_vcov = get_transform_iota_draws("iota_slope_vcov_eff", default = 1, n_components = sd$M_vcov_tf %||% ncol(vcov_coef_padded), n_iota = sd$estimate_iota_slope_vcov %||% 1L),
        iota_intercept_cv_mean = get_transform_iota_draws("iota_intercept_cv_mean_eff", default = 0, n_iota = sd$estimate_iota_intercept_cv_mean %||% 1L),
        iota_slope_cv_mean = get_transform_iota_draws("iota_slope_cv_mean_eff", default = 1, n_iota = sd$estimate_iota_slope_cv_mean %||% 1L),
        iota_intercept_cv_marker = get_transform_iota_draws("iota_intercept_cv_marker_eff", default = 0, n_iota = sd$estimate_iota_intercept_cv_marker %||% 1L),
        iota_slope_cv_marker = get_transform_iota_draws("iota_slope_cv_marker_eff", default = 1, n_iota = sd$estimate_iota_slope_cv_marker %||% 1L),
        iota_intercept_cs_mean = get_transform_iota_draws("iota_intercept_cs_mean_eff", default = 0, n_iota = sd$estimate_iota_intercept_cs_mean %||% 1L),
        iota_slope_cs_mean = get_transform_iota_draws("iota_slope_cs_mean_eff", default = 1, n_iota = sd$estimate_iota_slope_cs_mean %||% 1L),
        iota_intercept_cs_marker = get_transform_iota_draws("iota_intercept_cs_marker_eff", default = 0, n_iota = sd$estimate_iota_intercept_cs_marker %||% 1L),
        iota_slope_cs_marker = get_transform_iota_draws("iota_slope_cs_marker_eff", default = 1, n_iota = sd$estimate_iota_slope_cs_marker %||% 1L),
        coeff_cv = get_transform_coeff_draws("coeff_cv_eff", sd$coeff_cv, sd$n_coeff_cv),
        coeff_cs = get_transform_coeff_draws("coeff_cs_eff", sd$coeff_cs, sd$n_coeff_cs),
        coeff_corr = get_transform_coeff_draws("coeff_corr_eff", sd$coeff_corr, sd$n_coeff_corr, n_components = sd$M_corr_tf %||% ncol(corr_coef_padded)),
        coeff_vcov = get_transform_coeff_draws("coeff_vcov_eff", sd$coeff_vcov, sd$n_coeff_vcov, n_components = sd$M_vcov_tf %||% ncol(vcov_coef_padded)),
        coeff_cv_mean = get_transform_coeff_draws("coeff_cv_mean_eff", sd$coeff_cv_mean, sd$n_coeff_cv_mean),
        coeff_cv_marker = get_transform_coeff_draws("coeff_cv_marker_eff", sd$coeff_cv_marker, sd$n_coeff_cv_marker),
        coeff_cs_mean = get_transform_coeff_draws("coeff_cs_mean_eff", sd$coeff_cs_mean, sd$n_coeff_cs_mean),
        coeff_cs_marker = get_transform_coeff_draws("coeff_cs_marker_eff", sd$coeff_cs_marker, sd$n_coeff_cs_marker),
        cutpoints_ord = if (!is.null(sd$K_ord) && sd$K_ord > 1 && "cutpoints_ord[1]" %in% colnames(dmat)) {
            get_mat(paste0("cutpoints_ord[", 1:(sd$K_ord - 1), "]"))
        } else {
            matrix(0, n, max(1, (sd$K_ord %||% 2L) - 1))
        }
    )
}

.prepare_subject_standata <- function(dE, dL, object, tmax, knots, col_means, t_cond, t_grid, t_surv_grid, forms, draws_list, grainsize = NULL, control = NULL, degree = NULL) {
    id_var <- eval(object$call$id_var) %||% "id"
    time_var <- eval(object$call$time_var) %||% "time"
    marker_var <- eval(object$call$marker_var) %||% "marker"

    sd <- object$stan_data
    blueprints <- sd$design_blueprints %||% list()

    dL[[time_var]] <- dL[[time_var]] / tmax
    T_cond_scaled <- t_cond / tmax

    if (!.is_model_matrix_blueprint(blueprints$fixed)) {
        f_exp <- reformulas::expandDoubleVerts(forms$formulaLong)
        bars <- reformulas::findbars(f_exp)
        f_fix <- reformulas::nobars(f_exp)
        fixed_rhs <- stats::update(f_fix, . ~ .)
        fixed_rhs[[2]] <- NULL

        grp <- vapply(bars, function(b) .group_name_from_expr(b[[3]]), character(1))
        id_designs <- .bar_terms_to_rhs_list(bars[which(grp == id_var)])
        nested <- .extract_nested_marker_terms(forms$formulaLong, marker_var, id_var)
        marker_designs <- nested$mk_rhs_list
        idm_designs <- nested$idm_rhs_list
    } else {
        fixed_rhs <- blueprints$fixed
        id_designs <- blueprints$id %||% list()
        marker_designs <- blueprints$marker %||% list()
        idm_designs <- blueprints$idm %||% list()
    }

    mat_fixed_obs <- .mm(fixed_rhs, dL)
    mat_id_obs <- if (length(id_designs) > 0) {
        do.call(cbind, lapply(id_designs, function(rhs) .mm(rhs, dL)))
    } else {
        matrix(0, nrow(dL), sd$R_id)
    }
    mat_marker_obs <- if (length(marker_designs) > 0) {
        do.call(cbind, lapply(marker_designs, function(rhs) .mm(rhs, dL)))
    } else {
        matrix(0, nrow(dL), sd$R_mk)
    }
    mat_marker_id_obs <- if (length(idm_designs) > 0) {
        do.call(cbind, lapply(idm_designs, function(rhs) .mm(rhs, dL)))
    } else {
        matrix(0, nrow(dL), sd$Q_idm)
    }

    marker_levels_fit <- as.character(sd$marker_levels %||% levels(dL[[marker_var]]) %||% sort(unique(as.character(dL[[marker_var]]))))
    family_codes_fit <- as.integer(sd$family_long %||% integer(0))
    if (length(marker_levels_fit) != length(family_codes_fit)) {
        cli::cli_abort(c(
            x = "Fitted marker/family metadata is inconsistent.",
            i = "Cannot build family-scoped distributional prediction matrices."
        ))
    }

    dist_formulas <- .normalize_formula_dist(forms$formulaDist)
    family_by_row_obs <- .family_by_row_from_marker(
        marker_values = dL[[marker_var]],
        marker_levels = marker_levels_fit,
        family_codes = family_codes_fit
    )
    dist_sigma_obs <- .build_dist_matrix(dist_formulas$sigma, dL, family_by_row = family_by_row_obs)
    dist_nu_obs <- .build_dist_matrix(dist_formulas$nu, dL, family_by_row = family_by_row_obs)
    dist_phi_obs <- .build_dist_matrix(dist_formulas$phi, dL, family_by_row = family_by_row_obs)
    dist_alpha_obs <- .build_dist_matrix(dist_formulas$alpha, dL, family_by_row = family_by_row_obs)
    dist_phi_beta_obs <- .build_dist_matrix(dist_formulas$phi_beta, dL, family_by_row = family_by_row_obs)
    dist_tau_sde_obs <- .build_dist_matrix(dist_formulas$tau_sde, dL, family_by_row = family_by_row_obs)

    vec_cov_hazard <- .mm_event(forms$formulaEvent, dE)

    vcov_design <- .build_vcov_design(
        formulaVCov = forms$formulaVCov,
        dataEvent = dE,
        time_var = eval(object$call$time_var) %||% "time",
        context = "predict.JoiNMeFit()"
    )
    vec_cov_vcov <- vcov_design$Xcov

    quadrature_nodes <- control$quadrature_nodes %||% sd$quadrature_nodes %||% sd$n_gk %||% NULL
    quadrature_nodes_input <- quadrature_nodes %||% 15L
    quad_req <- .resolve_gk_request(nodes = quadrature_nodes_input)
    quad <- .gk_single_panel(rule = quad_req$rule)
    n_gk <- as.integer(quad$n_gk)
    gk_nodes <- quad$nodes

    u_cond <- T_cond_scaled * gk_nodes
    u_cond_fwd <- u_cond + sd$eps_fd

    .eval_on_times <- function(rhs_l, times) {
        if (length(rhs_l) == 0) {
            return(matrix(0, n_gk, 0))
        }
        dd <- dE[rep(1, n_gk), , drop = FALSE]
        dd[[time_var]] <- times
        do.call(cbind, lapply(rhs_l, function(rhs) .mm(rhs, dd)))
    }

    mat_fixed_gk_cond <- .eval_on_times(list(fixed_rhs), u_cond)
    mat_id_gk_cond <- if (length(id_designs) > 0) .eval_on_times(id_designs, u_cond) else matrix(0, n_gk, sd$R_id)
    mat_marker_gk_cond <- if (length(marker_designs) > 0) .eval_on_times(marker_designs, u_cond) else matrix(0, n_gk, sd$R_mk)
    mat_marker_id_gk_cond <- if (length(idm_designs) > 0) .eval_on_times(idm_designs, u_cond) else matrix(0, n_gk, sd$Q_idm)

    mat_fixed_gk_cond_fwd <- .eval_on_times(list(fixed_rhs), u_cond_fwd)
    mat_id_gk_cond_fwd <- if (length(id_designs) > 0) .eval_on_times(id_designs, u_cond_fwd) else matrix(0, n_gk, sd$R_id)
    mat_marker_gk_cond_fwd <- if (length(marker_designs) > 0) .eval_on_times(marker_designs, u_cond_fwd) else matrix(0, n_gk, sd$R_mk)
    mat_marker_id_gk_cond_fwd <- if (length(idm_designs) > 0) .eval_on_times(idm_designs, u_cond_fwd) else matrix(0, n_gk, sd$Q_idm)

    if (!is.null(object$config$Bs_obj)) {
        bs_basis <- as.matrix(predict(object$config$Bs_obj, newx = u_cond))
    } else {
        # Fallback for formula-based or if Bs_obj missing
        degree_val <- degree %||% 3
        bs_basis <- .make_basehaz_basis(u_cond, basis = "bs", knots = knots, degree = degree_val, boundary = c(0, 1))
    }
    mat_basis_gk_cond <- sweep(bs_basis, 2, col_means, "-")

    # Validate marker levels in newdataLong against fitted marker levels.
    #
    # Prediction requires marker design matrices and family vectors indexed on the
    # training marker space. New marker levels cannot be assigned random effects
    # or family parameters consistently, so we fail early with a clear message.
    marker_levels <- object$stan_data$marker_levels
    if (is.null(marker_levels)) {
        marker_levels <- levels(dL[[marker_var]]) %||% sort(unique(as.character(dL[[marker_var]])))
    }
    marker_levels <- as.character(marker_levels)
    marker_values <- as.character(dL[[marker_var]])
    marker_index_lookup <- setNames(seq_along(marker_levels), marker_levels)
    unknown_markers <- setdiff(unique(marker_values), marker_levels)
    if (length(unknown_markers) > 0) {
        cli::cli_abort(c(
            x = "Found marker level(s) in {.arg newdataLong} not seen during fitting: {.val {paste(unknown_markers, collapse = ', ')}}.",
            i = "Use only training marker levels: {.val {paste(marker_levels, collapse = ', ')}}."
        ))
    }

    # Get unique markers for this subject
    markers_subject <- unique(dL[[marker_var]])
    n_markers_subject <- length(markers_subject)
    
    # Create prediction grid: expand time grid to include all markers
    # - vectorized expansion avoids nested loops for speed
    if (length(t_grid) == 0 && length(t_surv_grid) > 0) {
        t_grid <- t_surv_grid
    }
    n_obs_pred_per_marker <- length(t_grid)
    n_obs_pred <- n_obs_pred_per_marker * n_markers_subject
    
    if (n_obs_pred > 0) {
        dl_pred <- dL[rep(1, n_obs_pred), , drop = FALSE]
        time_rep <- rep(t_grid, times = n_markers_subject)
        marker_rep <- rep(markers_subject, each = n_obs_pred_per_marker)
        marker_rep_chr <- as.character(marker_rep)
        idx_marker_pred_vec <- as.integer(marker_index_lookup[marker_rep_chr])
        if (any(is.na(idx_marker_pred_vec))) {
            cli::cli_abort(c(
                x = "Prediction marker indexing failed for marker(s): {.val {paste(sort(unique(marker_rep_chr[is.na(idx_marker_pred_vec)])), collapse = ', ')}}.",
                i = "Ensure prediction marker levels match fitted marker levels."
            ))
        }

        dl_pred[[time_var]] <- time_rep
        dl_pred[[marker_var]] <- factor(marker_rep_chr, levels = marker_levels)
        
        # Scale time for Stan model
        dl_pred_scaled <- dl_pred
        dl_pred_scaled[[time_var]] <- dl_pred_scaled[[time_var]] / tmax
        
        mat_fixed_pred <- .mm(fixed_rhs, dl_pred_scaled)
        mat_id_pred <- if (length(id_designs) > 0) {
            do.call(cbind, lapply(id_designs, function(rhs) .mm(rhs, dl_pred_scaled)))
        } else {
            matrix(0, n_obs_pred, sd$R_id)
        }
        mat_marker_pred <- if (length(marker_designs) > 0) {
            do.call(cbind, lapply(marker_designs, function(rhs) .mm(rhs, dl_pred_scaled)))
        } else {
            matrix(0, n_obs_pred, sd$R_mk)
        }
        mat_marker_id_pred <- if (length(idm_designs) > 0) {
            do.call(cbind, lapply(idm_designs, function(rhs) .mm(rhs, dl_pred_scaled)))
        } else {
            matrix(0, n_obs_pred, sd$Q_idm)
        }
        idx_marker_pred <- idx_marker_pred_vec

        family_by_row_pred <- .family_by_row_from_marker(
            marker_values = dl_pred_scaled[[marker_var]],
            marker_levels = marker_levels_fit,
            family_codes = family_codes_fit
        )
        dist_sigma_pred <- .build_dist_matrix(dist_formulas$sigma, dl_pred_scaled, family_by_row = family_by_row_pred)
        dist_nu_pred <- .build_dist_matrix(dist_formulas$nu, dl_pred_scaled, family_by_row = family_by_row_pred)
        dist_phi_pred <- .build_dist_matrix(dist_formulas$phi, dl_pred_scaled, family_by_row = family_by_row_pred)
        dist_alpha_pred <- .build_dist_matrix(dist_formulas$alpha, dl_pred_scaled, family_by_row = family_by_row_pred)
        dist_phi_beta_pred <- .build_dist_matrix(dist_formulas$phi_beta, dl_pred_scaled, family_by_row = family_by_row_pred)
        dist_tau_sde_pred <- .build_dist_matrix(dist_formulas$tau_sde, dl_pred_scaled, family_by_row = family_by_row_pred)
    } else {
        mat_fixed_pred <- matrix(0, 0, sd$P)
        mat_id_pred <- matrix(0, 0, sd$R_id)
        mat_marker_pred <- matrix(0, 0, sd$R_mk)
        mat_marker_id_pred <- matrix(0, 0, sd$Q_idm)
        idx_marker_pred <- integer(0)

        dist_sigma_pred <- list(P = dist_sigma_obs$P, X = matrix(0.0, 0, dist_sigma_obs$P), cols = dist_sigma_obs$cols)
        dist_nu_pred <- list(P = dist_nu_obs$P, X = matrix(0.0, 0, dist_nu_obs$P), cols = dist_nu_obs$cols)
        dist_phi_pred <- list(P = dist_phi_obs$P, X = matrix(0.0, 0, dist_phi_obs$P), cols = dist_phi_obs$cols)
        dist_alpha_pred <- list(P = dist_alpha_obs$P, X = matrix(0.0, 0, dist_alpha_obs$P), cols = dist_alpha_obs$cols)
        dist_phi_beta_pred <- list(P = dist_phi_beta_obs$P, X = matrix(0.0, 0, dist_phi_beta_obs$P), cols = dist_phi_beta_obs$cols)
        dist_tau_sde_pred <- list(P = dist_tau_sde_obs$P, X = matrix(0.0, 0, dist_tau_sde_obs$P), cols = dist_tau_sde_obs$cols)
    }

    n_times_surv <- length(t_surv_grid)
    mat_basis_gk_surv <- array(0, dim = c(n_times_surv, n_gk, sd$Kbs))
    mat_fixed_gk_surv <- array(0, dim = c(n_times_surv, n_gk, sd$P))
    mat_id_gk_surv <- array(0, dim = c(n_times_surv, n_gk, sd$R_id))
    mat_marker_gk_surv <- array(0, dim = c(n_times_surv, n_gk, sd$R_mk))
    mat_marker_id_gk_surv <- array(0, dim = c(n_times_surv, n_gk, sd$Q_idm))
    mat_fixed_gk_surv_fwd <- mat_fixed_gk_surv
    mat_id_gk_surv_fwd <- mat_id_gk_surv
    mat_marker_gk_surv_fwd <- mat_marker_gk_surv
    mat_marker_id_gk_surv_fwd <- mat_marker_id_gk_surv

    if (n_times_surv > 0) {
        for (s in 1:n_times_surv) {
            ts <- t_surv_grid[s] / tmax
            us <- ts * gk_nodes
            us_f <- us + sd$eps_fd

            if (!is.null(object$config$Bs_obj)) {
                bs <- as.matrix(predict(object$config$Bs_obj, newx = us))
            } else {
                degree_val <- degree %||% 3
                bs <- .make_basehaz_basis(us, basis = "bs", knots = knots, degree = degree_val, boundary = c(0, 1))
            }
            mat_basis_gk_surv[s, , ] <- sweep(bs, 2, col_means, "-")
            mat_fixed_gk_surv[s, , ] <- .eval_on_times(list(fixed_rhs), us)
            mat_id_gk_surv[s, , ] <- if (length(id_designs) > 0) .eval_on_times(id_designs, us) else matrix(0, n_gk, sd$R_id)
            mat_marker_gk_surv[s, , ] <- if (length(marker_designs) > 0) .eval_on_times(marker_designs, us) else matrix(0, n_gk, sd$R_mk)
            mat_marker_id_gk_surv[s, , ] <- if (length(idm_designs) > 0) .eval_on_times(idm_designs, us) else matrix(0, n_gk, sd$Q_idm)

            mat_fixed_gk_surv_fwd[s, , ] <- .eval_on_times(list(fixed_rhs), us_f)
            mat_id_gk_surv_fwd[s, , ] <- if (length(id_designs) > 0) .eval_on_times(id_designs, us_f) else matrix(0, n_gk, sd$R_id)
            mat_marker_gk_surv_fwd[s, , ] <- if (length(marker_designs) > 0) .eval_on_times(marker_designs, us_f) else matrix(0, n_gk, sd$R_mk)
            mat_marker_id_gk_surv_fwd[s, , ] <- if (length(idm_designs) > 0) .eval_on_times(idm_designs, us_f) else matrix(0, n_gk, sd$Q_idm)
        }
    }

    M_cov_val <- if (sd$Q_idm > 0) sd$Q_idm * (sd$Q_idm + 1) / 2 else 0
    if (sd$indep_idmarker_cov == 1) M_cov_val <- if (sd$Q_idm > 0) sd$Q_idm else 0

    ridx <- integer(M_cov_val)
    cidx <- integer(M_cov_val)
    if (M_cov_val > 0) {
        if (sd$indep_idmarker_cov == 1) {
            ridx <- 1:M_cov_val
            cidx <- 1:M_cov_val
        } else {
            idx <- 1
            for (r in 1:sd$Q_idm) {
                for (c in 1:r) {
                    ridx[idx] <- r
                    cidx[idx] <- c
                    idx <- idx + 1
                }
            }
        }
    }

    # Response and marker identification
    y_var <- .resolve_response_var(
        formulaLong = object$formulaLong,
        dataLong = dL,
        context = "predict.joinme()"
    )

    # Map markers carefully using fitted levels
    marker_levels <- object$stan_data$marker_levels
    if (is.null(marker_levels)) {
        marker_levels <- levels(dL[[marker_var]]) %||% sort(unique(as.character(dL[[marker_var]])))
    }

    # Marker weights (use fitted term-specific weights when available)
    marker_weights_by_term <- object$stan_data$marker_weights_by_term %||% list()
    marker_weights_base <- lapply(.weighted_assoc_term_keys(), function(term_key) {
        weights_term <- as.numeric(marker_weights_by_term[[term_key]] %||% object$stan_data$marker_weights %||% rep(1, length(marker_levels)))
        if (length(weights_term) != length(marker_levels) || any(!is.finite(weights_term))) {
            weights_term <- rep(1, length(marker_levels))
        }
        weights_term
    })
    names(marker_weights_base) <- .weighted_assoc_term_keys()

    marker_weights_draws <- draws_list$marker_weights_draws %||% list()
    n_pred_draws <- nrow(draws_list$beta_fixed)
    for (term_key in .weighted_assoc_term_keys()) {
        term_draws <- marker_weights_draws[[term_key]]
        if (is.null(term_draws) || nrow(term_draws) != n_pred_draws) {
            marker_weights_draws[[term_key]] <- matrix(
                rep(marker_weights_base[[term_key]], each = n_pred_draws),
                nrow = n_pred_draws,
                byrow = TRUE
            )
        }
    }

    marker_id_row_scale <- rep(1, sd$Q_idm)
    idx_time_idm <- as.integer(sd$idx_time_idm %||% sd$idx_time_widm %||% integer(0))
    idx_time_idm <- idx_time_idm[is.finite(idx_time_idm) & idx_time_idm >= 1L & idx_time_idm <= sd$Q_idm]
    if (length(idx_time_idm) > 0L) {
        marker_id_row_scale[idx_time_idm] <- tmax
    }

    dL[[marker_var]] <- factor(dL[[marker_var]], levels = marker_levels)
    marker_int <- as.integer(dL[[marker_var]])
    if (any(is.na(marker_int))) {
        unknown_markers <- sort(unique(as.character(dL[[marker_var]])[is.na(marker_int)]))
        cli::cli_abort(c(
            x = "Markers in {.arg newdataLong} do not match fitted model levels: {.val {paste(unknown_markers, collapse = ', ')}}.",
            i = "Ensure {.arg marker_var} levels match the training data."
        ))
    }

    # Trials (optional) for binomial predictions
    trials_obs <- rep.int(1L, nrow(dL))
    if ("trials" %in% colnames(dL)) {
        trials_obs <- as.integer(dL$trials)
    }
    trials_pred <- rep.int(1L, n_obs_pred)
    if (n_obs_pred > 0 && "trials" %in% colnames(dL)) {
        trial_map <- setNames(rep.int(1L, length(marker_levels)), marker_levels)
        for (m_val in unique(as.character(dL[[marker_var]]))) {
            m_trials <- dL$trials[as.character(dL[[marker_var]]) == m_val]
            if (length(m_trials) > 0) as.integer(m_trials[1]) else 1L
            trial_map[m_val] <- if (length(m_trials) > 0) as.integer(m_trials[1]) else 1L
        }
        trials_pred <- as.integer(trial_map[as.character(marker_levels[idx_marker_pred])])
    }

    .default_marker_family_map <- function(map_vec, n_family) {
        if (!is.null(map_vec)) {
            return(as.integer(map_vec))
        }
        if (sd$D <= 0) {
            return(integer(0))
        }
        if (n_family <= 0) {
            return(rep.int(0L, sd$D))
        }
        rep.int(1L, sd$D)
    }

    n_family_sigma_out <- as.integer(sd$n_family_sigma %||% ncol(draws_list$sigma_family))
    n_family_nu_out <- as.integer(sd$n_family_nu %||% ncol(draws_list$nu_family))
    n_family_phi_out <- as.integer(sd$n_family_phi %||% ncol(draws_list$phi_family))
    n_family_alpha_out <- as.integer(sd$n_family_alpha %||% ncol(draws_list$alpha_family))
    n_family_phi_beta_out <- as.integer(sd$n_family_phi_beta %||% ncol(draws_list$phi_beta_family))
    n_family_tau_sde_out <- as.integer(sd$n_family_tau_sde %||% ncol(draws_list$tau_sde_family))

    marker_to_sigma_family_out <- .default_marker_family_map(sd$marker_to_sigma_family, n_family_sigma_out)
    marker_to_nu_family_out <- .default_marker_family_map(sd$marker_to_nu_family, n_family_nu_out)
    marker_to_phi_family_out <- .default_marker_family_map(sd$marker_to_phi_family, n_family_phi_out)
    marker_to_alpha_family_out <- .default_marker_family_map(sd$marker_to_alpha_family, n_family_alpha_out)
    marker_to_phi_beta_family_out <- .default_marker_family_map(sd$marker_to_phi_beta_family, n_family_phi_beta_out)
    marker_to_tau_sde_family_out <- .default_marker_family_map(sd$marker_to_tau_sde_family, n_family_tau_sde_out)

    out <- list(
        n_draws = n_pred_draws,
        n_obs_long = nrow(dL), idx_marker_obs = marker_int, n_marker_types = sd$D,
        marker_weights_cv_total = as.numeric(marker_weights_base$cv_total),
        marker_weights_cs_total = as.numeric(marker_weights_base$cs_total),
        marker_weights_cv_marker = as.numeric(marker_weights_base$cv_marker),
        marker_weights_cs_marker = as.numeric(marker_weights_base$cs_marker),
        marker_weights_draws_cv_total = marker_weights_draws$cv_total,
        marker_weights_draws_cs_total = marker_weights_draws$cs_total,
        marker_weights_draws_cv_marker = marker_weights_draws$cv_marker,
        marker_weights_draws_cs_marker = marker_weights_draws$cs_marker,
        y_real = as.numeric(dL[[y_var]]),
        y_int = as.integer(dL[[y_var]]),
        trials_obs = as.integer(trials_obs),
        n_fixed_effects = sd$P, n_random_id = sd$R_id, n_random_marker = sd$R_mk, n_random_marker_id = sd$Q_idm,
        mat_fixed_obs = mat_fixed_obs, mat_id_obs = mat_id_obs, mat_marker_obs = mat_marker_obs, mat_marker_id_obs = mat_marker_id_obs,
        marker_id_row_scale = as.numeric(marker_id_row_scale),
        n_cov_vcov = as.integer(ncol(vec_cov_vcov)), vec_cov_vcov = array(as.numeric(vec_cov_vcov), dim = as.integer(ncol(vec_cov_vcov))),
        n_cov_hazard = as.integer(ncol(vec_cov_hazard)), vec_cov_hazard = array(as.numeric(vec_cov_hazard), dim = as.integer(ncol(vec_cov_hazard))),
        n_basehaz_basis = sd$Kbs, time_condition = T_cond_scaled,
        n_gk = as.integer(n_gk),
        mat_basis_gk_cond = mat_basis_gk_cond,
        mat_fixed_gk_cond = mat_fixed_gk_cond, mat_id_gk_cond = mat_id_gk_cond, mat_marker_gk_cond = mat_marker_gk_cond, mat_marker_id_gk_cond = mat_marker_id_gk_cond,
        mat_fixed_gk_cond_fwd = mat_fixed_gk_cond_fwd, mat_id_gk_cond_fwd = mat_id_gk_cond_fwd, mat_marker_gk_cond_fwd = mat_marker_gk_cond_fwd, mat_marker_id_gk_cond_fwd = mat_marker_id_gk_cond_fwd,
        eps_finite_diff = sd$eps_fd,
        P_sigma = as.integer(dist_sigma_obs$P),
        X_sigma_obs = dist_sigma_obs$X,
        X_sigma_pred = dist_sigma_pred$X,
        P_nu = as.integer(dist_nu_obs$P),
        X_nu_obs = dist_nu_obs$X,
        X_nu_pred = dist_nu_pred$X,
        P_phi = as.integer(dist_phi_obs$P),
        X_phi_obs = dist_phi_obs$X,
        X_phi_pred = dist_phi_pred$X,
        P_alpha = as.integer(dist_alpha_obs$P),
        X_alpha_obs = dist_alpha_obs$X,
        X_alpha_pred = dist_alpha_pred$X,
        P_phi_beta = as.integer(dist_phi_beta_obs$P),
        X_phi_beta_obs = dist_phi_beta_obs$X,
        X_phi_beta_pred = dist_phi_beta_pred$X,
        P_tau_sde = as.integer(dist_tau_sde_obs$P),
        X_tau_sde_obs = dist_tau_sde_obs$X,
        X_tau_sde_pred = dist_tau_sde_pred$X,
        n_obs_pred = n_obs_pred, idx_marker_pred = as.array(as.integer(idx_marker_pred)),
        mat_fixed_pred = mat_fixed_pred, mat_id_pred = mat_id_pred, mat_marker_pred = mat_marker_pred, mat_marker_id_pred = mat_marker_id_pred,
        trials_pred = as.integer(trials_pred),
        n_times_surv = n_times_surv, vec_time_surv = as.array(t_surv_grid / tmax),
        mat_basis_gk_surv = mat_basis_gk_surv,
        mat_fixed_gk_surv = mat_fixed_gk_surv, mat_id_gk_surv = mat_id_gk_surv, mat_marker_gk_surv = mat_marker_gk_surv, mat_marker_id_gk_surv = mat_marker_id_gk_surv,
        mat_fixed_gk_surv_fwd = mat_fixed_gk_surv_fwd, mat_id_gk_surv_fwd = mat_id_gk_surv_fwd, mat_marker_gk_surv_fwd = mat_marker_gk_surv_fwd, mat_marker_id_gk_surv_fwd = mat_marker_id_gk_surv_fwd,
        beta_fixed = draws_list$beta_fixed,
        tau_id = draws_list$tau_id, Lcorr_id = draws_list$Lcorr_id,
        tau_marker = draws_list$tau_marker, Lcorr_marker = draws_list$Lcorr_marker, B_cross = draws_list$B_cross,
        num_unique_cov_entries = M_cov_val,
        alpha_vcov_reg = draws_list$alpha_vcov_reg, beta_vcov_reg_flat = draws_list$beta_vcov_reg_flat,
        lambda_vcov_reg = draws_list$lambda_vcov_reg,
        K_event = sd$K_event %||% 1L,
        bs_gamma_c = draws_list$bs_gamma_c, gamma_hazard = draws_list$gamma_hazard,
        beta_sigma = draws_list$beta_sigma,
        beta_nu = draws_list$beta_nu,
        beta_phi = draws_list$beta_phi,
        beta_alpha = draws_list$beta_alpha,
        beta_phi_beta = draws_list$beta_phi_beta,
        beta_tau_sde = draws_list$beta_tau_sde,
        sigma_family = draws_list$sigma_family,
        nu_family = draws_list$nu_family,
        phi_family = draws_list$phi_family,
        alpha_family = draws_list$alpha_family,
        phi_beta_family = draws_list$phi_beta_family,
        tau_sde_family = draws_list$tau_sde_family,
        family_long = sd$family_long,
        link_long = sd$link_long,
        max_inv_link_ops = sd$max_inv_link_ops,
        inv_link_n_ops = sd$inv_link_n_ops,
        inv_link_ops = sd$inv_link_ops,
        max_inv_link_const = sd$max_inv_link_const,
        inv_link_n_const = sd$inv_link_n_const,
        inv_link_const = sd$inv_link_const,
        n_family_sigma = n_family_sigma_out,
        marker_to_sigma_family = marker_to_sigma_family_out,
        n_family_nu = n_family_nu_out,
        marker_to_nu_family = marker_to_nu_family_out,
        n_family_phi = n_family_phi_out,
        marker_to_phi_family = marker_to_phi_family_out,
        n_family_alpha = n_family_alpha_out,
        marker_to_alpha_family = marker_to_alpha_family_out,
        n_family_phi_beta = n_family_phi_beta_out,
        marker_to_phi_beta_family = marker_to_phi_beta_family_out,
        n_family_tau_sde = n_family_tau_sde_out,
        marker_to_tau_sde_family = marker_to_tau_sde_family_out,
        flag_resid_dim = sd$flag_resid_dim,
        vcov_diag_link = sd$vcov_diag_link,
        use_tau_sde_fixed = sd$use_tau_sde_fixed,
        tau_sde_fixed = sd$tau_sde_fixed,
        coeff_assoc_cv_total = draws_list$coeff_assoc_cv_total, coeff_assoc_cs_total = draws_list$coeff_assoc_cs_total,
        coeff_assoc_cv_mean = draws_list$coeff_assoc_cv_mean, coeff_assoc_cs_mean = draws_list$coeff_assoc_cs_mean,
        coeff_assoc_cv_marker = draws_list$coeff_assoc_cv_marker, coeff_assoc_cs_marker = draws_list$coeff_assoc_cs_marker,
        coeff_assoc_corr = draws_list$coeff_assoc_corr,
        coeff_assoc_vcov = draws_list$coeff_assoc_vcov,
        iota_intercept_cv = draws_list$iota_intercept_cv,
        iota_slope_cv = draws_list$iota_slope_cv,
        iota_intercept_cs = draws_list$iota_intercept_cs,
        iota_slope_cs = draws_list$iota_slope_cs,
        iota_intercept_corr = draws_list$iota_intercept_corr,
        iota_slope_corr = draws_list$iota_slope_corr,
        iota_intercept_vcov = draws_list$iota_intercept_vcov,
        iota_slope_vcov = draws_list$iota_slope_vcov,
        iota_intercept_cv_mean = draws_list$iota_intercept_cv_mean,
        iota_slope_cv_mean = draws_list$iota_slope_cv_mean,
        iota_intercept_cv_marker = draws_list$iota_intercept_cv_marker,
        iota_slope_cv_marker = draws_list$iota_slope_cv_marker,
        iota_intercept_cs_mean = draws_list$iota_intercept_cs_mean,
        iota_slope_cs_mean = draws_list$iota_slope_cs_mean,
        iota_intercept_cs_marker = draws_list$iota_intercept_cs_marker,
        iota_slope_cs_marker = draws_list$iota_slope_cs_marker,
        K_ord = sd$K_ord %||% 2L,
        cutpoints_ord = draws_list$cutpoints_ord,
        flag_assoc_cv_total = sd$assoc_cv_total, flag_assoc_cv_mean = sd$assoc_cv_mean, flag_assoc_cv_marker = sd$assoc_cv_marker,
        flag_assoc_cs_total = sd$assoc_cs_total, flag_assoc_cs_mean = sd$assoc_cs_mean, flag_assoc_cs_marker = sd$assoc_cs_marker,
        flag_assoc_corr = sd$assoc_corr, flag_assoc_vcov = sd$assoc_vcov,
        M_corr_tf = sd$M_corr_tf %||% ncol(draws_list$coeff_assoc_corr),
        M_vcov_tf = sd$M_vcov_tf %||% ncol(draws_list$coeff_assoc_vcov),
        estimate_iota_intercept_cv = sd$estimate_iota_intercept_cv %||% 0L,
        estimate_iota_slope_cv = sd$estimate_iota_slope_cv %||% 0L,
        estimate_iota_intercept_cs = sd$estimate_iota_intercept_cs %||% 0L,
        estimate_iota_slope_cs = sd$estimate_iota_slope_cs %||% 0L,
        estimate_iota_intercept_corr = sd$estimate_iota_intercept_corr %||% 0L,
        estimate_iota_slope_corr = sd$estimate_iota_slope_corr %||% 0L,
        estimate_iota_intercept_vcov = sd$estimate_iota_intercept_vcov %||% 0L,
        estimate_iota_slope_vcov = sd$estimate_iota_slope_vcov %||% 0L,
        estimate_iota_intercept_cv_mean = sd$estimate_iota_intercept_cv_mean %||% 0L,
        estimate_iota_slope_cv_mean = sd$estimate_iota_slope_cv_mean %||% 0L,
        estimate_iota_intercept_cv_marker = sd$estimate_iota_intercept_cv_marker %||% 0L,
        estimate_iota_slope_cv_marker = sd$estimate_iota_slope_cv_marker %||% 0L,
        estimate_iota_intercept_cs_mean = sd$estimate_iota_intercept_cs_mean %||% 0L,
        estimate_iota_slope_cs_mean = sd$estimate_iota_slope_cs_mean %||% 0L,
        estimate_iota_intercept_cs_marker = sd$estimate_iota_intercept_cs_marker %||% 0L,
        estimate_iota_slope_cs_marker = sd$estimate_iota_slope_cs_marker %||% 0L,
        tf_mode_cv_tot = sd$tf_mode_cv_tot, tf_mode_cs_tot = sd$tf_mode_cs_tot,
        tf_mode_cv_mean = sd$tf_mode_cv_mean, tf_mode_cv_marker = sd$tf_mode_cv_marker,
        tf_mode_cs_mean = sd$tf_mode_cs_mean, tf_mode_cs_marker = sd$tf_mode_cs_marker,
        tf_mode_corr = sd$tf_mode_corr, tf_mode_vcov = sd$tf_mode_vcov,
        n_functional_ops_cv = sd$n_functional_ops_cv, functional_ops_cv = sd$functional_ops_cv,
        functional_iota_intercept_idx_cv = sd$functional_iota_intercept_idx_cv %||% integer(sd$n_functional_ops_cv %||% 0L),
        functional_iota_slope_idx_cv = sd$functional_iota_slope_idx_cv %||% integer(sd$n_functional_ops_cv %||% 0L),
        n_const_cv = sd$n_const_cv, const_data_cv = sd$const_data_cv,
        n_functional_ops_cs = sd$n_functional_ops_cs, functional_ops_cs = sd$functional_ops_cs,
        functional_iota_intercept_idx_cs = sd$functional_iota_intercept_idx_cs %||% integer(sd$n_functional_ops_cs %||% 0L),
        functional_iota_slope_idx_cs = sd$functional_iota_slope_idx_cs %||% integer(sd$n_functional_ops_cs %||% 0L),
        n_const_cs = sd$n_const_cs, const_data_cs = sd$const_data_cs,
        n_functional_ops_corr = sd$n_functional_ops_corr, functional_ops_corr = sd$functional_ops_corr,
        functional_iota_intercept_idx_corr = sd$functional_iota_intercept_idx_corr %||% integer(sd$n_functional_ops_corr %||% 0L),
        functional_iota_slope_idx_corr = sd$functional_iota_slope_idx_corr %||% integer(sd$n_functional_ops_corr %||% 0L),
        n_const_corr = sd$n_const_corr, const_data_corr = sd$const_data_corr,
        n_functional_ops_vcov = sd$n_functional_ops_vcov, functional_ops_vcov = sd$functional_ops_vcov,
        functional_iota_intercept_idx_vcov = sd$functional_iota_intercept_idx_vcov %||% integer(sd$n_functional_ops_vcov %||% 0L),
        functional_iota_slope_idx_vcov = sd$functional_iota_slope_idx_vcov %||% integer(sd$n_functional_ops_vcov %||% 0L),
        n_const_vcov = sd$n_const_vcov, const_data_vcov = sd$const_data_vcov,
        n_knots_cv = sd$n_knots_cv, knots_cv = sd$knots_cv,
        n_coeff_cv = sd$n_coeff_cv, coeff_cv = draws_list$coeff_cv, spline_degree_cv = sd$spline_degree_cv,
        n_knots_cs = sd$n_knots_cs, knots_cs = sd$knots_cs,
        n_coeff_cs = sd$n_coeff_cs, coeff_cs = draws_list$coeff_cs, spline_degree_cs = sd$spline_degree_cs,
        n_knots_corr = sd$n_knots_corr, knots_corr = sd$knots_corr,
        n_coeff_corr = sd$n_coeff_corr, coeff_corr = draws_list$coeff_corr, spline_degree_corr = sd$spline_degree_corr,
        n_knots_vcov = sd$n_knots_vcov, knots_vcov = sd$knots_vcov,
        n_coeff_vcov = sd$n_coeff_vcov, coeff_vcov = draws_list$coeff_vcov, spline_degree_vcov = sd$spline_degree_vcov,
        n_functional_ops_cv_mean = sd$n_functional_ops_cv_mean, functional_ops_cv_mean = sd$functional_ops_cv_mean,
        functional_iota_intercept_idx_cv_mean = sd$functional_iota_intercept_idx_cv_mean %||% integer(sd$n_functional_ops_cv_mean %||% 0L),
        functional_iota_slope_idx_cv_mean = sd$functional_iota_slope_idx_cv_mean %||% integer(sd$n_functional_ops_cv_mean %||% 0L),
        n_const_cv_mean = sd$n_const_cv_mean, const_data_cv_mean = sd$const_data_cv_mean,
        n_knots_cv_mean = sd$n_knots_cv_mean, knots_cv_mean = sd$knots_cv_mean,
        n_coeff_cv_mean = sd$n_coeff_cv_mean, coeff_cv_mean = draws_list$coeff_cv_mean, spline_degree_cv_mean = sd$spline_degree_cv_mean,
        n_functional_ops_cv_marker = sd$n_functional_ops_cv_marker, functional_ops_cv_marker = sd$functional_ops_cv_marker,
        functional_iota_intercept_idx_cv_marker = sd$functional_iota_intercept_idx_cv_marker %||% integer(sd$n_functional_ops_cv_marker %||% 0L),
        functional_iota_slope_idx_cv_marker = sd$functional_iota_slope_idx_cv_marker %||% integer(sd$n_functional_ops_cv_marker %||% 0L),
        n_const_cv_marker = sd$n_const_cv_marker, const_data_cv_marker = sd$const_data_cv_marker,
        n_knots_cv_marker = sd$n_knots_cv_marker, knots_cv_marker = sd$knots_cv_marker,
        n_coeff_cv_marker = sd$n_coeff_cv_marker, coeff_cv_marker = draws_list$coeff_cv_marker, spline_degree_cv_marker = sd$spline_degree_cv_marker,
        n_functional_ops_cs_mean = sd$n_functional_ops_cs_mean, functional_ops_cs_mean = sd$functional_ops_cs_mean,
        functional_iota_intercept_idx_cs_mean = sd$functional_iota_intercept_idx_cs_mean %||% integer(sd$n_functional_ops_cs_mean %||% 0L),
        functional_iota_slope_idx_cs_mean = sd$functional_iota_slope_idx_cs_mean %||% integer(sd$n_functional_ops_cs_mean %||% 0L),
        n_const_cs_mean = sd$n_const_cs_mean, const_data_cs_mean = sd$const_data_cs_mean,
        n_knots_cs_mean = sd$n_knots_cs_mean, knots_cs_mean = sd$knots_cs_mean,
        n_coeff_cs_mean = sd$n_coeff_cs_mean, coeff_cs_mean = draws_list$coeff_cs_mean, spline_degree_cs_mean = sd$spline_degree_cs_mean,
        n_functional_ops_cs_marker = sd$n_functional_ops_cs_marker, functional_ops_cs_marker = sd$functional_ops_cs_marker,
        functional_iota_intercept_idx_cs_marker = sd$functional_iota_intercept_idx_cs_marker %||% integer(sd$n_functional_ops_cs_marker %||% 0L),
        functional_iota_slope_idx_cs_marker = sd$functional_iota_slope_idx_cs_marker %||% integer(sd$n_functional_ops_cs_marker %||% 0L),
        n_const_cs_marker = sd$n_const_cs_marker, const_data_cs_marker = sd$const_data_cs_marker,
        n_knots_cs_marker = sd$n_knots_cs_marker, knots_cs_marker = sd$knots_cs_marker,
        n_coeff_cs_marker = sd$n_coeff_cs_marker, coeff_cs_marker = draws_list$coeff_cs_marker, spline_degree_cs_marker = sd$spline_degree_cs_marker,
        flag_indep_id_re = sd$indep_id_re, flag_indep_marker_re = sd$indep_marker_re,
        flag_indep_idmarker_cov = sd$indep_idmarker_cov, flag_allow_marker_crosscorr = sd$allow_marker_crosscorr,
        idx_row_cov = as.array(as.integer(ridx)), idx_col_cov = as.array(as.integer(cidx))
    )

    if (!is.null(grainsize)) {
        out$grainsize <- as.integer(grainsize)
    }
    out
}

.summarize_pred_long <- function(mat, t_grid, marker_pred, levels, id) {
    # mat: n_samples x n_obs_pred
    # marker_pred: vector of marker indices (1-based)
    # levels: character vector of marker names
    sum_fun <- function(x) c(Mean = mean(x), Median = median(x), SD = sd(x), L95 = quantile(x, 0.025, names = FALSE), U95 = quantile(x, 0.975, names = FALSE))
    stats <- t(apply(mat, 2, sum_fun))
    
    if (is.null(levels) || length(levels) == 0) {
        levels <- as.character(sort(unique(marker_pred)))
    }
    marker_names <- if (length(marker_pred) > 0) levels[marker_pred] else NA_character_

    data.frame(
        id = id, 
        time = t_grid, 
        marker = marker_names,
        Estimate = stats[, "Mean"], 
        Median = stats[, "Median"], 
        Est.Error = stats[, "SD"], 
        L95 = stats[, "L95"], 
        U95 = stats[, "U95"]
    )
}

# Compute per-draw median survival time from S(t) trajectories.
#
# Algorithm:
# 1. Find the first index where S(t) <= threshold for each draw.
# 2. If no crossing occurs on the grid, median survival is NA for that draw.
# 3. If crossing is between adjacent points, linearly interpolate crossing time.
# 4. If first point is already below threshold, use first time point.
.median_survival_time_by_draw <- function(survival_draw_matrix, time_grid, threshold = 0.5) {
    if (is.null(survival_draw_matrix) || nrow(survival_draw_matrix) == 0 || ncol(survival_draw_matrix) == 0) {
        return(numeric(0))
    }
    if (length(time_grid) != ncol(survival_draw_matrix)) {
        return(rep(NA_real_, nrow(survival_draw_matrix)))
    }

    vapply(seq_len(nrow(survival_draw_matrix)), function(draw_index) {
        survival_values <- as.numeric(survival_draw_matrix[draw_index, ])
        if (all(!is.finite(survival_values))) return(NA_real_)

        crossing_index <- which(survival_values <= threshold)[1]
        if (is.na(crossing_index)) {
            finite_idx <- which(is.finite(survival_values) & survival_values > 0)
            if (length(finite_idx) < 2) return(NA_real_)

            tail_idx <- utils::tail(finite_idx, 2)
            t_left <- as.numeric(time_grid[tail_idx[1]])
            t_right <- as.numeric(time_grid[tail_idx[2]])
            s_left <- as.numeric(survival_values[tail_idx[1]])
            s_right <- as.numeric(survival_values[tail_idx[2]])
            if (!is.finite(t_left) || !is.finite(t_right) || t_right <= t_left) return(NA_real_)
            if (!is.finite(s_left) || !is.finite(s_right) || s_left <= 0 || s_right <= 0) return(NA_real_)
            if (s_right >= s_left) return(NA_real_)

            slope <- (log(s_right) - log(s_left)) / (t_right - t_left)
            if (!is.finite(slope) || slope >= 0 || abs(slope) < .Machine$double.eps) return(NA_real_)

            t_star <- t_right + (log(threshold) - log(s_right)) / slope
            if (!is.finite(t_star)) return(NA_real_)
            return(t_star)
        }
        if (crossing_index <= 1) return(as.numeric(time_grid[1]))

        time_left <- as.numeric(time_grid[crossing_index - 1])
        time_right <- as.numeric(time_grid[crossing_index])
        survival_left <- as.numeric(survival_values[crossing_index - 1])
        survival_right <- as.numeric(survival_values[crossing_index])

        if (!is.finite(survival_left) || !is.finite(survival_right) || abs(survival_right - survival_left) < .Machine$double.eps) {
            return(time_right)
        }

        time_left + (threshold - survival_left) * (time_right - time_left) / (survival_right - survival_left)
    }, numeric(1))
}

# Summarise one draw vector with fit-style diagnostics columns.
.summarize_draw_vector_with_diagnostics <- function(draw_values) {
    draw_values <- as.numeric(draw_values)
    draw_values <- draw_values[is.finite(draw_values)]

    if (length(draw_values) == 0) {
        return(data.frame(
            Estimate = NA_real_,
            Est.Error = NA_real_,
            Q2.5 = NA_real_,
            Q97.5 = NA_real_,
            Rhat = NA_real_,
            ess_bulk = NA_real_,
            ess_tail = NA_real_,
            stringsAsFactors = FALSE
        ))
    }

    q_bounds <- stats::quantile(draw_values, probs = c(0.025, 0.975), names = FALSE)
    rhat_val <- suppressWarnings(tryCatch(as.numeric(posterior::rhat(draw_values)), error = function(e) NA_real_))
    ess_bulk_val <- suppressWarnings(tryCatch(as.numeric(posterior::ess_basic(draw_values)), error = function(e) NA_real_))
    ess_tail_val <- suppressWarnings(tryCatch(as.numeric(posterior::ess_tail(draw_values)), error = function(e) NA_real_))

    data.frame(
        Estimate = mean(draw_values),
        Est.Error = stats::sd(draw_values),
        Q2.5 = q_bounds[1],
        Q97.5 = q_bounds[2],
        Rhat = rhat_val,
        ess_bulk = ess_bulk_val,
        ess_tail = ess_tail_val,
        stringsAsFactors = FALSE
    )
}

# Reconstruct realized id random effects u_id from latent z_u draws.
.reconstruct_subject_u_id_draws <- function(draws_matrix, standata_subject, n_draws_target) {
    n_random_id <- as.integer(standata_subject$n_random_id %||% 0L)
    if (n_random_id <= 0) return(NULL)

    z_u_draws <- .extract_matrix_from_stan(draws_matrix, "z_u", n_random_id, n_draws_target)
    if (is.null(z_u_draws) || nrow(z_u_draws) == 0 || ncol(z_u_draws) != n_random_id) return(NULL)

    tau_id_draws <- standata_subject$tau_id
    lcorr_id_draws <- standata_subject$Lcorr_id
    if (is.null(tau_id_draws) || is.null(lcorr_id_draws)) return(NULL)

    is_indep <- as.integer(standata_subject$flag_indep_id_re %||% standata_subject$indep_id_re %||% 0L) == 1L
    u_id_draws <- matrix(NA_real_, nrow = n_draws_target, ncol = n_random_id)

    for (draw_index in seq_len(n_draws_target)) {
        tau_vec <- as.numeric(tau_id_draws[draw_index, ])
        z_vec <- as.numeric(z_u_draws[draw_index, ])

        if (is_indep) {
            u_id_draws[draw_index, ] <- tau_vec * z_vec
        } else {
            lcorr_mat <- lcorr_id_draws[draw_index, , ]
            l_u <- diag(tau_vec, nrow = n_random_id, ncol = n_random_id) %*% lcorr_mat
            u_id_draws[draw_index, ] <- as.numeric(l_u %*% z_vec)
        }
    }

    u_id_draws
}

.reconstruct_subject_marker_id_draws <- function(draws_matrix, standata_subject, n_draws_target, marker_levels = NULL) {
    n_random_marker_id <- as.integer(standata_subject$n_random_marker_id %||% 0L)
    n_marker_types <- as.integer(standata_subject$n_marker_types %||% 0L)
    if (n_random_marker_id <= 0L || n_marker_types <= 0L) {
        return(NULL)
    }

    alpha_vcov_reg <- standata_subject$alpha_vcov_reg
    beta_vcov_reg_flat <- standata_subject$beta_vcov_reg_flat
    lambda_vcov_reg <- standata_subject$lambda_vcov_reg
    vec_cov_vcov <- as.numeric(standata_subject$vec_cov_vcov %||% numeric(0))
    idx_row_cov <- as.integer(standata_subject$idx_row_cov %||% integer(0))
    idx_col_cov <- as.integer(standata_subject$idx_col_cov %||% integer(0))

    if (is.null(alpha_vcov_reg) || is.null(beta_vcov_reg_flat) || is.null(lambda_vcov_reg) ||
        length(idx_row_cov) == 0 || length(idx_col_cov) == 0) {
        return(NULL)
    }

    n_random_marker <- as.integer(standata_subject$n_random_marker %||% 0L)
    tau_marker_draws <- standata_subject$tau_marker
    lcorr_marker_draws <- standata_subject$Lcorr_marker
    b_cross_draws <- standata_subject$B_cross

    is_indep_marker <- as.integer(standata_subject$flag_indep_marker_re %||% 0L) == 1L
    allow_marker_crosscorr <- as.integer(standata_subject$flag_allow_marker_crosscorr %||% 0L) == 1L

    n_post_rows <- nrow(draws_matrix)
    row_idx <- ((seq_len(n_draws_target) - 1L) %% max(1L, n_post_rows)) + 1L

    get_draw_scalar <- function(var_name, draw_index) {
        if (!(var_name %in% colnames(draws_matrix))) return(NA_real_)
        as.numeric(draws_matrix[row_idx[draw_index], var_name])
    }

    marker_labels <- as.character(marker_levels %||% seq_len(n_marker_types))
    if (length(marker_labels) != n_marker_types) {
        marker_labels <- as.character(seq_len(n_marker_types))
    }
    marker_id_terms <- standata_subject$zidm_cols %||% paste0("w_idm[", seq_len(n_random_marker_id), "]")
    marker_id_row_scale <- as.numeric(standata_subject$marker_id_row_scale %||% rep(1, n_random_marker_id))
    if (length(marker_id_row_scale) != n_random_marker_id || any(!is.finite(marker_id_row_scale))) {
        marker_id_row_scale <- rep(1, n_random_marker_id)
    }

    out_matrix <- matrix(NA_real_, nrow = n_draws_target, ncol = n_marker_types * n_random_marker_id)
    out_corr <- array(NA_real_, dim = c(n_draws_target, n_random_marker_id, n_random_marker_id))
    col_names <- unlist(lapply(marker_labels, function(marker_name) {
        paste0(marker_name, "::", marker_id_terms)
    }), use.names = FALSE)
    colnames(out_matrix) <- col_names

    m_cov <- length(idx_row_cov)
    z_l_draws <- .extract_matrix_from_stan(draws_matrix, "z_L", m_cov, n_draws_target)
    if (is.null(z_l_draws) || nrow(z_l_draws) == 0L || ncol(z_l_draws) != m_cov) {
        z_l_draws <- matrix(0, nrow = n_draws_target, ncol = m_cov)
    }

    for (draw_index in seq_len(n_draws_target)) {
        l_v <- NULL
        if (n_random_marker > 0 && !is.null(tau_marker_draws) && !is.null(lcorr_marker_draws)) {
            tau_v <- as.numeric(tau_marker_draws[draw_index, ])
            if (is_indep_marker) {
                l_v <- diag(tau_v, nrow = n_random_marker, ncol = n_random_marker)
            } else {
                l_v <- diag(tau_v, nrow = n_random_marker, ncol = n_random_marker) %*%
                    lcorr_marker_draws[draw_index, , ]
            }
        }

        lp_vec <- vapply(seq_len(m_cov), function(m) {
            b_slice_start <- (m - 1L) * length(vec_cov_vcov) + 1L
            b_slice_end <- m * length(vec_cov_vcov)
            b_vec <- if (length(vec_cov_vcov) > 0) {
                as.numeric(beta_vcov_reg_flat[draw_index, b_slice_start:b_slice_end])
            } else {
                numeric(0)
            }
            as.numeric(alpha_vcov_reg[draw_index, m]) +
                if (length(vec_cov_vcov) > 0) sum(b_vec * vec_cov_vcov) else 0 +
                as.numeric(lambda_vcov_reg[draw_index, m]) * as.numeric(z_l_draws[draw_index, m])
        }, numeric(1))

        l_i <- .cov_lp_to_chol(
            lp_vec = lp_vec,
            q_idm = n_random_marker_id,
            idx_row = idx_row_cov,
            idx_col = idx_col_cov,
            diag_link = standata_subject$vcov_diag_link
        )

        l_i_eff <- sweep(l_i, 1, marker_id_row_scale, `*`)

        draw_values <- numeric(n_marker_types * n_random_marker_id)
        col_offset <- 0L
        for (marker_index in seq_len(n_marker_types)) {
            z_w_lat <- vapply(seq_len(n_random_marker_id), function(q_index) {
                get_draw_scalar(
                    paste0("z_w_lat[", draw_index, ",", marker_index, ",", q_index, "]"),
                    draw_index
                )
            }, numeric(1))
            z_w_lat[!is.finite(z_w_lat)] <- 0

            cross <- rep(0, n_random_marker_id)
            if (allow_marker_crosscorr && n_random_marker > 0 && !is.null(l_v) && !is.null(b_cross_draws)) {
                z_v <- vapply(seq_len(n_random_marker), function(r_index) {
                    get_draw_scalar(
                        paste0("z_v[", draw_index, ",", marker_index, ",", r_index, "]"),
                        draw_index
                    )
                }, numeric(1))
                z_v[!is.finite(z_v)] <- 0
                v_marker <- as.numeric(l_v %*% z_v)
                b_cross_mat <- matrix(
                    as.numeric(b_cross_draws[draw_index, , , drop = FALSE]),
                    nrow = n_random_marker_id,
                    ncol = n_random_marker,
                    byrow = FALSE
                )
                if (ncol(b_cross_mat) == length(v_marker)) {
                    cross <- as.numeric(b_cross_mat %*% v_marker)
                }
            }

            z_w <- cross + z_w_lat
            w_idm <- as.numeric(l_i_eff %*% z_w)

            idx <- (col_offset + 1L):(col_offset + n_random_marker_id)
            draw_values[idx] <- w_idm
            col_offset <- col_offset + n_random_marker_id
        }

        out_matrix[draw_index, ] <- draw_values
        out_corr[draw_index, , ] <- l_i %*% t(l_i)
    }

    list(
        matrix = out_matrix,
        corr = out_corr,
        markers = marker_labels,
        terms = marker_id_terms,
        n_random_marker_id = n_random_marker_id
    )
}

.marker_corr_depends_on_id <- function(object, newdataEvent = NULL) {
    sd <- object$stan_data
    q_idm <- as.integer(sd$Q_idm %||% 0L)
    if (is.null(sd) || q_idm <= 0L) {
        return(FALSE)
    }

    # If marker-by-id random-effects are present in the fitted model,
    # prediction-time marker covariance varies by subject through L_i.
    # This is the primary availability condition for ranef/corr extraction.
    TRUE
}

#' Summarise dynamic prediction outputs
#'
#' @description
#' Summarises a `JoiNMeDynPred` object with detailed subject-level prediction
#' summaries and diagnostics.
#'
#' The summary includes:
#' - an overview table of row/subject counts by process,
#' - median survival time per subject (draw-wise crossing of `S(t)=0.5`,
#'   summarised with estimate, uncertainty, and diagnostics),
#' - predicted id-level random effects per subject (estimate, uncertainty,
#'   interval, and diagnostics),
#' - predicted marker-by-id random effects and covariance per subject when
#'   marker covariance depends on id,
#' - a compact diagnostics table counting potential convergence/ESS issues.
#'
#' @param object A `JoiNMeDynPred` object produced by [predict].
#' @param ... Unused.
#'
#' @return A `summary_JoiNMeDynPred` object containing tabular summaries.
#' @method summary JoiNMeDynPred
#' @export
summary.JoiNMeDynPred <- function(object, ...) {
    if (!inherits(object, "JoiNMeDynPred")) {
        cli::cli_abort(c(
            x = "Object must be a {.cls JoiNMeDynPred} prediction.",
            i = "Call predict() on a JoiNMeFit object first."
        ))
    }
    cached <- object$cache_get("summary")
    if (!is.null(cached)) return(cached)

    n_long_rows <- if (!is.null(object$predictions$longitudinal)) nrow(object$predictions$longitudinal) else 0
    n_surv_rows <- if (!is.null(object$predictions$survival)) nrow(object$predictions$survival) else 0
    n_cumhaz_rows <- if (!is.null(object$predictions$cumhaz)) nrow(object$predictions$cumhaz) else 0
    n_long_subjects <- if (!is.null(object$predictions$longitudinal)) length(unique(object$predictions$longitudinal$id)) else 0
    n_surv_subjects <- if (!is.null(object$predictions$survival)) length(unique(object$predictions$survival$id)) else 0
    n_cumhaz_subjects <- if (!is.null(object$predictions$cumhaz)) length(unique(object$predictions$cumhaz$id)) else 0

    overview_table <- data.frame(
        metric = c(
            "longitudinal_rows",
            "longitudinal_subjects",
            "survival_rows",
            "survival_subjects",
            "cumhaz_rows",
            "cumhaz_subjects",
            "prediction_scale",
            "posterior_draws"
        ),
        value = c(
            n_long_rows,
            n_long_subjects,
            n_surv_rows,
            n_surv_subjects,
            n_cumhaz_rows,
            n_cumhaz_subjects,
            if (!is.null(object$metadata$scales)) paste(object$metadata$scales, collapse = ",") else object$metadata$scale %||% NA_character_,
            object$n_samples %||% object$metadata$n_samples %||% NA_real_
        ),
        stringsAsFactors = FALSE
    )

    median_survival_table <- NULL
    if (!is.null(object$draws$survival) && length(object$draws$survival) > 0) {
        median_survival_rows <- lapply(names(object$draws$survival), function(subject_id) {
            subject_survival <- object$draws$survival[[subject_id]]
            draw_median_times <- .median_survival_time_by_draw(
                survival_draw_matrix = subject_survival$matrix,
                time_grid = subject_survival$time,
                threshold = 0.5
            )
            draw_summary <- .summarize_draw_vector_with_diagnostics(draw_median_times)
            cbind(
                id = subject_id,
                term = "median_survival_time",
                n_reached = sum(is.finite(draw_median_times)),
                n_total = length(draw_median_times),
                draw_summary,
                stringsAsFactors = FALSE
            )
        })
        median_survival_table <- do.call(rbind, median_survival_rows)
    }

    random_effects_id_table <- NULL
    if (!is.null(object$draws$random_effects_id) && length(object$draws$random_effects_id) > 0) {
        random_effect_rows <- lapply(names(object$draws$random_effects_id), function(subject_id) {
            subject_effects <- object$draws$random_effects_id[[subject_id]]
            u_id_matrix <- subject_effects$matrix
            u_id_terms <- subject_effects$terms %||% paste0("u_id[", seq_len(ncol(u_id_matrix)), "]")

            per_term <- lapply(seq_len(ncol(u_id_matrix)), function(term_index) {
                draw_summary <- .summarize_draw_vector_with_diagnostics(u_id_matrix[, term_index])
                cbind(
                    id = subject_id,
                    term = u_id_terms[term_index],
                    draw_summary,
                    stringsAsFactors = FALSE
                )
            })
            do.call(rbind, per_term)
        })
        random_effects_id_table <- do.call(rbind, random_effect_rows)
    }

    marker_corr_depends_on_id <- isTRUE(object$metadata$marker_corr_depends_on_id)
    any_re_indep <- any(as.integer(c(
        object$metadata$indep_id_re %||% 0L,
        object$metadata$indep_marker_re %||% 0L,
        object$metadata$indep_idmarker_cov %||% 0L
    )) == 1L)

    random_effects_marker_id_table <- NULL
    if (marker_corr_depends_on_id && !is.null(object$draws$random_effects_marker_id) && length(object$draws$random_effects_marker_id) > 0) {
        marker_id_rows <- lapply(names(object$draws$random_effects_marker_id), function(subject_id) {
            subject_effects <- object$draws$random_effects_marker_id[[subject_id]]
            w_idm_matrix <- subject_effects$matrix
            if (is.null(w_idm_matrix) || ncol(w_idm_matrix) == 0) return(NULL)

            marker_terms <- colnames(w_idm_matrix)
            if (is.null(marker_terms)) {
                marker_terms <- paste0("marker_id[", seq_len(ncol(w_idm_matrix)), "]")
            }

            per_term <- lapply(seq_len(ncol(w_idm_matrix)), function(term_index) {
                draw_summary <- .summarize_draw_vector_with_diagnostics(w_idm_matrix[, term_index])
                term_label <- as.character(marker_terms[term_index])
                marker_label <- if (grepl("::", term_label, fixed = TRUE)) {
                    sub("::.*$", "", term_label)
                } else {
                    NA_character_
                }
                basis_term <- if (grepl("::", term_label, fixed = TRUE)) {
                    sub("^.*::", "", term_label)
                } else {
                    term_label
                }
                cbind(
                    id = subject_id,
                    marker = marker_label,
                    term = basis_term,
                    draw_summary,
                    stringsAsFactors = FALSE
                )
            })
            do.call(rbind, per_term)
        })
        marker_id_rows <- Filter(Negate(is.null), marker_id_rows)
        if (length(marker_id_rows) > 0) {
            random_effects_marker_id_table <- do.call(rbind, marker_id_rows)
        }
    }

    corr_marker_id_table <- NULL
    if (marker_corr_depends_on_id && !is.null(object$draws$random_effects_marker_id) && length(object$draws$random_effects_marker_id) > 0) {
        corr_rows <- lapply(names(object$draws$random_effects_marker_id), function(subject_id) {
            subject_effects <- object$draws$random_effects_marker_id[[subject_id]]
            corr_draws <- subject_effects$corr
            if (is.null(corr_draws) || length(dim(corr_draws)) != 3) return(NULL)

            terms <- as.character(subject_effects$terms %||% paste0("w_idm[", seq_len(dim(corr_draws)[2]), "]"))
            q_dim <- dim(corr_draws)[2]

            per_cell <- lapply(seq_len(q_dim), function(r_idx) {
                lapply(seq_len(q_dim), function(c_idx) {
                    draw_summary <- .summarize_draw_vector_with_diagnostics(corr_draws[, r_idx, c_idx])
                    cbind(
                        id = subject_id,
                        row = terms[r_idx],
                        col = terms[c_idx],
                        draw_summary,
                        stringsAsFactors = FALSE
                    )
                })
            })
            do.call(rbind, unlist(per_cell, recursive = FALSE))
        })
        corr_rows <- Filter(Negate(is.null), corr_rows)
        if (length(corr_rows) > 0) {
            corr_marker_id_table <- do.call(rbind, corr_rows)
            if (isTRUE(any_re_indep)) {
                corr_marker_id_table <- corr_marker_id_table[corr_marker_id_table$row == corr_marker_id_table$col, , drop = FALSE]
            }
        }
    }

    diagnostics_table <- NULL
    bind_rows_safe <- function(...) {
        row_list <- Filter(Negate(is.null), list(...))
        if (length(row_list) == 0) return(NULL)
        all_columns <- unique(unlist(lapply(row_list, names), use.names = FALSE))
        row_list <- lapply(row_list, function(df) {
            missing_columns <- setdiff(all_columns, names(df))
            if (length(missing_columns) > 0) {
                for (column_name in missing_columns) df[[column_name]] <- NA
            }
            df[, all_columns, drop = FALSE]
        })
        do.call(rbind, row_list)
    }
    diagnostics_source <- bind_rows_safe(
        if (!is.null(random_effects_id_table)) transform(random_effects_id_table, section = "random_effects_id") else NULL,
        if (!is.null(random_effects_marker_id_table)) transform(random_effects_marker_id_table, section = "random_effects_marker_id") else NULL,
        if (!is.null(median_survival_table)) transform(median_survival_table, section = "median_survival_time") else NULL
    )
    term_diag <- .term_diagnostics_from_tables(diagnostics_source)
    sampler_diag <- object$metadata$sampler_diagnostics %||% list()

    diagnostics_table <- .build_common_diagnostics_table(
        draws = as.numeric(object$n_samples %||% object$metadata$n_samples %||% NA_real_),
        divergences = as.numeric(sampler_diag$divergences %||% 0),
        treedepth_hits = as.numeric(sampler_diag$treedepth_hits %||% 0),
        ebfmi_min = as.numeric(sampler_diag$ebfmi_min %||% NA_real_),
        max_rhat = as.numeric(term_diag$max_rhat %||% sampler_diag$max_rhat %||% NA_real_),
        min_ess_bulk = as.numeric(term_diag$min_ess_bulk %||% sampler_diag$min_ess_bulk %||% NA_real_),
        min_ess_tail = as.numeric(term_diag$min_ess_tail %||% sampler_diag$min_ess_tail %||% NA_real_),
        n_terms_total = as.numeric(term_diag$n_terms_total %||% 0),
        n_terms_bad_rhat = as.numeric(term_diag$n_terms_bad_rhat %||% 0),
        n_terms_low_ess_bulk = as.numeric(term_diag$n_terms_low_ess_bulk %||% 0),
        n_terms_low_ess_tail = as.numeric(term_diag$n_terms_low_ess_tail %||% 0)
    )

    tables <- list(
        diagnostics = diagnostics_table,
        overview = overview_table,
        median_survival_time = median_survival_table,
        random_effects_id = random_effects_id_table,
        random_effects_marker_id = random_effects_marker_id_table,
        corr_marker_id = corr_marker_id_table
    )

    summary_obj <- SummaryJoiNMeDynPred$new(tables = tables, metadata = object$metadata)
    object$cache_set("summary", summary_obj)
    summary_obj
}

#' Extract predicted random effects from dynamic predictions
#'
#' @description
#' Returns predicted random effects from a `JoiNMeDynPred` object. Marker-by-id
#' random effects are available only when marker covariance is configured to be
#' subject-dependent.
#'
#' @param object A `JoiNMeDynPred` object.
#' @param ... Unused.
#'
#' @return A named list containing random-effects summary tables.
#' @importFrom lme4 ranef
#' @export
ranef.JoiNMeDynPred <- function(object, ...) {
    assertthat::assert_that(inherits(object, "JoiNMeDynPred"), msg = "Object must be a JoiNMeDynPred instance.")

    if (!isTRUE(object$metadata$marker_corr_depends_on_id)) {
        cli::cli_abort(c(
            x = "Predicted marker-by-id random effects are only available when marker covariance depends on id.",
            i = "Refit with subject-dependent covariance structure in {.arg formulaVCov} and marker-by-id random effects (Q_idm > 0)."
        ))
    }

    sum_obj <- summary(object)
    list(
        formulaLong = list(
            marker_by_id = sum_obj$tables$random_effects_marker_id
        )
    )
}

#' Extract predicted covariance summaries
#'
#' @description
#' Returns per-subject covariance summaries for marker-by-id random effects from
#' a `JoiNMeDynPred` object. This method is available only when marker covariance
#' depends on id.
#'
#' @param object A `JoiNMeDynPred` object.
#' @param ... Unused.
#'
#' @return A named list containing covariance summary tables.
#' @method vcov JoiNMeDynPred
#' @export
vcov.JoiNMeDynPred <- function(object, ...) {
    assertthat::assert_that(inherits(object, "JoiNMeDynPred"), msg = "Object must be a JoiNMeDynPred instance.")

    if (!isTRUE(object$metadata$marker_corr_depends_on_id)) {
        cli::cli_abort(c(
            x = "Predicted marker-by-id covariance is only available when marker covariance depends on id.",
            i = "Refit with subject-dependent covariance structure in {.arg formulaVCov} and marker-by-id random effects (Q_idm > 0)."
        ))
    }

    sum_obj <- summary(object)
    list(
        formulaLong = list(
            marker_by_id = sum_obj$tables$corr_marker_id
        )
    )
}

#' Print dynamic prediction results
#'
#' @param x A dynamic prediction object.
#' @param ... Unused.
#'
#' @return Invisibly returns the object.
#' @export
print.JoiNMeDynPred <- function(x, ...) {
    assertthat::assert_that(inherits(x, "JoiNMeDynPred"), msg = "Object must be a JoiNMeDynPred instance.")

    cat("JoiNMe dynamic prediction\n")
    cat("=======================\n")
    if (!is.null(x$call)) {
        cat("Call:\n")
        print(x$call)
    }
    if (!is.null(x$metadata$pred_type)) {
        cat("Prediction type: ", x$metadata$pred_type, "\n", sep = "")
    }
    scales <- x$metadata$scales %||% x$metadata$scale
    if (!is.null(scales)) {
        if (length(scales) > 1L) {
            cat("Scales: ", paste(scales, collapse = ", "), "\n", sep = "")
        } else {
            cat("Scale: ", scales, "\n", sep = "")
        }
    }
    if (!is.null(x$metadata$n_subjects)) {
        cat("Subjects: ", x$metadata$n_subjects, "\n", sep = "")
    }
    if (!is.null(x$n_samples)) {
        cat("Posterior draws: ", x$n_samples, "\n", sep = "")
    }
    if (!is.null(x$tmax)) {
        cat("tmax: ", x$tmax, "\n", sep = "")
    }
    cat("Use summary() for prediction summaries.\n")
    cat("Use plot() for trajectory and interval visualisation.\n")
    invisible(x)
}

#' @export
print.summary_JoiNMeDynPred <- function(x, ...) {
    .cli_summary_heading("Prediction summary", level = 1L)
    meta_lines <- character(0)
    if (!is.null(x$metadata$pred_type)) {
        meta_lines <- c(meta_lines, paste0("Prediction type: ", x$metadata$pred_type))
    }
    if (!is.null(x$metadata$n_subjects)) {
        meta_lines <- c(meta_lines, paste0("Subjects: ", x$metadata$n_subjects))
    }
    .cli_print_bullets(meta_lines)
    .cli_print_table_section("Diagnostics", x$tables$diagnostics, level = 2L, formatter = .format_common_diagnostics_for_print)
    .cli_print_table_section("Overview", x$tables$overview, level = 2L)
    .cli_print_table_section("Median survival time", x$tables$median_survival_time, level = 2L)
    .cli_print_table_section("Predicted id random effects", x$tables$random_effects_id, level = 2L)
    .cli_print_table_section("Predicted marker-by-id random effects", x$tables$random_effects_marker_id, level = 2L)
    .cli_print_table_section("Predicted marker-by-id covariance", x$tables$corr_marker_id, level = 2L)
    invisible(x)
}

.summarize_pred_surv <- function(mat, t_grid, id) {
    # mat: n_samples x n_times_surv
    sum_fun <- function(x) c(Mean = mean(x), Median = median(x), SD = sd(x), L95 = quantile(x, 0.025, names = FALSE), U95 = quantile(x, 0.975, names = FALSE))
    stats <- t(apply(mat, 2, sum_fun))

    data.frame(id = id, time = t_grid, Survival = stats[, "Mean"], Median = stats[, "Median"], Est.Error = stats[, "SD"], L95 = stats[, "L95"], U95 = stats[, "U95"])
}

.summarize_pred_cumhaz <- function(mat, t_grid, id) {
    # mat: n_samples x n_times_surv
    sum_fun <- function(x) c(Mean = mean(x), Median = median(x), SD = sd(x), L95 = quantile(x, 0.025, names = FALSE), U95 = quantile(x, 0.975, names = FALSE))
    stats <- t(apply(mat, 2, sum_fun))

    data.frame(id = id, time = t_grid, Cumhaz = stats[, "Mean"], Median = stats[, "Median"], Est.Error = stats[, "SD"], L95 = stats[, "L95"], U95 = stats[, "U95"])
}

# Compute rich quantile information for plotting flexibility
.compute_quantiles_long <- function(mat, t_grid, marker_pred, levels, id, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)) {
    # mat: n_samples x n_obs_pred
    # Returns data frame with quantiles for each time point
    q_labels <- .quantile_colnames(probs)
    quantile_mat <- t(apply(mat, 2, quantile, probs = probs, names = FALSE))
    quant_df <- as.data.frame(quantile_mat)
    names(quant_df) <- q_labels

    if (length(marker_pred) == 0) {
        marker_pred <- rep(NA_integer_, length(t_grid))
    }
    if (is.null(levels) || length(levels) == 0) {
        levels <- as.character(sort(unique(marker_pred)))
    }

    data.frame(
        id = id,
        time = t_grid,
        marker_idx = marker_pred,
        marker = if (all(is.na(marker_pred))) rep(NA_character_, length(t_grid)) else levels[marker_pred],
        quant_df,
        mean = colMeans(mat),
        sd = apply(mat, 2, sd)
    )
}

.compute_quantiles_surv <- function(mat, t_grid, id, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)) {
    # mat: n_samples x n_times_surv
    q_labels <- .quantile_colnames(probs)
    quantile_mat <- t(apply(mat, 2, quantile, probs = probs, names = FALSE))
    quant_df <- as.data.frame(quantile_mat)
    names(quant_df) <- q_labels

    data.frame(
        id = id,
        time = t_grid,
        quant_df,
        mean = colMeans(mat),
        sd = apply(mat, 2, sd)
    )
}

.quantile_colnames <- function(probs) {
    pct <- 100 * probs
    lbl <- formatC(pct, format = "fg", digits = 6)
    lbl <- gsub("\\s+", "", lbl)
    lbl <- sub("\\.?0+$", "", lbl)
    paste0("q", lbl)
}

.quantile_probs_from_ci <- function(ci_levels) {
    bounds <- unlist(lapply(ci_levels, function(level) c((1 - level) / 2, (1 + level) / 2)))
    sort(unique(c(0.5, bounds)))
}

.validate_ci_levels <- function(ci_levels) {
    if (is.null(ci_levels) || length(ci_levels) == 0) {
        cli::cli_abort(c(
            x = "{.arg ci_levels} must be a non-empty numeric vector.",
            i = "Example: ci_levels = c(0.5, 0.95)."
        ))
    }
    if (!is.numeric(ci_levels)) {
        cli::cli_abort(c(
            x = "{.arg ci_levels} must be numeric.",
            i = "Provide values between 0 and 1."
        ))
    }
    if (any(ci_levels <= 0 | ci_levels >= 1)) {
        cli::cli_abort(c(
            x = "{.arg ci_levels} must be strictly between 0 and 1.",
            i = "Example: ci_levels = c(0.5, 0.8, 0.95)."
        ))
    }
    sort(unique(ci_levels))
}

.resolve_time_grid <- function(times, id, t_cond, tmax_val, time_horizon, default_n = 50, min_points = 50, kind = "longitudinal") {
    grid <- NULL
    horizon_is_default <- is.finite(tmax_val) && isTRUE(all.equal(as.numeric(time_horizon), as.numeric(tmax_val)))
    raw_upper_time <- if (horizon_is_default) tmax_val else (t_cond + time_horizon)
    upper_support <- if (is.finite(tmax_val)) max(t_cond, tmax_val) else raw_upper_time
    upper_time <- min(raw_upper_time, upper_support)

    if (raw_upper_time > upper_support) {
        cli::cli_warn(c(
            x = "Requested {kind} horizon exceeds training time support for subject {id}.",
            i = "Truncating to times <= {format(signif(upper_support, 6), scientific = FALSE)} to avoid extrapolation artifacts."
        ))
    }

    if (is.null(times)) {
            if (default_n < min_points) {
                cli::cli_warn(c(
                    x = "Requested {kind} grid size {default_n} is below the minimum of {min_points}.",
                    i = "Set {.arg n_times} = {min_points} or ensure {.arg times} has at least {min_points} time points."
                ))
            }
            grid <- seq(t_cond, upper_time, length.out = max(default_n, min_points))
    } else if (is.list(times)) {
        key <- as.character(id)
        grid <- times[[key]]
        if (is.null(grid)) {
            cli::cli_warn(c(
                x = "No {kind} time grid found for subject {id}.",
                i = "Falling back to the default time grid."
            ))
            grid <- seq(t_cond, upper_time, length.out = max(default_n, min_points))
        }
    } else {
        grid <- times
    }

    if (any(grid > upper_time, na.rm = TRUE)) {
        cli::cli_warn(c(
            x = "{kind} grid exceeds configured horizon/support for subject {id}.",
            i = "Truncating to times <= {format(signif(upper_time, 6), scientific = FALSE)}."
        ))
    }

    grid <- grid[grid >= t_cond]
    grid <- grid[grid <= upper_time]
    grid <- unique(sort(grid))

    if (length(grid) == 0) {
        grid <- seq(t_cond, upper_time, length.out = max(default_n, min_points))
    }

    if (length(grid) < min_points && default_n == min_points) {
        cli::cli_warn(c(
            x = "{kind} prediction grid has fewer than {min_points} time points for subject {id}.",
            i = "Consider supplying a denser time grid for smoother intervals."
        ))
    }
    grid
}

.cmdstan_threads_enabled <- function(mod) {
    if (is.null(mod)) return(FALSE)
    cpp_opts <- tryCatch(mod$cpp_options, error = function(e) NULL)
    if (is.function(cpp_opts)) {
        cpp_opts <- tryCatch(cpp_opts(), error = function(e) NULL)
    }
    isTRUE(cpp_opts$stan_threads)
}

`%||%` <- function(x, y) if (is.null(x)) y else x
