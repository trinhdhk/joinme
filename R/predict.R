#' Dynamic Prediction for Joint Models
#'
#' @importFrom stats predict median sd quantile na.omit optim terms
#' @importFrom utils head modifyList
#' @importFrom dplyr %>% all_of
#' @keywords internal
#' @name predict.JoinMeFit
NULL

#' Dynamic prediction for a fitted `JoinMeFit` model
#'
#' @rdname predict.JoinMeFit
#' @description
#' Performs dynamic prediction for a new subject using a fitted `JoinMeFit` joint model.
#' It calculates the posterior predictive distribution of longitudinal trajectories
#' and survival probabilities conditional on the subject's observed history up to a specific time point.
#'
#' @param object A fitted object of class `JoinMeFit`.
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
#' @param scale Character. Longitudinal prediction scale:
#'   - "epred": expected response (inverse-link),
#'   - "linpred": linear predictor,
#'   - "predict": predictive draw (includes noise).
#' @param times Numeric vector of times at which to predict the longitudinal and survival trajectories.
#' Can also be a named list of numeric vectors (one per subject id). If NULL, a grid
#' from `Tstart` to `Tstart + 5` (scaled) is generated. If fewer than 50 points are
#' supplied for a subject, a warning is emitted when `n_times = 50`.
#' @param Tstart Numeric scalar. The conditioning time (last observation time).
#' If NULL, defaults to the maximum observed time in `newdataLong` for the subject.
#' @param ci_levels Numeric vector of credible interval levels for plotting.
#' Must be strictly between 0 and 1.
#' @param tmax Numeric scalar. The maximum time used for scaling during model fitting.
#' **CRITICAL**: If this is not provided and cannot be inferred from the object, predictions will be on the wrong time scale.
#' @param n_samples Integer. Number of posterior draws to use for prediction (Monte Carlo integration). Default 200.
#' @param n_times Integer. Number of time points for the prediction grid when `times` is NULL. Default 50.
#' @param control Named list for prediction configuration.
#'   - cmdstanr::model$sample() arguments (e.g., `chains`, `parallel_chains`,
#'     `iter_warmup`, `iter_sampling`, `seed`, `refresh`, `adapt_delta`).
#'     Values override defaults except `data`.
#'   - engine: "cmdstanr" or "rstan". Defaults to options(stan_preferred_engine).
#'   - threads_per_chain: integer; if > 1 uses `joinme_dynpred_threading.stan`.
#'     For engine = "rstan", threading uses options(stan.thread = threads_per_chain).
#'   - grainsize: integer; reduce_sum grainsize for threaded prediction
#'     (default max(1, ceiling(n_id/(4*threads_per_chain)))).
#' @param seed Integer. Random seed for reproducibility of random effect sampling.
#' @importFrom stats predict median sd quantile na.omit optim terms
#' @param ... Additional arguments (unused).
#' 
#' @return A list with two components:
#' \item{longitudinal}{A data.frame containing longitudinal predictions (Mean, Median, SD, 95% CrI) for each time point in `times`. Columns: id, time, marker, Estimate, Median, Est.Error, L95, U95.}
#' \item{longitudinal_fitted}{A data.frame containing fitted values for observed history on the linpred and epred scales. Includes a `scale` column.}
#' \item{survival}{A data.frame containing survival probabilities (S(t|Tstart)) for each time point in `times`. Columns: id, time, Survival, Median, Est.Error, L95, U95.}
#' \item{cumhaz}{A data.frame containing conditional cumulative hazards H(t|Tstart). Columns: id, time, Cumhaz, Median, Est.Error, L95, U95.}
#' \item{draws}{A list containing raw posterior draws (`longitudinal`, `survival`, `cumhaz`).}
#'
#' @details
#' The function uses a Bayesian approach. For each of the `n_samples` posterior draws from the fitted model:
#' 1.  It samples subject-specific random effects (ID-level and marker-level) from their posterior distribution *conditional* on the observed history in `newdataLong` and the fixed parameters.
#' 2.  It calculates the expected longitudinal value and survival probability at the requested future times.
#' 3.  Results are pooled to provide the marginal predictive distribution.
#'
#' @export
predict.JoinMeFit <- function(object,
                           newdataLong,
                           newdataEvent,
                           process = c("longitudinal", "event"),
                           pred_type = c("per_marker_id", "marginal_marker", "marginal_id", "marker_subject", "subject_marker"),
                           scale = c("epred", "linpred", "predict"),
                           times = NULL,
                           Tstart = NULL,
                           tmax = NULL,
                           ci_levels = c(0.5, 0.95),
                           n_samples = 200,
                           n_times = 50,
                           control = list(),
                           seed = 123,
                           ...) {
    if (!inherits(object, "JoinMeFit")) {
        cli::cli_abort(c(
            x = "Object must be a {.cls JoinMeFit} fit.",
            i = "Fit the model with joinme() before predicting."
        ))
    }
    process <- match.arg(process, several.ok = TRUE)
    pred_type <- match.arg(pred_type)
    pred_type <- switch(pred_type,
                        marker_subject = "marginal_marker",
                        subject_marker = "marginal_id",
                        pred_type)
    scale <- match.arg(scale)

    ci_levels <- .validate_ci_levels(ci_levels)
    quantile_probs <- .quantile_probs_from_ci(ci_levels)

    # 1. Recover Metadata
    meta <- .recover_metadata(object, tmax)
    tmax_val <- meta$tmax
    knots <- meta$knots
    col_means <- meta$col_means

    # 2. Parse Formulas
    forms <- .parse_formulas(object)

    # 3. Identify Subjects
    id_var <- eval(object$call$id_var) %||% "id"
    ids <- unique(newdataEvent[[id_var]])
    if (length(ids) == 0) {
        cli::cli_abort("No subjects found in {.code newdataEvent}.")
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

    stan_candidates <- if (threads_per_chain > 1) {
        c(
            system.file("stan/joinme_dynpred_threading.stan", package = "joinme"),
            file.path("inst", "stan", "joinme_dynpred_threading.stan"),
            file.path("..", "inst", "stan", "joinme_dynpred_threading.stan"),
            file.path("..", "..", "inst", "stan", "joinme_dynpred_threading.stan")
        )
    } else {
        c(
            system.file("stan/joinme_dynpred.stan", package = "joinme"),
            file.path("inst", "stan", "joinme_dynpred.stan"),
            file.path("..", "inst", "stan", "joinme_dynpred.stan"),
            file.path("..", "..", "inst", "stan", "joinme_dynpred.stan")
        )
    }
    stan_file <- stan_candidates[file.exists(stan_candidates)][1]
    if (is.na(stan_file) || !nzchar(stan_file)) {
        cli::cli_abort(c(
            x = "Prediction Stan model file not found.",
            i = "Reinstall the package or restore inst/stan files."
        ))
    }

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

    cpp_opts <- if (threads_per_chain > 1) list(stan_threads = TRUE) else NULL
    force_recompile <- isTRUE(control$force_recompile %||% getOption("joinme.force_recompile", FALSE))
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

    allowed_data <- if (engine == "rstan" && isTRUE(getOption("joinme.filter_rstan_data", FALSE))) {
        .stan_data_names(stan_file)
    } else {
        NULL
    }

    # 5. Extract Draws
    draws_list <- .extract_draws_for_pred(object, n_samples, seed)
    n_samples_extracted <- if (ncol(draws_list$alpha_vcov_reg) > 0) {
        nrow(draws_list$alpha_vcov_reg)
    } else {
        nrow(draws_list$beta_fixed)
    }

    n_cores <- parallel::detectCores(logical = FALSE) %||% 1L
    max_threads <- min(n_samples_extracted, n_cores)
    if (threads_per_chain > max_threads) {
        cli::cli_warn(c(
            x = "Requested {threads_per_chain} threads exceeds max {max_threads}.",
            i = "Capping threads_per_chain to {max_threads}."
        ))
        threads_per_chain <- max_threads
    }

    grainsize <- control$grainsize
    if (is.null(grainsize)) {
        n_id_pred <- length(ids)
        n_chains <- control$chains %||% 1L
        denom <- 4L * as.integer(threads_per_chain) * as.integer(n_chains)
        denom <- max(1L, denom)
        grainsize <- max(1L, as.integer(ceiling(n_id_pred / denom)))
        grainsize <- min(as.integer(n_cores), grainsize)
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
    sample_control <- control[setdiff(names(control), c(
        "threads_per_chain",
        "threads",
        "mc.cores",
        "grainsize"
    ))]
    results_long <- list()
    results_surv <- list()
    results_cumhaz <- list()
    results_long_fit <- list()
    draws_long_list <- list()
    draws_surv_list <- list()
    draws_cumhaz_list <- list()
    draws_long_fit <- list()
    quantiles_long_list <- list()
    quantiles_surv_list <- list()
    quantiles_cumhaz_list <- list()
    quantiles_long_fit <- list()

    time_var <- eval(object$call$time_var) %||% "time"
    marker_var <- eval(object$call$marker_var) %||% "marker"

    for (id in ids) {
        dE <- newdataEvent[newdataEvent[[id_var]] == id, , drop = FALSE]
        dL <- newdataLong[newdataLong[[id_var]] == id, , drop = FALSE]

        if (nrow(dL) == 0) {
            cli::cli_warn("No longitudinal data for ID {id}.")
            next
        }

        # T_cond
        if (!is.null(Tstart)) {
            t_cond <- Tstart
        } else {
            t_cond <- max(dL[[time_var]], na.rm = TRUE)
        }

        # Grids
        t_grid <- numeric(0)
        if ("longitudinal" %in% process) {
            t_grid <- .resolve_time_grid(
                times,
                id,
                t_cond,
                tmax_val,
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
                default_n = n_times,
                min_points = 50,
                kind = "survival"
            )
        }

        # SD Prep
        grainsize_data <- if (threads_per_chain > 1) grainsize else NULL
        sd_pred <- .prepare_subject_standata(
            dE, dL, object, tmax_val, knots, col_means,
            t_cond, t_grid, t_surv_grid,
            forms, draws_list, grainsize_data
        )
        # Ensure time index arrays are preserved for cmdstanr JSON (avoid auto-unbox)
        sd_pred <- .coerce_rstan_time_indices(sd_pred)
        sd_pred <- .coerce_rstan_dist_arrays(sd_pred)
        sd_pred <- .coerce_rstan_vectors(sd_pred, c(
            "vec_cov_vcov",
            "const_data_cv",
            "const_data_cs",
            "const_data_vcov"
        ))
        if (!is.null(allowed_data) && length(allowed_data) > 0) {
            missing <- setdiff(allowed_data, names(sd_pred))
            if (length(missing) > 0) {
                cli::cli_abort(c(
                    x = "Prediction Stan data missing required fields for rstan: {paste(missing, collapse = ', ')}.",
                    i = "Check prediction standata construction and Stan data block."
                ))
            }
            sd_pred <- sd_pred[names(sd_pred) %in% allowed_data]
        }

        # Run Stan
        sample_args <- list(
            chains = 1,
            iter_warmup = 100,
            iter_sampling = 1, # Single pass over n_samples array
            fixed_param = FALSE,
            refresh = 0,
            show_messages = FALSE,
            seed = seed,
            adapt_delta = 0.95
        )
        sample_args <- utils::modifyList(sample_args, sample_control)
        if (threads_per_chain > 1) sample_args$threads_per_chain <- threads_per_chain
        sample_args$data <- sd_pred
        sample_args <- sample_args[!vapply(sample_args, is.null, logical(1))]
        if (engine == "cmdstanr") {
            allowed <- names(formals(mod$sample))
            sample_args <- sample_args[names(sample_args) %in% allowed]
            fit_pred <- do.call(mod$sample, sample_args)
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
        }

        long_var <- switch(scale,
               linpred = "y_pred_linpred",
               epred = "y_pred_epred",
               predict = "y_pred",
               "y_pred_epred")
        draws_mat <- suppressMessages(
            .get_draws_matrix(
                fit_pred,
                variables = c(long_var, "surv_prob", "cumhaz_cond", "y_fit_linpred", "y_fit_epred")
            )
        )

        if (nrow(dL) > 0) {
            marker_levels <- object$stan_data$marker_levels %||% levels(dL[[marker_var]]) %||% sort(unique(dL[[marker_var]]))
            marker_int_obs <- match(as.character(dL[[marker_var]]), marker_levels)
            if (any(is.na(marker_int_obs))) {
                cli::cli_abort(c(
                    x = "Markers in newdataLong don't match fitted model levels.",
                    i = "Ensure marker levels match the fitted model."
                ))
            }

            fit_lin_mat <- .extract_matrix_from_stan(draws_mat, "y_fit_linpred", nrow(dL), n_samples_extracted)
            fit_epred_mat <- .extract_matrix_from_stan(draws_mat, "y_fit_epred", nrow(dL), n_samples_extracted)

            draws_long_fit[[as.character(id)]] <- list(
                linpred = fit_lin_mat,
                epred = fit_epred_mat,
                time = dL[[time_var]],
                marker_idx = marker_int_obs
            )

            fit_lin_sum <- .summarize_pred_long(fit_lin_mat, dL[[time_var]], marker_int_obs, marker_levels, id)
            fit_lin_sum$scale <- "linpred"
            fit_epred_sum <- .summarize_pred_long(fit_epred_mat, dL[[time_var]], marker_int_obs, marker_levels, id)
            fit_epred_sum$scale <- "epred"
            results_long_fit[[as.character(id)]] <- rbind(fit_lin_sum, fit_epred_sum)

            fit_lin_quant <- .compute_quantiles_long(
                fit_lin_mat, dL[[time_var]], marker_int_obs, marker_levels, id, probs = quantile_probs
            )
            fit_lin_quant$scale <- "linpred"
            fit_epred_quant <- .compute_quantiles_long(
                fit_epred_mat, dL[[time_var]], marker_int_obs, marker_levels, id, probs = quantile_probs
            )
            fit_epred_quant$scale <- "epred"
            quantiles_long_fit[[as.character(id)]] <- rbind(fit_lin_quant, fit_epred_quant)
        }

        if (length(t_grid) > 0) {
            # Extract predictions for all (time, marker) combinations
            n_obs_pred_total <- nrow(sd_pred$mat_fixed_pred)  # This is n_time_points * n_markers
            n_markers_actual <- length(unique(dL[[marker_var]]))
            long_mat <- .extract_matrix_from_stan(draws_mat, long_var, n_obs_pred_total, n_samples_extracted)

            draws_long_list[[as.character(id)]] <- list(
                matrix = long_mat,
                marker_idx = sd_pred$idx_marker_pred,
                time = rep(t_grid, n_markers_actual)  # Repeat time grid for each marker
            )

            marker_levels <- object$stan_data$marker_levels %||% levels(dL[[marker_var]]) %||% sort(unique(dL[[marker_var]]))
            results_long[[as.character(id)]] <- .summarize_pred_long(long_mat, rep(t_grid, n_markers_actual), sd_pred$idx_marker_pred, marker_levels, id)
            
            # Also compute rich quantiles for plotting
            quantiles_long_list[[as.character(id)]] <- .compute_quantiles_long(
                long_mat, rep(t_grid, n_markers_actual), sd_pred$idx_marker_pred,
                marker_levels, id, probs = quantile_probs
            )
        }
        if (length(t_surv_grid) > 0) {
            surv_mat <- .extract_matrix_from_stan(draws_mat, "surv_prob", length(t_surv_grid), n_samples_extracted)
            cumhaz_mat <- .extract_matrix_from_stan(draws_mat, "cumhaz_cond", length(t_surv_grid), n_samples_extracted)

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
    }

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
                # Average across subjects for each marker-time combination
                results_long_aggregated <- results_long_combined |>
                    dplyr::group_by(time, marker) |>
                    dplyr::summarise(
                        Estimate = mean(Estimate, na.rm = TRUE),
                        Median = mean(Median, na.rm = TRUE),
                        Est.Error = mean(Est.Error, na.rm = TRUE),
                        L95 = mean(L95, na.rm = TRUE),
                        U95 = mean(U95, na.rm = TRUE),
                        .groups = "drop"
                    ) |>
                    dplyr::mutate(id = "Marginal_Marker", .before = time)
            } else if (pred_type == "marginal_id") {
                # Average across markers for each subject-time combination
                results_long_aggregated <- results_long_combined |>
                    dplyr::group_by(id, time) |>
                    dplyr::summarise(
                        Estimate = mean(Estimate, na.rm = TRUE),
                        Median = mean(Median, na.rm = TRUE),
                        Est.Error = mean(Est.Error, na.rm = TRUE),
                        L95 = mean(L95, na.rm = TRUE),
                        U95 = mean(U95, na.rm = TRUE),
                        .groups = "drop"
                    )
            }
            # else per_marker_id: keep as is (results_long will be used)
        }
        
        # Group by time and marker, average quantiles across subjects
        quantiles_long_marker_pop <- all_marker_preds |>
            dplyr::group_by(time, marker) |>
            dplyr::summarise(
                dplyr::across(dplyr::all_of(q_cols), ~ mean(.x, na.rm = TRUE)),
                mean = mean(mean, na.rm = TRUE),
                sd = mean(sd, na.rm = TRUE),
                .groups = "drop"
            ) |>
            dplyr::mutate(id = "Population_Marker", marker_idx = NA_integer_)
        
        # Overall population trajectory (average across subjects and markers)
        quantiles_long_overall_pop <- all_marker_preds |>
            dplyr::group_by(time) |>
            dplyr::summarise(
                dplyr::across(dplyr::all_of(q_cols), ~ mean(.x, na.rm = TRUE)),
                mean = mean(mean, na.rm = TRUE),
                sd = mean(sd, na.rm = TRUE),
                .groups = "drop"
            ) |>
            dplyr::mutate(id = "Population_Overall", marker = "All", marker_idx = NA_integer_)
    }

    # Prepare metadata
    # Extract variable names from forms for use in plotting
    formula_long <- object$call$formulaLong %||% forms$formulaLong
    resp_var <- tryCatch(
        all.vars(formula_long)[1],
        error = function(e) NA_character_
    )
    
    metadata <- list(
        conditioning_time = if (!is.null(Tstart)) Tstart else max(newdataLong[[eval(object$call$time_var) %||% "time"]], na.rm = TRUE),
        n_samples = n_samples_extracted,
        n_subjects = length(ids),
        pred_type = pred_type,
        scale = scale,
        fitted_scales = c("linpred", "epred"),
        ci_levels = ci_levels,
        quantile_levels = quantile_probs,
        time_max = tmax_val,
        formula_long = formula_long,
        response_var = resp_var,
        id_var = eval(object$call$id_var) %||% "id",
        time_var = eval(object$call$time_var) %||% "time",
        marker_var = eval(object$call$marker_var) %||% "marker"
    )

    JoinMeDynPred$new(
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
        n_samples = n_samples_extracted
    )
}

#' Posterior Linear Predictor
#'
#' @rdname predict.JoinMeFit
#' @description
#' Convenience wrapper around the [predict] method returning the posterior
#' linear predictor (linpred scale).
#'
#' @param object A fitted object of class `JoinMeFit`.
#' @param ... Additional arguments passed to [predict].
#' 
#' @importFrom rstantools posterior_linpred
#' @return A `JoinMeDynPred` object with `metadata$scale = "linpred"`.
#'
#' @export
posterior_linpred.JoinMeFit <- function(object, ...) {
    predict(object, scale = "linpred", ...)
}

#' Posterior Expected Predictor
#'
#' @rdname predict.JoinMeFit
#' @description
#' Convenience wrapper around the [predict] method returning the posterior
#' expected predictor (epred scale).
#'
#' @param object A fitted object of class `JoinMeFit`.
#' @param ... Additional arguments passed to [predict].
#'
#' @importFrom rstantools posterior_linpred
#' @return A `JoinMeDynPred` object with `metadata$scale = "epred"`
#' @export
posterior_epred.JoinMeFit <- function(object, ...) {
    predict(object, scale = "epred", ...)
}

#' Posterior Predictive Draws
#' 
#' @rdname predict.JoinMeFit
#' @description
#' Convenience wrapper around the [predict] method returning posterior
#' predictive draws (includes observation noise).
#'
#' @param object A fitted object of class `JoinMeFit`.
#' @param ... Additional arguments passed to [predict].
#'
#' @importFrom rstantools posterior_predict
#' @return A `JoinMeDynPred` object with `metadata$scale = "predict"`.
#' @export
posterior_predict.JoinMeFit <- function(object, ...) {
    predict(object, scale = "predict", ...)
}

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

.extract_matrix_from_stan <- function(draws_mat, var_name, N_cols, N_rows) {
    # Stan output flattened: var_name[k, n]
    # We want matrix (N_rows x N_cols)
    mat <- matrix(0, N_rows, N_cols)
    for (n in 1:N_cols) {
        # Construct names like "y_pred_new[1,1]", "y_pred_new[2,1]" ...
        # But wait, CmdStanR output format depends.
        # If we requested matrix format from fit$draws(..., format="matrix")
        # Columns are variable names.
        # "y_pred_new[1,1]", "y_pred_new[2,1]" is likely provided.
        # Actually, simpler:
        # loop k=1..N_rows? No that's slow.
        # We need ALL rows [1..N_rows, n]
        # Regex or paste loop?
        # Paste is cleaner if we assume 1-based indexing for k.
        nms <- paste0(var_name, "[", 1:N_rows, ",", n, "]")
        if (all(nms %in% colnames(draws_mat))) {
            mat[, n] <- as.numeric(draws_mat[1, nms])
        }
    }
    mat
}

.recover_metadata <- function(object, tmax_arg) {
    tmax <- tmax_arg
    if (is.null(tmax)) {
        if (!is.null(object$tmax)) tmax <- object$tmax
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
                splines::bs(S_event, knots = knots, Boundary.knots = c(0, 1), degree = degree, intercept = FALSE)
            )
        } else if (basehaz_type == "ns") {
            B_raw <- suppressWarnings(
                splines::ns(S_event, knots = knots, Boundary.knots = c(0, 1), intercept = FALSE)
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

    list(tmax = tmax, knots = knots, col_means = col_means)
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
        formulaVcov = object$formulaVcov,
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
    if (isTRUE(sd$assoc_vcov) && sd$Q_idm < 1) {
        cli::cli_abort(c(
            x = "Association {.arg vcov} requires marker-by-id random effects (Q_idm > 0).",
            i = "Refit with an inner ( ... | id ) term inside the marker block, or drop {.arg vcov} from {.arg assoc}."
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

    if (sd$Q_idm > 0) {
        tau_marker_id <- get_mat(paste0("tau_w[", 1:sd$Q_idm, "]"))
        Lcorr_marker_id <- array(0, dim = c(n, sd$Q_idm, sd$Q_idm))
        for (r in 1:sd$Q_idm) {
            for (c in 1:sd$Q_idm) {
                nm <- paste0("Lcorr_w[", r, ",", c, "]")
                if (nm %in% colnames(dmat)) Lcorr_marker_id[, r, c] <- dmat[, nm]
            }
        }
    } else {
        tau_marker_id <- matrix(0, n, 0)
        Lcorr_marker_id <- array(0, dim = c(n, 0, 0))
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

    log_h0_intercept <- matrix(0, n, sd$K_event %||% 1L)
    for (k_ev in 1:(sd$K_event %||% 1L)) {
        nm <- paste0("log_h0_intercept[", k_ev, "]")
        if (nm %in% colnames(dmat)) log_h0_intercept[, k_ev] <- dmat[, nm]
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

    marker_weights_draws <- NULL
    if (all(paste0("marker_weights_eff[", 1:sd$D, "]") %in% colnames(dmat))) {
        marker_weights_draws <- get_mat(paste0("marker_weights_eff[", 1:sd$D, "]"))
    } else if (all(paste0("marker_weights[", 1:sd$D, "]") %in% colnames(dmat))) {
        marker_weights_draws <- get_mat(paste0("marker_weights[", 1:sd$D, "]"))
    } else if (!is.null(sd$marker_weights)) {
        marker_weights_draws <- matrix(rep(as.numeric(sd$marker_weights), each = n), nrow = n, byrow = TRUE)
    }

    list(
        n_samples = n,
        beta_fixed = beta_fixed,
        tau_id = tau_id, Lcorr_id = Lcorr_id,
        tau_marker = tau_marker, Lcorr_marker = Lcorr_marker,
        tau_marker_id = tau_marker_id, Lcorr_marker_id = Lcorr_marker_id,
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
        tau_vcov_reg = get_col("tau_L"),
        lambda_vcov_reg = if (M_cov > 0) get_mat(paste0("lambda_L[", 1:M_cov, "]")) else matrix(0, n, 0),
        log_h0_intercept = log_h0_intercept,
        bs_gamma_c = bs_gamma_c,
        gamma_hazard = gamma_hazard,
        marker_weights_draws = marker_weights_draws,
        sigma_y_shared = get_col("sigma_y", default = 0),
        sigma_marker_specific = if (sd$D > 0 && "sigma_marker[1]" %in% colnames(dmat)) {
            get_mat(paste0("sigma_marker[", 1:sd$D, "]"))
        } else {
            matrix(0, n, 0)
        },
        nu_marker = if (sd$D > 0 && "nu_marker[1]" %in% colnames(dmat)) {
            get_mat(paste0("nu_marker[", 1:sd$D, "]"))
        } else {
            matrix(0, n, 0)
        },
        phi_nb_marker = if (sd$D > 0 && "phi_nb_marker[1]" %in% colnames(dmat)) {
            get_mat(paste0("phi_nb_marker[", 1:sd$D, "]"))
        } else {
            matrix(0, n, 0)
        },
        alpha_skew_marker = if (sd$D > 0 && "alpha_skew_marker[1]" %in% colnames(dmat)) {
            get_mat(paste0("alpha_skew_marker[", 1:sd$D, "]"))
        } else {
            matrix(0, n, 0)
        },
        phi_beta_marker = if (sd$D > 0 && "phi_beta_marker[1]" %in% colnames(dmat)) {
            get_mat(paste0("phi_beta_marker[", 1:sd$D, "]"))
        } else {
            matrix(0, n, 0)
        },
        tau_sde_marker = if (sd$D > 0 && "tau_sde_marker[1]" %in% colnames(dmat)) {
            get_mat(paste0("tau_sde_marker[", 1:sd$D, "]"))
        } else {
            matrix(0, n, 0)
        },
        coeff_assoc_cv_total = get_col("alpha_cv_total"),
        coeff_assoc_cs_total = get_col("alpha_cs_total"),
        coeff_assoc_cv_mean = get_col("alpha_cv_mean"),
        coeff_assoc_cs_mean = get_col("alpha_cs_mean"),
        coeff_assoc_cv_marker = get_col("alpha_cv_marker"),
        coeff_assoc_cs_marker = get_col("alpha_cs_marker"),
        coeff_assoc_vcov_var = if (sd$Q_idm > 0) get_mat(paste0("alpha_vcov_var[", 1:sd$Q_idm, "]")) else matrix(0, n, 0),
        cutpoints_ord = if (!is.null(sd$K_ord) && sd$K_ord > 1 && "cutpoints_ord[1]" %in% colnames(dmat)) {
            get_mat(paste0("cutpoints_ord[", 1:(sd$K_ord - 1), "]"))
        } else {
            matrix(0, n, max(1, (sd$K_ord %||% 2L) - 1))
        }
    )
}

.prepare_subject_standata <- function(dE, dL, object, tmax, knots, col_means, t_cond, t_grid, t_surv_grid, forms, draws_list, grainsize = NULL) {
    id_var <- eval(object$call$id_var) %||% "id"
    time_var <- eval(object$call$time_var) %||% "time"
    marker_var <- eval(object$call$marker_var) %||% "marker"

    sd <- object$stan_data

    dL[[time_var]] <- dL[[time_var]] / tmax
    T_cond_scaled <- t_cond / tmax

    f_exp <- reformulas::expandDoubleVerts(forms$formulaLong)
    bars <- reformulas::findbars(f_exp)
    f_fix <- reformulas::nobars(f_exp)
    fixed_rhs <- stats::update(f_fix, . ~ .)
    fixed_rhs[[2]] <- NULL

    grp <- vapply(bars, function(b) as.character(b[[3]]), character(1))
    id_rhs_list <- .bar_terms_to_rhs_list(bars[which(grp == id_var)])
    nested <- .extract_nested_marker_terms(forms$formulaLong, marker_var, id_var)

    mat_fixed_obs <- .mm(fixed_rhs, dL)
    mat_id_obs <- if (length(id_rhs_list) > 0) {
        do.call(cbind, lapply(id_rhs_list, function(rhs) .mm(rhs, dL)))
    } else {
        matrix(0, nrow(dL), sd$R_id)
    }
    mat_marker_obs <- if (length(nested$mk_rhs_list) > 0) {
        do.call(cbind, lapply(nested$mk_rhs_list, function(rhs) .mm(rhs, dL)))
    } else {
        matrix(0, nrow(dL), sd$R_mk)
    }
    mat_marker_id_obs <- if (length(nested$idm_rhs_list) > 0) {
        do.call(cbind, lapply(nested$idm_rhs_list, function(rhs) .mm(rhs, dL)))
    } else {
        matrix(0, nrow(dL), sd$Q_idm)
    }

    dist_formulas <- .normalize_formula_dist(forms$formulaDist)
    dist_sigma_obs <- .build_dist_matrix(dist_formulas$sigma, dL)
    dist_nu_obs <- .build_dist_matrix(dist_formulas$nu, dL)
    dist_phi_obs <- .build_dist_matrix(dist_formulas$phi, dL)
    dist_alpha_obs <- .build_dist_matrix(dist_formulas$alpha, dL)
    dist_phi_beta_obs <- .build_dist_matrix(dist_formulas$phi_beta, dL)
    dist_tau_sde_obs <- .build_dist_matrix(dist_formulas$tau_sde, dL)

    fe_rhs <- stats::update(forms$formulaEvent, . ~ .)
    fe_rhs[[2]] <- NULL
    vec_cov_hazard <- .mm(fe_rhs, dE)

    fv_rhs <- stats::update(forms$formulaVcov, . ~ .)
    fv_rhs[[2]] <- NULL
    vec_cov_vcov <- .mm(fv_rhs, dE)
    if (ncol(vec_cov_vcov) == 1 && colnames(vec_cov_vcov)[1] == "(Intercept)" && sd$K_cov == 1) vec_cov_vcov <- matrix(0, 1, 1)

    gk <- .gk15_nodes()
    u_cond <- 0.5 * T_cond_scaled * (gk + 1)
    u_cond_fwd <- u_cond + sd$eps_fd

    .eval_on_times <- function(rhs_l, times) {
        if (length(rhs_l) == 0) {
            return(matrix(0, 15, 0))
        }
        dd <- dE[rep(1, 15), , drop = FALSE]
        dd[[time_var]] <- times
        do.call(cbind, lapply(rhs_l, function(rhs) .mm(rhs, dd)))
    }

    mat_fixed_gk_cond <- .eval_on_times(list(fixed_rhs), u_cond)
    mat_id_gk_cond <- if (length(id_rhs_list) > 0) .eval_on_times(id_rhs_list, u_cond) else matrix(0, 15, sd$R_id)
    mat_marker_gk_cond <- if (length(nested$mk_rhs_list) > 0) .eval_on_times(nested$mk_rhs_list, u_cond) else matrix(0, 15, sd$R_mk)
    mat_marker_id_gk_cond <- if (length(nested$idm_rhs_list) > 0) .eval_on_times(nested$idm_rhs_list, u_cond) else matrix(0, 15, sd$Q_idm)

    mat_fixed_gk_cond_fwd <- .eval_on_times(list(fixed_rhs), u_cond_fwd)
    mat_id_gk_cond_fwd <- if (length(id_rhs_list) > 0) .eval_on_times(id_rhs_list, u_cond_fwd) else matrix(0, 15, sd$R_id)
    mat_marker_gk_cond_fwd <- if (length(nested$mk_rhs_list) > 0) .eval_on_times(nested$mk_rhs_list, u_cond_fwd) else matrix(0, 15, sd$R_mk)
    mat_marker_id_gk_cond_fwd <- if (length(nested$idm_rhs_list) > 0) .eval_on_times(nested$idm_rhs_list, u_cond_fwd) else matrix(0, 15, sd$Q_idm)

    if (!is.null(object$config$Bs_obj)) {
        bs_basis <- as.matrix(predict(object$config$Bs_obj, newx = u_cond))
    } else {
        # Fallback for formula-based or if Bs_obj missing
        bs_basis <- .make_basehaz_basis(u_cond, basis = "bs", knots = knots, degree = 3, boundary = c(0, 1))
    }
    mat_basis_gk_cond <- sweep(bs_basis, 2, col_means, "-")

    # Get unique markers for this subject
    markers_subject <- unique(dL[[marker_var]])
    n_markers_subject <- length(markers_subject)
    
    # Create prediction grid: expand time grid to include all markers
    # Each (time, marker) pair gets a row in the prediction data
    n_obs_pred_per_marker <- length(t_grid)
    n_obs_pred <- n_obs_pred_per_marker * n_markers_subject
    
    if (n_obs_pred > 0) {
        # Create expanded data frame with (time, marker) combinations
        # Start with copies of first row to preserve structure
        dl_pred <- dL[rep(1, n_obs_pred), , drop = FALSE]
        
        # Fill in prediction grid
        pred_idx <- 1
        idx_marker_pred_vec <- integer(n_obs_pred)
        for (m_idx in seq_along(markers_subject)) {
            for (t_idx in seq_along(t_grid)) {
                dl_pred[pred_idx, time_var] <- t_grid[t_idx]
                dl_pred[pred_idx, marker_var] <- as.character(markers_subject[m_idx])
                idx_marker_pred_vec[pred_idx] <- m_idx
                pred_idx <- pred_idx + 1
            }
        }
        
        # Scale time for Stan model
        dl_pred_scaled <- dl_pred
        dl_pred_scaled[[time_var]] <- dl_pred_scaled[[time_var]] / tmax
        
        mat_fixed_pred <- .mm(fixed_rhs, dl_pred_scaled)
        mat_id_pred <- if (length(id_rhs_list) > 0) {
            do.call(cbind, lapply(id_rhs_list, function(rhs) .mm(rhs, dl_pred_scaled)))
        } else {
            matrix(0, n_obs_pred, sd$R_id)
        }
        mat_marker_pred <- if (length(nested$mk_rhs_list) > 0) {
            do.call(cbind, lapply(nested$mk_rhs_list, function(rhs) .mm(rhs, dl_pred_scaled)))
        } else {
            matrix(0, n_obs_pred, sd$R_mk)
        }
        mat_marker_id_pred <- if (length(nested$idm_rhs_list) > 0) {
            do.call(cbind, lapply(nested$idm_rhs_list, function(rhs) .mm(rhs, dl_pred_scaled)))
        } else {
            matrix(0, n_obs_pred, sd$Q_idm)
        }
        idx_marker_pred <- idx_marker_pred_vec

        dist_sigma_pred <- .build_dist_matrix(dist_formulas$sigma, dl_pred_scaled)
        dist_nu_pred <- .build_dist_matrix(dist_formulas$nu, dl_pred_scaled)
        dist_phi_pred <- .build_dist_matrix(dist_formulas$phi, dl_pred_scaled)
        dist_alpha_pred <- .build_dist_matrix(dist_formulas$alpha, dl_pred_scaled)
        dist_phi_beta_pred <- .build_dist_matrix(dist_formulas$phi_beta, dl_pred_scaled)
        dist_tau_sde_pred <- .build_dist_matrix(dist_formulas$tau_sde, dl_pred_scaled)
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
    mat_basis_gk_surv <- array(0, dim = c(n_times_surv, 15, sd$Kbs))
    mat_fixed_gk_surv <- array(0, dim = c(n_times_surv, 15, sd$P))
    mat_id_gk_surv <- array(0, dim = c(n_times_surv, 15, sd$R_id))
    mat_marker_gk_surv <- array(0, dim = c(n_times_surv, 15, sd$R_mk))
    mat_marker_id_gk_surv <- array(0, dim = c(n_times_surv, 15, sd$Q_idm))
    mat_fixed_gk_surv_fwd <- mat_fixed_gk_surv
    mat_id_gk_surv_fwd <- mat_id_gk_surv
    mat_marker_gk_surv_fwd <- mat_marker_gk_surv
    mat_marker_id_gk_surv_fwd <- mat_marker_id_gk_surv

    if (n_times_surv > 0) {
        for (s in 1:n_times_surv) {
            ts <- t_surv_grid[s] / tmax
            us <- 0.5 * ts * (gk + 1)
            us_f <- us + sd$eps_fd

            if (!is.null(object$config$Bs_obj)) {
                bs <- as.matrix(predict(object$config$Bs_obj, newx = us))
            } else {
                bs <- .make_basehaz_basis(us, basis = "bs", knots = knots, degree = 3, boundary = c(0, 1))
            }
            mat_basis_gk_surv[s, , ] <- sweep(bs, 2, col_means, "-")
            mat_fixed_gk_surv[s, , ] <- .eval_on_times(list(fixed_rhs), us)
            mat_id_gk_surv[s, , ] <- if (length(id_rhs_list) > 0) .eval_on_times(id_rhs_list, us) else matrix(0, 15, sd$R_id)
            mat_marker_gk_surv[s, , ] <- if (length(nested$mk_rhs_list) > 0) .eval_on_times(nested$mk_rhs_list, us) else matrix(0, 15, sd$R_mk)
            mat_marker_id_gk_surv[s, , ] <- if (length(nested$idm_rhs_list) > 0) .eval_on_times(nested$idm_rhs_list, us) else matrix(0, 15, sd$Q_idm)

            mat_fixed_gk_surv_fwd[s, , ] <- .eval_on_times(list(fixed_rhs), us_f)
            mat_id_gk_surv_fwd[s, , ] <- if (length(id_rhs_list) > 0) .eval_on_times(id_rhs_list, us_f) else matrix(0, 15, sd$R_id)
            mat_marker_gk_surv_fwd[s, , ] <- if (length(nested$mk_rhs_list) > 0) .eval_on_times(nested$mk_rhs_list, us_f) else matrix(0, 15, sd$R_mk)
            mat_marker_id_gk_surv_fwd[s, , ] <- if (length(nested$idm_rhs_list) > 0) .eval_on_times(nested$idm_rhs_list, us_f) else matrix(0, 15, sd$Q_idm)
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
    y_var <- eval(object$call$y_var) %||% "y"
    if (!(y_var %in% colnames(dL))) stop(paste("Response variable", y_var, "not found in dL."))

    # Map markers carefully using fitted levels
    marker_levels <- object$stan_data$marker_levels
    if (is.null(marker_levels)) {
        marker_levels <- levels(dL[[marker_var]]) %||% sort(unique(dL[[marker_var]]))
    }

    # Marker weights (use fitted weights when available)
    marker_weights <- object$stan_data$marker_weights
    if (is.null(marker_weights)) {
        marker_weights <- rep(1 / length(marker_levels), length(marker_levels))
    } else {
        if (!is.null(names(marker_weights))) {
            marker_weights <- marker_weights[marker_levels]
        }
        if (length(marker_weights) != length(marker_levels) || any(!is.finite(marker_weights)) || any(marker_weights < 0)) {
            marker_weights <- rep(1 / length(marker_levels), length(marker_levels))
        } else {
            sw <- sum(marker_weights)
            if (!is.finite(sw) || sw <= 0) {
                marker_weights <- rep(1 / length(marker_levels), length(marker_levels))
            } else {
                marker_weights <- marker_weights / sw
            }
        }
    }

    marker_weights_draws <- draws_list$marker_weights_draws
    if (is.null(marker_weights_draws) || nrow(marker_weights_draws) != length(draws_list$sigma_y_shared)) {
        marker_weights_draws <- matrix(rep(marker_weights, each = length(draws_list$sigma_y_shared)),
                                        nrow = length(draws_list$sigma_y_shared), byrow = TRUE)
    }

    dL[[marker_var]] <- factor(dL[[marker_var]], levels = marker_levels)
    marker_int <- as.integer(dL[[marker_var]])
    if (any(is.na(marker_int))) stop("Markers in newdataLong don't match fitted model levels.")

    # Trials (optional) for binomial predictions
    trials_obs <- rep.int(1L, nrow(dL))
    if ("trials" %in% colnames(dL)) {
        trials_obs <- as.integer(dL$trials)
    }
    trials_pred <- rep.int(1L, n_obs_pred)
    if (n_obs_pred > 0 && "trials" %in% colnames(dL)) {
        for (m_idx in seq_along(markers_subject)) {
            m_val <- markers_subject[m_idx]
            m_trials <- dL$trials[dL[[marker_var]] == m_val]
            tr_val <- if (length(m_trials) > 0) m_trials[1] else 1L
            trials_pred[idx_marker_pred == m_idx] <- as.integer(tr_val)
        }
    }

    out <- list(
        n_draws = length(draws_list$sigma_y_shared),
        n_obs_long = nrow(dL), idx_marker_obs = marker_int, n_marker_types = sd$D,
        marker_weights = as.numeric(marker_weights),
        marker_weights_draws = marker_weights_draws,
        y_real = as.numeric(dL[[y_var]]),
        y_int = as.integer(dL[[y_var]]),
        trials_obs = as.integer(trials_obs),
        n_fixed_effects = sd$P, n_random_id = sd$R_id, n_random_marker = sd$R_mk, n_random_marker_id = sd$Q_idm,
        mat_fixed_obs = mat_fixed_obs, mat_id_obs = mat_id_obs, mat_marker_obs = mat_marker_obs, mat_marker_id_obs = mat_marker_id_obs,
        n_cov_vcov = sd$K_cov, vec_cov_vcov = as.vector(vec_cov_vcov),
        n_cov_hazard = sd$p_w, vec_cov_hazard = as.vector(vec_cov_hazard),
        n_basehaz_basis = sd$Kbs, time_condition = T_cond_scaled,
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
        tau_marker_id = draws_list$tau_marker_id, Lcorr_marker_id = draws_list$Lcorr_marker_id,
        num_unique_cov_entries = M_cov_val,
        alpha_vcov_reg = draws_list$alpha_vcov_reg, beta_vcov_reg_flat = draws_list$beta_vcov_reg_flat,
        tau_vcov_reg = draws_list$tau_vcov_reg, lambda_vcov_reg = draws_list$lambda_vcov_reg,
        K_event = sd$K_event %||% 1L,
        log_h0_intercept = draws_list$log_h0_intercept, bs_gamma_c = draws_list$bs_gamma_c, gamma_hazard = draws_list$gamma_hazard,
        beta_sigma = draws_list$beta_sigma,
        beta_nu = draws_list$beta_nu,
        beta_phi = draws_list$beta_phi,
        beta_alpha = draws_list$beta_alpha,
        beta_phi_beta = draws_list$beta_phi_beta,
        beta_tau_sde = draws_list$beta_tau_sde,
        sigma_y_shared = draws_list$sigma_y_shared, sigma_marker_specific = draws_list$sigma_marker_specific,
        nu_marker = draws_list$nu_marker, phi_nb_marker = draws_list$phi_nb_marker, alpha_skew_marker = draws_list$alpha_skew_marker,
        phi_beta_marker = draws_list$phi_beta_marker, tau_sde_marker = draws_list$tau_sde_marker,
        family_long = sd$family_long, flag_resid_dim = sd$flag_resid_dim,
        coeff_assoc_cv_total = draws_list$coeff_assoc_cv_total, coeff_assoc_cs_total = draws_list$coeff_assoc_cs_total,
        coeff_assoc_cv_mean = draws_list$coeff_assoc_cv_mean, coeff_assoc_cs_mean = draws_list$coeff_assoc_cs_mean,
        coeff_assoc_cv_marker = draws_list$coeff_assoc_cv_marker, coeff_assoc_cs_marker = draws_list$coeff_assoc_cs_marker,
        coeff_assoc_vcov_var = draws_list$coeff_assoc_vcov_var,
        marker_weights_draws = draws_list$marker_weights_draws,
        marker_weights_draws = draws_list$marker_weights_draws,
        K_ord = sd$K_ord %||% 2L,
        cutpoints_ord = draws_list$cutpoints_ord,
        flag_assoc_cv_total = sd$assoc_cv_total, flag_assoc_cv_mean = sd$assoc_cv_mean, flag_assoc_cv_marker = sd$assoc_cv_marker,
        flag_assoc_cs_total = sd$assoc_cs_total, flag_assoc_cs_mean = sd$assoc_cs_mean, flag_assoc_cs_marker = sd$assoc_cs_marker,
        flag_assoc_vcov = sd$assoc_vcov,
        tf_mode_cv_tot = sd$tf_mode_cv_tot, tf_mode_cs_tot = sd$tf_mode_cs_tot,
        tf_mode_cv_mean = sd$tf_mode_cv_mean, tf_mode_cv_marker = sd$tf_mode_cv_marker,
        tf_mode_cs_mean = sd$tf_mode_cs_mean, tf_mode_cs_marker = sd$tf_mode_cs_marker,
        tf_mode_vcov = sd$tf_mode_vcov,
        n_functional_ops_cv = sd$n_functional_ops_cv, functional_ops_cv = sd$functional_ops_cv,
        n_const_cv = sd$n_const_cv, const_data_cv = sd$const_data_cv,
        n_functional_ops_cs = sd$n_functional_ops_cs, functional_ops_cs = sd$functional_ops_cs,
        n_const_cs = sd$n_const_cs, const_data_cs = sd$const_data_cs,
        n_functional_ops_vcov = sd$n_functional_ops_vcov, functional_ops_vcov = sd$functional_ops_vcov,
        n_const_vcov = sd$n_const_vcov, const_data_vcov = sd$const_data_vcov,
        n_knots_cv = sd$n_knots_cv, knots_cv = sd$knots_cv,
        n_coeff_cv = sd$n_coeff_cv, coeff_cv = sd$coeff_cv, spline_degree_cv = sd$spline_degree_cv,
        n_knots_cs = sd$n_knots_cs, knots_cs = sd$knots_cs,
        n_coeff_cs = sd$n_coeff_cs, coeff_cs = sd$coeff_cs, spline_degree_cs = sd$spline_degree_cs,
        n_knots_vcov = sd$n_knots_vcov, knots_vcov = sd$knots_vcov,
        n_coeff_vcov = sd$n_coeff_vcov, coeff_vcov = sd$coeff_vcov, spline_degree_vcov = sd$spline_degree_vcov,
        n_functional_ops_cv_mean = sd$n_functional_ops_cv_mean, functional_ops_cv_mean = sd$functional_ops_cv_mean,
        n_const_cv_mean = sd$n_const_cv_mean, const_data_cv_mean = sd$const_data_cv_mean,
        n_knots_cv_mean = sd$n_knots_cv_mean, knots_cv_mean = sd$knots_cv_mean,
        n_coeff_cv_mean = sd$n_coeff_cv_mean, coeff_cv_mean = sd$coeff_cv_mean, spline_degree_cv_mean = sd$spline_degree_cv_mean,
        n_functional_ops_cv_marker = sd$n_functional_ops_cv_marker, functional_ops_cv_marker = sd$functional_ops_cv_marker,
        n_const_cv_marker = sd$n_const_cv_marker, const_data_cv_marker = sd$const_data_cv_marker,
        n_knots_cv_marker = sd$n_knots_cv_marker, knots_cv_marker = sd$knots_cv_marker,
        n_coeff_cv_marker = sd$n_coeff_cv_marker, coeff_cv_marker = sd$coeff_cv_marker, spline_degree_cv_marker = sd$spline_degree_cv_marker,
        n_functional_ops_cs_mean = sd$n_functional_ops_cs_mean, functional_ops_cs_mean = sd$functional_ops_cs_mean,
        n_const_cs_mean = sd$n_const_cs_mean, const_data_cs_mean = sd$const_data_cs_mean,
        n_knots_cs_mean = sd$n_knots_cs_mean, knots_cs_mean = sd$knots_cs_mean,
        n_coeff_cs_mean = sd$n_coeff_cs_mean, coeff_cs_mean = sd$coeff_cs_mean, spline_degree_cs_mean = sd$spline_degree_cs_mean,
        n_functional_ops_cs_marker = sd$n_functional_ops_cs_marker, functional_ops_cs_marker = sd$functional_ops_cs_marker,
        n_const_cs_marker = sd$n_const_cs_marker, const_data_cs_marker = sd$const_data_cs_marker,
        n_knots_cs_marker = sd$n_knots_cs_marker, knots_cs_marker = sd$knots_cs_marker,
        n_coeff_cs_marker = sd$n_coeff_cs_marker, coeff_cs_marker = sd$coeff_cs_marker, spline_degree_cs_marker = sd$spline_degree_cs_marker,
        flag_indep_id_re = sd$indep_id_re, flag_indep_marker_re = sd$indep_marker_re,
        flag_indep_marker_byid_latent_re = sd$indep_marker_byid_latent_re,
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

#' @export
summary.JoinMeDynPred <- function(object, ...) {
    if (!inherits(object, "JoinMeDynPred")) {
        cli::cli_abort(c(
            x = "Object must be a {.cls JoinMeDynPred} prediction.",
            i = "Call predict() on a JoinMeFit object first."
        ))
    }
    cached <- object$cache_get("summary")
    if (!is.null(cached)) return(cached)

    n_long <- if (!is.null(object$predictions$longitudinal)) nrow(object$predictions$longitudinal) else 0
    n_surv <- if (!is.null(object$predictions$survival)) nrow(object$predictions$survival) else 0
    ids_long <- if (!is.null(object$predictions$longitudinal)) length(unique(object$predictions$longitudinal$id)) else 0
    ids_surv <- if (!is.null(object$predictions$survival)) length(unique(object$predictions$survival$id)) else 0
    n_cumhaz <- if (!is.null(object$predictions$cumhaz)) nrow(object$predictions$cumhaz) else 0
    ids_cumhaz <- if (!is.null(object$predictions$cumhaz)) length(unique(object$predictions$cumhaz$id)) else 0

    tables <- list(
        longitudinal = data.frame(
            metric = c("rows", "subjects", "scale"),
            value = c(n_long, ids_long, object$metadata$scale %||% NA_character_),
            stringsAsFactors = FALSE
        ),
        survival = data.frame(
            metric = c("rows", "subjects"),
            value = c(n_surv, ids_surv),
            stringsAsFactors = FALSE
        ),
        cumhaz = data.frame(
            metric = c("rows", "subjects"),
            value = c(n_cumhaz, ids_cumhaz),
            stringsAsFactors = FALSE
        )
    )

    summary_obj <- SummaryJoinMeDynPred$new(tables = tables, metadata = object$metadata)
    object$cache_set("summary", summary_obj)
    summary_obj
}

#' @export
print.summary_JoinMeDynPred <- function(x, ...) {
    cat("Prediction summary\n")
    cat("==================\n")
    if (!is.null(x$tables$longitudinal)) {
        cat("Longitudinal\n")
        cat("------------\n")
        print(x$tables$longitudinal, row.names = FALSE)
    }
    if (!is.null(x$tables$survival)) {
        cat("\nSurvival\n")
        cat("--------\n")
        print(x$tables$survival, row.names = FALSE)
    }
    if (!is.null(x$tables$cumhaz)) {
        cat("\nCumulative hazard\n")
        cat("-----------------\n")
        print(x$tables$cumhaz, row.names = FALSE)
    }
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

.resolve_time_grid <- function(times, id, t_cond, tmax_val, default_n = 50, min_points = 50, kind = "longitudinal") {
    grid <- NULL
    if (is.null(times)) {
            if (default_n < min_points) {
                cli::cli_warn(c(
                    x = "Requested {kind} grid size {default_n} is below the minimum of {min_points}.",
                    i = "Using {min_points} time points instead."
                ))
            }
            grid <- seq(t_cond, min(t_cond + 5, tmax_val), length.out = max(default_n, min_points))
    } else if (is.list(times)) {
        key <- as.character(id)
        grid <- times[[key]]
        if (is.null(grid)) {
            cli::cli_warn(c(
                x = "No {kind} time grid found for subject {id}.",
                i = "Falling back to the default time grid."
            ))
            grid <- seq(t_cond, min(t_cond + 5, tmax_val), length.out = max(default_n, min_points))
        }
    } else {
        grid <- times
    }

    grid <- grid[grid >= t_cond]
    if (length(grid) < min_points && default_n == min_points) {
        cli::cli_warn(c(
            x = "{kind} prediction grid has fewer than {min_points} time points for subject {id}.",
            i = "Consider supplying a denser time grid for smoother intervals."
        ))
    }
    grid
}

`%||%` <- function(x, y) if (is.null(x)) y else x
