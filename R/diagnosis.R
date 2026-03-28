#' @name joinme_diagnosis
#' @title Model Diagnosis
#'
#' @description
#' Posterior predictive checks and information criteria.
#'
#' @keywords internal
NULL

# File overview:
# - Extract log-likelihood and compute LOO/WAIC/ELPD summaries.
# - Provide posterior predictive checks and Bayes factor utilities.

# ---- diagnosis generic ---------------------------------------------------

#' Diagnostic summary for joinme objects
#'
#' @description
#' Returns a compact diagnostics table for fitted (`JoinMeFit`) and dynamic
#' prediction (`JoinMeDynPred`) objects.
#'
#' For `JoinMeFit`, diagnostics summarise the Stan sampler run. For
#' `JoinMeDynPred`, diagnostics summarise posterior-draw quality for predicted
#' quantities and are aligned to the same metric schema used for `JoinMeFit`.
#'
#' @param object A joinme object.
#' @param ... Additional arguments passed to class-specific methods.
#'
#' @return A data frame with `metric` and `value` columns.
#' @export
diagnosis <- function(object, ...) {
	UseMethod("diagnosis")
}

#' @export
diagnosis.JoinMeFit <- function(object, ...) {
	assertthat::assert_that(inherits(object, "JoinMeFit"), msg = "Object must be a JoinMeFit instance.")
	diag <- .joinme_sampler_diagnostics(object$fit)
	.diagnostics_table_from_sampler(diag)
}

#' @export
diagnosis.JoinMeDynPred <- function(object, ...) {
	assertthat::assert_that(inherits(object, "JoinMeDynPred"), msg = "Object must be a JoinMeDynPred instance.")
	sum_obj <- summary(object)
	sum_obj$tables$diagnostics %||% .build_common_diagnostics_table()
}

# ---- log-likelihood extraction -------------------------------------------

#' Log-likelihood summary
#'
#' @rdname log_lik.JoinMeFit
#' @param object A fitted object of class `JoinMeFit`.
#' @param what Character; which component to return: `"long"`, `"surv"`, or `"total"`.
#' @param draws Optional number of posterior draws to subset.
#' @param seed Random seed for draw subsetting.
#' @param ... Unused.
#'
#' @importFrom rstantools log_lik
#' @return A matrix
#' @seealso [rstantools::log_lik()]
#' @export
log_lik.JoinMeFit <- function(object, what = c("long", "surv", "total"), draws = NULL, seed = 1, ...) {
	assertthat::assert_that(inherits(object, "JoinMeFit"), msg = "Object must be a JoinMeFit instance.")
	what <- match.arg(what)
	sd <- object$stan_data

	vars_long <- paste0("log_lik_long[", seq_len(sd$N), "]")
	vars_surv <- paste0("log_lik_surv[", seq_len(sd$n_id), "]")
	vars <- switch(what,
								 long = vars_long,
								 surv = vars_surv,
								 total = c(vars_long, vars_surv))

	ddf <- .get_draws_df(object$fit, variables = vars, draws = draws, seed = seed)
	mat <- as.matrix(ddf[, vars, drop = FALSE])
	attr(mat, "component") <- what
	mat
}

#' Backup, unused
#' @noRd
#' @keywords internal
unused_loo.JoinMeFit <- function(object, what = c("total", "long", "surv"), draws = NULL, seed = 1, ...) {
	what <- match.arg(what)
	ll_mat <- log_lik.JoinMeFit(object, what = what, draws = draws, seed = seed)

	log_mean_exp <- function(x) {
		m <- max(x)
		m + log(mean(exp(x - m)))
	}

	log_colMeansExp <- function(mat) {
		S <- nrow(mat)
		matrixStats::colLogSumExps(mat, na.rm=TRUE) - log(S)
	}

	# pointwise <- apply(ll_mat, 2, log_mean_exp)
	pointwise <- matrix(log_colMeansExp(ll_mat), nrow = 1)
	elpd <- sum(pointwise)
	se_elpd <- sqrt(ncol(ll_mat) * var(pointwise))
	output <- structure(
		list(
			pointwise = structure(pointwise, component = 'elpd_loo'),
			estimates = structure(
				cbind(elpd, se_elpd),
				.dimnames = list('elpd_loo', c("Estimate", "SE")))
		), class = c('psis_loo', 'loo'))
	attr(output, 'what') <- what
	output
}

# ---- LOO / WAIC / ELPD ---------------------------------------------------

#' LOO-CV for JoinMe models
#' 
#' @rdname loo.JoinMeFit
#' @param object A fitted object of class `JoinMeFit`.
#' @param what Character; which component to return: `"long"`, `"surv"`, or `"total"`.
#' @param draws Optional number of posterior draws to subset.
#' @param seed Random seed for draw subsetting.
#' @param ... Additional arguments passed to `loo::loo()`.
#'
#' @return A `loo` object.
#' @importFrom loo loo
#' @export
loo.JoinMeFit <- function(object, what = c("total", "long", "surv"), draws = NULL, seed = 1, ...) {
	what <- match.arg(what)
	ll_mat <- log_lik.JoinMeFit(object, what = what, draws = draws, seed = seed)
	loo::loo(ll_mat, ...)
}

#' WAIC for JoinMe models
#' @rdname waic.JoinMeFit
#' @param object A fitted object of class `JoinMeFit`.
#' @param what Character; which component to return: `"long"`, `"surv"`, or `"total"`.
#' @param draws Optional number of posterior draws to subset.
#' @param seed Random seed for draw subsetting.
#' @param ... Additional arguments passed to `loo::waic()`.
#'
#' @importFrom loo waic
#' 
#' @return A `waic` object.
#' @export
waic.JoinMeFit <- function(object, what = c("total", "long", "surv"), draws = NULL, seed = 1, ...) {
	if (!requireNamespace("loo", quietly = TRUE)) {
		cli::cli_abort("Package {.pkg loo} is required for waic(). Install it with install.packages('loo').")
	}
	what <- match.arg(what)
	ll_mat <- log_lik.JoinMeFit(object, what = what, draws = draws, seed = seed)
	loo::waic(ll_mat, ...)
}


#' ELPD summary for JoinMe models
#' 
#' @param object A fitted object of class `JoinMeFit`.
#' @param what Character; which component to return: `"long"`, `"surv"`, or `"total"`.
#' @param draws Optional number of posterior draws to subset.
#' @param seed Random seed for draw subsetting.
#' @param ... Additional arguments passed to `loo::elpd()`.
#' 
#' @importFrom loo elpd
#' @rdname elpd.JoinMeFit
#' @return A data frame with ELPD and standard error.
#' @export
elpd.JoinMeFit <- function(object, what = c("total", "long", "surv"), draws = NULL, seed = 1, ...) {
	loo_obj <- loo.JoinMeFit(object, what = what, draws = draws, seed = seed, ...)
	structure(
		data.frame(
			elpd_loo = loo_obj$estimates["elpd_loo", "Estimate"],
			se_elpd_loo = loo_obj$estimates["elpd_loo", "SE"],
			row.names = NULL
		),
		class = c("elpd_generic", "loo")
	)
}

# ---- Bayes factor (bridge sampling) --------------------------------------

#' Bayes factor between two JoinMe fits
#'
#' @param fit1 First fitted `JoinMeFit` object.
#' @param fit2 Second fitted `JoinMeFit` object.
#' @param ... Additional arguments passed to `bridgesampling::bridge_sampler()`.
#'
#' @return A list with bridge sampling results and Bayes factor.
#' @export
bayes_factor <- function(fit1, fit2, ...) {
	assertthat::assert_that(inherits(fit1, "JoinMeFit"), msg = "fit1 must be a JoinMeFit instance.")
	assertthat::assert_that(inherits(fit2, "JoinMeFit"), msg = "fit2 must be a JoinMeFit instance.")

	if (!requireNamespace("bridgesampling", quietly = TRUE)) {
		cli::cli_abort("Package {.pkg bridgesampling} is required for bayes_factor(). Install it with install.packages('bridgesampling').")
	}

	bs1 <- tryCatch(
		bridgesampling::bridge_sampler(fit1$fit, ...),
		error = function(e) {
			cli::cli_abort(c(
				x = "Bridge sampling failed for fit1.",
				i = "Ensure the cmdstanr fit is compatible with bridgesampling.",
				i = "Error: {e$message}"
			))
		}
	)

	bs2 <- tryCatch(
		bridgesampling::bridge_sampler(fit2$fit, ...),
		error = function(e) {
			cli::cli_abort(c(
				x = "Bridge sampling failed for fit2.",
				i = "Ensure the cmdstanr fit is compatible with bridgesampling.",
				i = "Error: {e$message}"
			))
		}
	)

	bf <- bridgesampling::bf(bs1, bs2)
	list(bridge1 = bs1, bridge2 = bs2, bayes_factor = bf)
}

# ---- Posterior predictive checks -----------------------------------------

#' Posterior predictive check for longitudinal outcomes
#'
#' @param object A fitted object of class `JoinMeFit`.
#' @param newdataLong New longitudinal data frame for predictions.
#' @param newdataEvent New event data frame for predictions.
#' @param ci_level Numeric; credible interval level (default 0.95).
#' @param n_samples Integer; number of posterior samples to use (default 200).	
#' @param ... Additional arguments.
#'
#' @importFrom bayesplot pp_check
#' @rdname pp_check.JoinMeFit
#' @details 
#' This crap is still under development.
#' @return A list with observation-level summaries and overall diagnostics.
#' @export
pp_check.JoinMeFit <- function(object, newdataLong = NULL, newdataEvent = NULL, ci_level = 0.95, n_samples = 200, seed = 123, plot = FALSE,...) {
	# Workflow: run posterior_epred -> summarise coverage and RMSE
	assertthat::assert_that(inherits(object, "JoinMeFit"), msg = "Object must be a JoinMeFit instance.")
	if (is.null(newdataLong) && is.null(newdataEvent)) {
		if (!is.null(object$dataLong) && !is.null(object$dataEvent)) {
			newdataLong <- object$dataLong
			newdataEvent <- object$dataEvent
		} else {
			cli::cli_abort(c(
				x = "Training data not found on the fitted object.",
				i = "Provide {.arg newdataLong} and {.arg newdataEvent} explicitly."
			))
		}
	}
	if (is.null(newdataLong) != is.null(newdataEvent)) {
		cli::cli_abort(c(
			x = "Both {.arg newdataLong} and {.arg newdataEvent} are required.",
			i = "Provide both or neither to use the training data."
		))
	}

	ci_level <- .validate_ci_levels(ci_level)
	probs <- .quantile_probs_from_ci(ci_level)

	pred <- posterior_epred(object,
				newdataLong = newdataLong,
				newdataEvent = newdataEvent,
				control = list(n_samples = n_samples),
				seed = seed,
				...)

	draws_fit <- pred$draws$longitudinal_fitted
	if (is.null(draws_fit) || length(draws_fit) == 0) {
		cli::cli_abort("No fitted draws found. Ensure longitudinal data is provided.")
	}

	id_var <- pred$metadata$id_var %||% (eval(object$call$id_var) %||% "id")
	time_var <- pred$metadata$time_var %||% (eval(object$call$time_var) %||% "time")
	marker_var <- pred$metadata$marker_var %||% (eval(object$call$marker_var) %||% "marker")
	y_var <- pred$metadata$response_var %||% .resolve_response_var(
		formulaLong = object$formulaLong,
		dataLong = newdataLong,
		context = "joinme_diagnosis()"
	)

	out_list <- list()
	for (id in names(draws_fit)) {
		# Per-subject summary of fitted draws vs observed outcomes
		dL <- newdataLong[newdataLong[[id_var]] == id, , drop = FALSE]
		if (nrow(dL) == 0) next

		epred_mat <- draws_fit[[id]]$epred
		if (is.null(epred_mat) || nrow(epred_mat) == 0) next

		mean_epred <- colMeans(epred_mat)
		sd_epred <- apply(epred_mat, 2, stats::sd)
		lwr <- apply(epred_mat, 2, stats::quantile, probs = probs[1], names = FALSE)
		upr <- apply(epred_mat, 2, stats::quantile, probs = probs[2], names = FALSE)

		df <- data.frame(
			id = dL[[id_var]],
			time = dL[[time_var]],
			marker = dL[[marker_var]],
			observed = dL[[y_var]],
			pred_mean = mean_epred,
			pred_sd = sd_epred,
			pred_lwr = lwr,
			pred_upr = upr
		)
		df$covered <- with(df, observed >= pred_lwr & observed <= pred_upr)
		out_list[[id]] <- df
	}

	out_df <- do.call(rbind, out_list)
	overall <- list(
		coverage = mean(out_df$covered, na.rm = TRUE),
		rmse = sqrt(mean((out_df$observed - out_df$pred_mean)^2, na.rm = TRUE)),
		n = nrow(out_df)
	)

	list(summary = overall, data = out_df)
}

#' @keywords internal
.resolve_train_data <- function(object, newdataLong, newdataEvent, purpose = "analysis") {
	if (is.null(newdataLong) && is.null(newdataEvent)) {
		if (!is.null(object$dataLong) && !is.null(object$dataEvent)) {
			newdataLong <- object$dataLong
			newdataEvent <- object$dataEvent
		} else {
			cli::cli_abort(c(
				x = "Training data not found on the fitted object.",
				i = "Provide {.arg newdataLong} and {.arg newdataEvent} for {purpose}."
			))
		}
	}
	if (is.null(newdataLong) != is.null(newdataEvent)) {
		cli::cli_abort(c(
			x = "Both {.arg newdataLong} and {.arg newdataEvent} are required.",
			i = "Provide both or neither to use the training data."
		))
	}
	list(newdataLong = newdataLong, newdataEvent = newdataEvent)
}

# ---- Concordance ----------------------------------------------------------

#' Time-varying concordance for the survival component
#'
#' @rdname concordance.JoinMeFit
#' @param object A fitted object of class `JoinMeFit`.
#' @param newdataLong Longitudinal data for evaluation (defaults to training data).
#' @param newdataEvent Event data for evaluation (defaults to training data).
#' @param time_start Numeric landmark time(s) or a column name in `newdataEvent`.
#' @param time_horizon Numeric horizon time(s) or a column name in `newdataEvent`.
#' @param Dt Numeric horizon width; used when `time_horizon` is not supplied.
#' @param cause Integer; cause index for competing risks (default 1).
#' @param n_samples Number of posterior draws for prediction (default 200).
#' @param seed Random seed for draw subsetting.
#' @param type_weights Time-weighting for concordance; passed to `survival::concordance`.
#' @param ... Additional arguments passed to `predict.JoinMeFit()`.
#'
#' @importFrom survival Surv concordance
#' @return A data frame with time-varying concordance at each landmark time.
#' @export
concordance.JoinMeFit <- function(object, newdataLong = NULL, newdataEvent = NULL, time_start,
                                  time_horizon = NULL, Dt = NULL, cause = 1, n_samples = 200,
                                  seed = 123, type_weights = "none", ...) {
	assertthat::assert_that(inherits(object, "JoinMeFit"), msg = "Object must be a JoinMeFit instance.")
	if (missing(time_start)) {
		cli::cli_abort("{.arg time_start} is required for time-varying concordance.")
	}
	data_in <- .resolve_train_data(object, newdataLong, newdataEvent, purpose = "time-varying concordance")
	newdataLong <- data_in$newdataLong
	newdataEvent <- data_in$newdataEvent

	if (!is.numeric(cause) || length(cause) != 1) {
		cli::cli_abort("{.arg cause} must be a single integer.")
	}

	id_var <- eval(object$call$id_var) %||% "id"
	time_var <- eval(object$call$time_var) %||% "time"
	event_vars <- .resolve_event_model_vars(
		formulaEvent = object$formulaEvent,
		dataEvent = newdataEvent,
		context = "concordance.JoinMeFit()"
	)
	event_time <- as.numeric(event_vars$event_time)
	event_status <- event_vars$event_status

	if (!(id_var %in% names(newdataEvent))) {
		cli::cli_abort("{.arg newdataEvent} must include id column {id_var}.")
	}
	if (length(event_time) != nrow(newdataEvent) || length(event_status) != nrow(newdataEvent)) {
		cli::cli_abort("Could not derive event time/status vectors aligned to {.arg newdataEvent}.")
	}
	if (!(id_var %in% names(newdataLong)) || !(time_var %in% names(newdataLong))) {
		cli::cli_abort("{.arg newdataLong} must include id and time columns used in the model.")
	}

	ids <- unique(newdataEvent[[id_var]])
	ids <- ids[!is.na(ids)]
	if (length(ids) == 0) {
		cli::cli_abort("No subjects found in {.arg newdataEvent}.")
	}

	build_time_map <- function(val, label) {
		if (is.character(val)) {
			if (length(val) != 1) {
				cli::cli_abort("{label} must be a single column name when character.")
			}
			if (!(val %in% names(newdataEvent))) {
				cli::cli_abort("{label} column {val} not found in newdataEvent.")
			}
			return(stats::setNames(as.numeric(newdataEvent[[val]]), newdataEvent[[id_var]]))
		}
		if (!is.numeric(val)) {
			cli::cli_abort("{label} must be numeric or a column name.")
		}
		if (length(val) == 1) {
			return(stats::setNames(rep(as.numeric(val), length(ids)), ids))
		}
		if (!is.null(names(val))) {
			return(stats::setNames(as.numeric(val), names(val)))
		}
		if (length(val) == length(ids)) {
			return(stats::setNames(as.numeric(val), ids))
		}
		NULL
	}

	use_grid <- is.numeric(time_start) && length(time_start) > 1 && is.null(names(time_start))
	if (use_grid) {
		time_start_grid <- as.numeric(time_start)
		if (!is.null(time_horizon) && length(time_horizon) > 1 && length(time_horizon) != length(time_start_grid)) {
			cli::cli_abort("{.arg time_horizon} must be length 1 or the same length as {.arg time_start}.")
		}
		if (is.null(time_horizon) && is.null(Dt)) {
			cli::cli_abort("Provide {.arg time_horizon} or {.arg Dt} when {.arg time_start} is a vector.")
		}
		results <- lapply(seq_along(time_start_grid), function(i) {
			th <- if (!is.null(time_horizon)) {
				if (length(time_horizon) == 1) time_horizon else time_horizon[i]
			} else {
				time_start_grid[i] + Dt
			}
			.time_varying_concordance_single(object, newdataLong, newdataEvent, time_start_grid[i], th, cause, n_samples, seed, type_weights, ...)
		})
		out <- do.call(rbind, results)
		class(out) <- c("tvConcordance_JoinMeFit", "data.frame")
		return(out)
	}

	time_start_map <- build_time_map(time_start, "time_start")
	if (is.null(time_start_map)) {
		cli::cli_abort("{.arg time_start} must be scalar, per-subject, or a column name.")
	}

	if (is.null(time_horizon) && is.null(Dt)) {
		cli::cli_abort("Provide {.arg time_horizon} or {.arg Dt}.")
	}
	if (!is.null(time_horizon)) {
		time_horizon_map <- build_time_map(time_horizon, "time_horizon")
		if (is.null(time_horizon_map)) {
			cli::cli_abort("{.arg time_horizon} must be scalar, per-subject, or a column name.")
		}
	} else {
		if (!is.numeric(Dt) || length(Dt) != 1) {
			cli::cli_abort("{.arg Dt} must be a single numeric value.")
		}
		time_horizon_map <- time_start_map + as.numeric(Dt)
	}

	out <- .time_varying_concordance_single(object, newdataLong, newdataEvent, time_start_map, time_horizon_map, cause, n_samples, seed, type_weights, ...)
	class(out) <- c("tvConcordance_JoinMeFit", "data.frame")
	out
}

#' @keywords internal
.time_varying_concordance_single <- function(object, newdataLong, newdataEvent, time_start, time_horizon, cause, n_samples, seed, type_weights, ...) {
	id_var <- eval(object$call$id_var) %||% "id"
	time_var <- eval(object$call$time_var) %||% "time"
	event_vars <- .resolve_event_model_vars(
		formulaEvent = object$formulaEvent,
		dataEvent = newdataEvent,
		context = "concordance.JoinMeFit()"
	)
	event_time <- as.numeric(event_vars$event_time)
	event_status <- event_vars$event_status

	ids <- unique(newdataEvent[[id_var]])
	ids <- ids[!is.na(ids)]

	if (length(time_start) == 1 && !is.null(names(time_start))) time_start <- unname(time_start)
	if (length(time_horizon) == 1 && !is.null(names(time_horizon))) time_horizon <- unname(time_horizon)

	  if (length(time_start) == 1) {
			time_start_map <- stats::setNames(rep(as.numeric(time_start), length(ids)), ids)
	} else {
			time_start_map <- time_start
	}
	if (length(time_horizon) == 1) {
			time_horizon_map <- stats::setNames(rep(as.numeric(time_horizon), length(ids)), ids)
	} else {
			time_horizon_map <- time_horizon
	}

		if (any(time_horizon_map <= time_start_map, na.rm = TRUE)) {
		cli::cli_abort("{.arg time_horizon} must be greater than {.arg time_start} for all subjects.")
	}

	# Restrict histories to each subject's landmark time
	dL_split <- split(newdataLong, newdataLong[[id_var]])
	dL_filtered <- lapply(names(dL_split), function(id) {
		dd <- dL_split[[id]]
		dd[dd[[time_var]] <= time_start_map[[id]], , drop = FALSE]
	})
	newdataLong <- do.call(rbind, dL_filtered)
	if (is.null(newdataLong) || nrow(newdataLong) == 0) {
		cli::cli_abort("No longitudinal history available at or before the landmark time(s).")
	}

	# Prepare per-subject time grids
	times_map <- stats::setNames(lapply(ids, function(id) {
		t_start_i <- as.numeric(time_start_map[[id]])
		t_horizon_i <- as.numeric(time_horizon_map[[id]])
		seq(t_start_i, t_horizon_i, length.out = 50L)
	}), ids)

	pred <- predict(
		object,
		newdataLong = newdataLong,
		newdataEvent = newdataEvent,
		process = "event",
		times = times_map,
		time_start = time_start_map,
		control = list(n_samples = n_samples),
		seed = seed,
		...
	)

	surv_df <- pred$predictions$survival
	if (is.null(surv_df) || nrow(surv_df) == 0) {
		cli::cli_abort("No survival predictions returned for time-varying concordance.")
	}

	key_df <- data.frame(
		id = ids,
		time_horizon = as.numeric(time_horizon_map[ids]),
		stringsAsFactors = FALSE
	)
	colnames(key_df)[1] <- id_var

	merge_df <- merge(key_df, surv_df, by.x = id_var, by.y = "id", all.x = TRUE)
	merge_df <- merge_df[abs(merge_df$time - merge_df$time_horizon) <= sqrt(.Machine$double.eps), , drop = FALSE]
	if (nrow(merge_df) == 0) {
		cli::cli_abort("Failed to align survival predictions at the requested horizon.")
	}

	merge_df$risk <- 1 - merge_df$Survival

	event_df <- data.frame(
		event_time = as.numeric(event_time),
		event_status = event_status,
		stringsAsFactors = FALSE
	)
	event_df[[id_var]] <- newdataEvent[[id_var]]
	event_df <- event_df[, c(id_var, "event_time", "event_status"), drop = FALSE]
	event_decoded <- .derive_event_outcomes(event_df$event_status, context = "concordance.JoinMeFit()")
	event_df$event_status <- event_decoded$d_event
	event_df$event_type <- event_decoded$event_type
	merge_df <- merge(merge_df, event_df, by = id_var, all.x = TRUE)
	merge_df$time_start <- as.numeric(time_start_map[merge_df[[id_var]]])

	merge_df <- merge_df[merge_df$event_time > merge_df$time_start, , drop = FALSE]
	if (nrow(merge_df) == 0) {
		return(data.frame(time_start = as.numeric(time_start)[1], time_horizon = as.numeric(time_horizon)[1],
			concordance = NA_real_, n_cases = 0, n_controls = 0, n_pairs = 0))
	}

	status <- as.integer(merge_df$event_status)
	status <- ifelse(status == 1 & merge_df$event_type == cause, 1L, 0L)

	merge_df$event_window <- as.integer(status == 1 & merge_df$event_time <= merge_df$time_horizon)
	merge_df$time_window <- pmin(merge_df$event_time, merge_df$time_horizon) - merge_df$time_start

	n_cases <- sum(merge_df$event_window == 1, na.rm = TRUE)
	n_controls <- sum(merge_df$event_window == 0 & merge_df$event_time > merge_df$time_horizon, na.rm = TRUE)

	if (n_cases == 0 || n_controls == 0) {
		return(data.frame(time_start = as.numeric(time_start)[1], time_horizon = as.numeric(time_horizon)[1],
			concordance = NA_real_, n_cases = n_cases, n_controls = n_controls, n_pairs = 0))
	}

	concord_data <- merge_df[merge_df$event_time > merge_df$time_start, , drop = FALSE]

	weight_map <- list(none = "n")
	timewt <- weight_map[[type_weights]] %||% type_weights

	conc <- survival::concordance(
		survival::Surv(concord_data$time_window, concord_data$event_window) ~ concord_data$risk,
		timewt = timewt
	)
	comparable_pairs <- if (!is.null(conc$count)) {
		sum(as.numeric(conc$count), na.rm = TRUE)
	} else {
		NA_real_
	}

	data.frame(
		time_start = as.numeric(time_start)[1],
		time_horizon = as.numeric(time_horizon)[1],
		concordance = conc$concordance,
		n_cases = n_cases,
		n_controls = n_controls,
		n_pairs = comparable_pairs
	)
}

# ---- Stan diagnostics ----------------------------------------------------

#' Stan diagnostics for JoinMe models
#'
#' @name stan_diagnostics.JoinMeFit
#' @rdname stan_diagnostics.JoinMeFit
#' @param object A fitted object of class `JoinMeFit`.
#' @param pars Optional character vector of parameter names to include.
#' @param regex_pars Optional regular expression for parameter selection.
#' @param draws Optional number of posterior draws to subset.
#' @param seed Random seed for draw subsetting.
#' @param type Diagnostic type (for ESS and MCSE).
#' @param ... Unused.
#'
#' @return A data frame of diagnostics by parameter.
#' @export
stan_rhat.JoinMeFit <- function(object, pars = NULL, regex_pars = NULL, draws = NULL, seed = 1, ...) {
	assertthat::assert_that(inherits(object, "JoinMeFit"), msg = "Object must be a JoinMeFit instance.")
	draws_obj <- .get_draws_obj(object$fit)
	vars <- posterior::variables(draws_obj)
	vars <- .filter_diag_vars(vars, pars, regex_pars)
	if (length(vars) == 0) {
		return(data.frame(variable = character(0), rhat = numeric(0)))
	}
	draws_obj <- .get_draws_obj(object$fit, variables = vars, draws = draws, seed = seed)
	.diag_summary_df(draws_obj, metric = "rhat", vars = vars)
}

#' @rdname stan_diagnostics.JoinMeFit
#' @export
stan_ess.JoinMeFit <- function(object, pars = NULL, regex_pars = NULL, draws = NULL, seed = 1,
                               type = c("bulk", "tail"), ...) {
	assertthat::assert_that(inherits(object, "JoinMeFit"), msg = "Object must be a JoinMeFit instance.")
	type <- match.arg(type)
	draws_obj <- .get_draws_obj(object$fit)
	vars <- posterior::variables(draws_obj)
	vars <- .filter_diag_vars(vars, pars, regex_pars)
	if (length(vars) == 0) {
		return(data.frame(variable = character(0), ess = numeric(0)))
	}
	draws_obj <- .get_draws_obj(object$fit, variables = vars, draws = draws, seed = seed)
	metric_name <- if (type == "bulk") "ess_bulk" else "ess_tail"
	metric_df <- .diag_summary_df(draws_obj, metric = metric_name, vars = vars)
	names(metric_df)[2] <- "ess"
	metric_df
}

#' @rdname stan_diagnostics.JoinMeFit
#' @export
stan_mcse.JoinMeFit <- function(object, pars = NULL, regex_pars = NULL, draws = NULL, seed = 1,
                                type = c("mean", "sd", "median"), ...) {
	assertthat::assert_that(inherits(object, "JoinMeFit"), msg = "Object must be a JoinMeFit instance.")
	type <- match.arg(type)
	draws_obj <- .get_draws_obj(object$fit)
	vars <- posterior::variables(draws_obj)
	vars <- .filter_diag_vars(vars, pars, regex_pars)
	if (length(vars) == 0) {
		return(data.frame(variable = character(0), mcse = numeric(0)))
	}
	draws_obj <- .get_draws_obj(object$fit, variables = vars, draws = draws, seed = seed)
	metric_name <- switch(type,
		mean = "mcse_mean",
		sd = "mcse_sd",
		median = "mcse_median"
	)
	metric_df <- .diag_summary_df(draws_obj, metric = metric_name, vars = vars)
	names(metric_df)[2] <- "mcse"
	metric_df
}

#' @keywords internal

.diag_summary_df <- function(draws_obj, metric, vars) {
	metric_df <- suppressWarnings(tryCatch(
		posterior::summarise_draws(draws_obj, metric),
		error = function(e) NULL
	))
	if (is.null(metric_df) || !all(c("variable", metric) %in% names(metric_df))) {
		cli::cli_abort("Failed to summarise posterior diagnostics for plotting.")
	}
	metric_df <- as.data.frame(metric_df[, c("variable", metric), drop = FALSE], row.names = NULL, stringsAsFactors = FALSE)
	if (nrow(metric_df) == 0L) {
		out <- data.frame(variable = character(0), value = numeric(0), row.names = NULL)
		names(out)[2] <- metric
		return(out)
	}
	if (!is.null(vars)) {
		metric_df <- metric_df[match(vars, metric_df$variable), , drop = FALSE]
		if (nrow(metric_df) != length(vars) || anyNA(metric_df$variable)) {
			cli::cli_abort("Posterior diagnostic output does not align with selected parameters.")
		}
	}
	out <- data.frame(variable = metric_df$variable, value = metric_df[[metric]], row.names = NULL)
	names(out)[2] <- metric
	out
}

#' @keywords internal
.filter_diag_vars <- function(vars, pars = NULL, regex_pars = NULL) {
	if (!is.null(pars)) {
		vars <- intersect(vars, pars)
	}
	if (!is.null(regex_pars)) {
		vars <- vars[grepl(regex_pars, vars)]
	}
	vars
}
