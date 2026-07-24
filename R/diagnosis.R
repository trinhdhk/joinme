#' @name joinme_diagnosis
#' @title Model Diagnosis
#'
#' @description
#' Posterior predictive checks and information criteria.
#'
#' @keywords internal
#' @noRd
NULL

# File overview:
# - Extract log-likelihood and compute LOO/WAIC/ELPD summaries.
# - Provide posterior predictive checks and Bayes factor utilities.

# ---- diagnosis generic ---------------------------------------------------

#' Diagnostic summary for JoiNMe objects
#'
#' @description
#' Returns diagnostic summaries for fitted (`JoiNMeFit`) and dynamic
#' prediction (`JoiNMeDynPred`) objects.
#'
#' For `JoiNMeFit`, diagnostics include both an overall summary table and a
#' parameter-level diagnostics table built from the same cached posterior
#' summaries used by [summary.JoiNMeFit()]. This avoids the slower backend-wide
#' diagnostic pass and keeps the reported metrics aligned with the summary
#' sections users already inspect.
#'
#' For `JoiNMeDynPred`, diagnostics summarise posterior-draw quality for
#' predicted quantities and are aligned to the same metric schema used for
#' `JoiNMeFit`.
#'
#' @param object A JoiNMe object.
#' @param ... Additional arguments passed to class-specific methods.
#'
#' @return
#' For `JoiNMeFit`, a `JoiNMe_diagnosis` object with components `summary` and
#' `by_parameter`. For `JoiNMeDynPred`, a data frame with `metric` and `value`
#' columns.
#' @export
diagnosis <- function(object, ...) {
	UseMethod("diagnosis")
}

#' @keywords internal
#' @noRd
.diagnosis_parameter_label <- function(df) {
	label <- rep("", nrow(df))
	if ("term" %in% names(df)) {
		label <- as.character(df$term)
	}
	if ("parameter" %in% names(df)) {
		param_label <- as.character(df$parameter)
		empty_idx <- !nzchar(trimws(label))
		label[empty_idx] <- param_label[empty_idx]
		label[!empty_idx] <- paste0(
			param_label[!empty_idx],
			": ",
			label[!empty_idx]
		)
	}
	if ("channel" %in% names(df)) {
		label <- ifelse(
			nzchar(trimws(label)),
			paste0(df$channel, ": ", label),
			as.character(df$channel)
		)
	}
	if ("block" %in% names(df)) {
		label <- ifelse(
			nzchar(trimws(label)),
			paste0(df$block, ": ", label),
			as.character(df$block)
		)
	}
	if (all(c("row", "col") %in% names(df))) {
		rc_label <- ifelse(
			is.finite(df$row) & is.finite(df$col),
			paste0("[", df$row, ",", df$col, "]"),
			""
		)
		label <- ifelse(
			nzchar(rc_label),
			ifelse(nzchar(trimws(label)), paste0(label, " ", rc_label), rc_label),
			label
		)
	}
	label[!nzchar(trimws(label))] <- paste0(
		"parameter_",
		seq_len(sum(!nzchar(trimws(label))))
	)
	trimws(label)
}

#' @keywords internal
#' @noRd
.diagnosis_parameter_table_from_tables <- function(tables) {
	metric_cols <- c(
		"Estimate",
		"Est.Error",
		"Q2.5",
		"Q97.5",
		"Rhat",
		"ess_bulk",
		"ess_tail"
	)
	drop_cols <- c(metric_cols, "Hazard.Ratio", "HR.Q2.5", "HR.Q97.5")

	collect_tables <- function(x, path = character()) {
		if (is.null(x)) {
			return(list())
		}
		if (is.data.frame(x)) {
			if (!any(c("Rhat", "ess_bulk", "ess_tail") %in% names(x))) {
				return(list())
			}
			id_cols <- setdiff(names(x), drop_cols)
			out <- x[, c(id_cols, intersect(metric_cols, names(x))), drop = FALSE]
			out$section <- if (length(path) >= 1L) path[[1L]] else NA_character_
			out$subsection <- if (length(path) > 1L) {
				paste(path[-1L], collapse = " / ")
			} else {
				NA_character_
			}
			out$parameter_label <- .diagnosis_parameter_label(out)
			keep_id_cols <- setdiff(
				id_cols,
				c("section", "subsection", "parameter_label")
			)
			out <- out[,
				c(
					"section",
					"subsection",
					"parameter_label",
					keep_id_cols,
					intersect(metric_cols, names(out))
				),
				drop = FALSE
			]
			return(list(out))
		}
		if (is.list(x)) {
			nm <- names(x)
			if (is.null(nm)) {
				nm <- paste0("item", seq_along(x))
			}
			out <- list()
			for (idx in seq_along(x)) {
				child_name <- nm[[idx]]
				if (is.null(child_name) || !nzchar(child_name)) {
					child_name <- paste0("item", idx)
				}
				out <- c(out, collect_tables(x[[idx]], c(path, child_name)))
			}
			return(out)
		}
		list()
	}

	tables <- tables[setdiff(names(tables), "diagnostics")]
	out <- collect_tables(tables)
	if (length(out) == 0L) {
		return(data.frame())
	}
	all_cols <- unique(unlist(lapply(out, names), use.names = FALSE))
	out <- lapply(out, function(df) {
		missing_cols <- setdiff(all_cols, names(df))
		if (length(missing_cols) > 0L) {
			for (col_name in missing_cols) {
				df[[col_name]] <- NA
			}
		}
		df[, all_cols, drop = FALSE]
	})
	out <- do.call(rbind, out)
	rownames(out) <- NULL
	out
}

#' @export
diagnosis.JoiNMeFit <- function(
	object,
	draws = NULL,
	seed = 1,
	digits = 3,
	include_corr = TRUE,
	...
) {
	
	cache_key <- paste0(
		"diagnosis_",
		draws %||% "default",
		"_",
		digits,
		"_",
		include_corr
	)
	cached <- object$cache_get(cache_key)
	if (!is.null(cached)) {
		return(cached)
	}
	sum_obj <- summary(
		object,
		draws = draws,
		seed = seed,
		digits = digits,
		include_corr = include_corr,
		...
	)
	result <- structure(
		list(
			summary = sum_obj$tables$diagnostics %||%
				.diagnostics_table_from_sampler(sum_obj$diagnostics %||% list()),
			by_parameter = .diagnosis_parameter_table_from_tables(sum_obj$tables),
			sampler = sum_obj$diagnostics %||% list(),
			metadata = c(
				sum_obj$metadata %||% list(),
				list(include_corr = include_corr, digits = digits)
			)
		),
		class = "JoiNMe_diagnosis"
	)
	object$cache_set(cache_key, result)
	result
}

#' @export
diagnosis.JoiNMeDynPred <- function(object, ...) {
	assertthat::assert_that(
		inherits(object, "JoiNMeDynPred"),
		msg = "Object must be a JoiNMeDynPred instance."
	)
	sum_obj <- summary(object)
	sum_obj$tables$diagnostics %||% .build_common_diagnostics_table()
}

#' @export
print.JoiNMe_diagnosis <- function(x, max_rows = 20, ...) {
	
	.cli_print_table_section(
		"Diagnostics Summary",
		x$summary,
		level = 1L,
		formatter = .format_common_diagnostics_for_print
	)
	if (is.data.frame(x$by_parameter) && nrow(x$by_parameter) > 0L) {
		.tbl <- x$by_parameter
		if (nrow(.tbl) > max_rows) {
			.tbl <- utils::head(.tbl, max_rows)
		}
		.cli_print_table_section("Per-Parameter Diagnostics", .tbl, level = 2L)
		if (nrow(x$by_parameter) > max_rows) {
			cat(
				"... truncated to ",
				max_rows,
				" rows; inspect $by_parameter for the full table.\n",
				sep = ""
			)
		}
	}
	invisible(x)
}

# ---- log-likelihood extraction -------------------------------------------

#' Log-likelihood summary
#'
#' @rdname log_lik.JoiNMeFit
#' @param object A fitted object of class `JoiNMeFit`.
#' @param what Character; which component to return: `"long"`, `"surv"`, or `"total"`.
#' @param draws Optional number of posterior draws to subset.
#' @param seed Random seed for draw subsetting.
#' @param ... Unused.
#'
#' @details The stored pointwise likelihood is evaluated by Stan with the full
#'   fitted association contribution. This includes the posterior ordinates of
#'   every ordered piecewise-linear transform; legacy `y` values are not used.
#'
#' @return A matrix
#' @seealso [log_lik()]
#' @export
log_lik.JoiNMeFit <- function(
	object,
	what = c("long", "surv", "total"),
	draws = NULL,
	seed = 1,
	...
) {
	
	what <- match.arg(what)
	sd <- object$stan_data

	vars_long <- paste0("log_lik_long[", seq_len(sd$N), "]")
	vars_surv <- paste0("log_lik_surv[", seq_len(sd$n_id), "]")
	vars <- switch(
		what,
		long = vars_long,
		surv = vars_surv,
		total = c(vars_long, vars_surv)
	)

	ddf <- .get_draws_df(object$fit, variables = vars, draws = draws, seed = seed)
	mat <- as.matrix(ddf[, vars, drop = FALSE])
	attr(mat, "component") <- what
	mat
}

#' Backup, unused
#' @noRd
#' @keywords internal
unused_loo.JoiNMeFit <- function(
	object,
	what = c("total", "long", "surv"),
	draws = NULL,
	seed = 1,
	...
) {
	what <- match.arg(what)
	ll_mat <- log_lik.JoiNMeFit(object, what = what, draws = draws, seed = seed)

	log_mean_exp <- function(x) {
		m <- max(x)
		m + log(mean(exp(x - m)))
	}

	log_colMeansExp <- function(mat) {
		S <- nrow(mat)
		matrixStats::colLogSumExps(mat, na.rm = TRUE) - log(S)
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
				.dimnames = list('elpd_loo', c("Estimate", "SE"))
			)
		),
		class = c('psis_loo', 'loo')
	)
	attr(output, 'what') <- what
	output
}

# ---- LOO / WAIC / ELPD ---------------------------------------------------

#' LOO-CV for JoiNMe models
#'
#' @rdname loo.JoiNMeFit
#' @param object A fitted object of class `JoiNMeFit`.
#' @param what Character; which component to return: `"long"`, `"surv"`, or `"total"`.
#' @param draws Optional number of posterior draws to subset.
#' @param seed Random seed for draw subsetting.
#' @param ... Additional arguments passed to `loo::loo()`.
#'
#' @return A `loo` object.
#' @importFrom loo loo
#' @export
loo.JoiNMeFit <- function(
	object,
	what = c("total", "long", "surv"),
	draws = NULL,
	seed = 1,
	...
) {
	what <- match.arg(what)
	ll_mat <- log_lik.JoiNMeFit(object, what = what, draws = draws, seed = seed)
	loo::loo(ll_mat, ...)
}

#' WAIC for JoiNMe models
#' @rdname waic.JoiNMeFit
#' @param object A fitted object of class `JoiNMeFit`.
#' @param what Character; which component to return: `"long"`, `"surv"`, or `"total"`.
#' @param draws Optional number of posterior draws to subset.
#' @param seed Random seed for draw subsetting.
#' @param ... Additional arguments passed to `loo::waic()`.
#'
#' @importFrom loo waic
#'
#' @return A `waic` object.
#' @export
waic.JoiNMeFit <- function(
	object,
	what = c("total", "long", "surv"),
	draws = NULL,
	seed = 1,
	...
) {
	if (!requireNamespace("loo", quietly = TRUE)) {
		cli::cli_abort(
			"Package {.pkg loo} is required for waic(). Install it with install.packages('loo')."
		)
	}
	what <- match.arg(what)
	ll_mat <- log_lik.JoiNMeFit(object, what = what, draws = draws, seed = seed)
	loo::waic(ll_mat, ...)
}


#' ELPD summary for JoiNMe models
#'
#' @param object A fitted object of class `JoiNMeFit`.
#' @param what Character; which component to return: `"long"`, `"surv"`, or `"total"`.
#' @param draws Optional number of posterior draws to subset.
#' @param seed Random seed for draw subsetting.
#' @param ... Additional arguments passed to `loo::elpd()`.
#'
#' @importFrom loo elpd
#' @rdname elpd.JoiNMeFit
#' @return A data frame with ELPD and standard error.
#' @export
elpd.JoiNMeFit <- function(
	object,
	what = c("total", "long", "surv"),
	draws = NULL,
	seed = 1,
	...
) {
	loo_obj <- loo.JoiNMeFit(object, what = what, draws = draws, seed = seed, ...)
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

#' Bayes factor between two JoiNMe fits
#'
#' @param fit1 First fitted `JoiNMeFit` object.
#' @param fit2 Second fitted `JoiNMeFit` object.
#' @param ... Additional arguments passed to `bridgesampling::bridge_sampler()`.
#'
#' @return A list with bridge sampling results and Bayes factor.
#' @export
bayes_factor <- function(fit1, fit2, ...) {
	assertthat::assert_that(
		inherits(fit1, "JoiNMeFit"),
		msg = "fit1 must be a JoiNMeFit instance."
	)
	assertthat::assert_that(
		inherits(fit2, "JoiNMeFit"),
		msg = "fit2 must be a JoiNMeFit instance."
	)

	if (!requireNamespace("bridgesampling", quietly = TRUE)) {
		cli::cli_abort(
			"Package {.pkg bridgesampling} is required for bayes_factor(). Install it with install.packages('bridgesampling')."
		)
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
#' @param object A fitted object of class `JoiNMeFit`.
#' @param newdataLong New longitudinal data frame for predictions.
#' @param newdataEvent New event data frame for predictions.
#' @param ci_level Numeric; credible interval level (default 0.95).
#' @param n_samples Integer; number of posterior samples to use (default 200).
#' @param ... Additional arguments.
#'
#' @importFrom bayesplot pp_check
#' @rdname pp_check.JoiNMeFit
#' @details
#' This crap is still under development.
#' @return A list with observation-level summaries and overall diagnostics.
#' @export
pp_check.JoiNMeFit <- function(
	object,
	newdataLong = NULL,
	newdataEvent = NULL,
	ci_level = 0.95,
	n_samples = 200,
	seed = 123,
	plot = FALSE,
	...
) {
	# Workflow: run posterior_epred -> summarise coverage and RMSE
	
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

	pred <- posterior_epred(
		object,
		newdataLong = newdataLong,
		newdataEvent = newdataEvent,
		control = list(n_samples = n_samples),
		seed = seed,
		...
	)

	draws_fit <- pred$draws$longitudinal_fitted
	if (is.null(draws_fit) || length(draws_fit) == 0) {
		cli::cli_abort(
			"No fitted draws found. Ensure longitudinal data is provided."
		)
	}

	id_var <- pred$metadata$id_var %||% (eval(object$call$id_var) %||% "id")
	time_var <- pred$metadata$time_var %||%
		(eval(object$call$time_var) %||% "time")
	marker_var <- pred$metadata$marker_var %||%
		(eval(object$call$marker_var) %||% "marker")
	y_var <- pred$metadata$response_var %||%
		.resolve_response_var(
			formulaLong = object$formulaLong,
			dataLong = newdataLong,
			context = "joinme_diagnosis()"
		)

	out_list <- list()
	for (id in names(draws_fit)) {
		# Per-subject summary of fitted draws vs observed outcomes
		dL <- newdataLong[newdataLong[[id_var]] == id, , drop = FALSE]
		if (nrow(dL) == 0) {
			next
		}

		epred_mat <- draws_fit[[id]]$epred
		if (is.null(epred_mat) || nrow(epred_mat) == 0) {
			next
		}

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
#' @noRd
.resolve_train_data <- function(
	object,
	newdataLong,
	newdataEvent,
	purpose = "analysis"
) {
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

#' Identify one finite value within a subject
#'
#' @description
#' Return the unique finite value in a vector, or `NA_real_` when no such value
#' exists or values vary.  This small helper permits `tapply()` to validate
#' interval-split time columns without defining functions inside the mapping
#' routine.
#'
#' @param x Numeric vector observed within one subject.
#'
#' @return One numeric value or `NA_real_`.
#' @keywords internal
#' @noRd
.constant_finite_value <- function(x) {
	x <- as.numeric(x)
	values <- unique(x[is.finite(x)])
	if (length(values) == 1L) values else NA_real_
}

#' Construct one dense dynamic-prediction time grid
#'
#' @param start Numeric landmark time.
#' @param horizon Numeric prediction horizon greater than `start`.
#'
#' @return A numeric vector of 50 equally spaced times including both bounds.
#' @keywords internal
#' @noRd
.discrimination_time_grid <- function(start, horizon) {
	seq(as.numeric(start), as.numeric(horizon), length.out = 50L)
}

#' Determine the default horizon for dynamic discrimination
#'
#' @description
#' Use the fitted time support when the caller omits both `time_horizon` and
#' `Dt`.  Fitted objects ordinarily retain `tmax`; the maximum observed event
#' or censoring time is a compatibility fallback for older objects.
#'
#' @param object A fitted `JoiNMeFit` object.
#' @param data_event Event-process data used when fitted time metadata are not
#'   available.
#'
#' @return One positive finite numeric horizon on the original time scale.
#' @keywords internal
#' @noRd
.default_discrimination_horizon <- function(object, data_event) {
	horizon <- object$tmax %||% object$config$tmax %||% object$stan_data$tmax
	if (
		is.numeric(horizon) &&
			length(horizon) == 1L &&
			is.finite(horizon) &&
			horizon > 0
	) {
		return(as.numeric(horizon))
	}
	event_vars <- .resolve_event_model_vars(
		formulaEvent = object$formulaEvent,
		dataEvent = data_event,
		context = "dynamic discrimination"
	)
	horizon <- max(as.numeric(event_vars$event_time), na.rm = TRUE)
	if (!is.finite(horizon) || horizon <= 0) {
		cli::cli_abort(
			"Could not determine a positive default prediction horizon from the fitted data."
		)
	}
	horizon
}

#' Validate a concordance time-weighting rule
#'
#' @param type_weights Character rule requested by the user.  `"none"` is the
#'   JoiNMe spelling of the unweighted `"n"` estimator.
#'
#' @return A character value accepted by the `timewt` argument of
#'   [survival::concordance()].
#' @keywords internal
#' @noRd
.concordance_time_weight <- function(type_weights) {
	valid_weights <- c("none", "n", "S", "S/G", "n/G2", "I")
	if (
		!is.character(type_weights) ||
			length(type_weights) != 1L ||
			!(type_weights %in% valid_weights)
	) {
		cli::cli_abort("{.arg type_weights} must be one of {.val {valid_weights}}.")
	}
	if (identical(type_weights, "none")) "n" else type_weights
}

#' Construct a subject-indexed landmark or horizon vector
#'
#' @description
#' Convert one of the time specifications accepted by the discrimination
#' methods into a finite numeric vector indexed by subject identifier.  A
#' scalar is repeated across subjects, a named vector is matched by subject,
#' and a character value identifies a column in the event data.  Keeping this
#' conversion in one helper ensures that concordance and AUC apply identical
#' subject alignment rules.
#'
#' @param value Numeric time specification or the name of a column in
#'   `data_event`.
#' @param ids Vector of distinct subject identifiers, in the required output
#'   order.
#' @param data_event Event-process data containing `id_var` and, when needed,
#'   the requested time column.
#' @param id_var Character name of the subject identifier column.
#' @param argument Character argument label used in diagnostic messages.
#'
#' @return A finite numeric vector with one element per subject and names equal
#'   to `as.character(ids)`.
#' @keywords internal
#' @noRd
.subject_time_map <- function(value, ids, data_event, id_var, argument) {
	id_keys <- as.character(ids)

	if (is.character(value)) {
		if (length(value) != 1L || is.na(value)) {
			cli::cli_abort(
				"{.arg {argument}} must name one column in {.arg newdataEvent}."
			)
		}
		if (!(value %in% names(data_event))) {
			cli::cli_abort(
				"Column {.val {value}} named by {.arg {argument}} was not found in {.arg newdataEvent}."
			)
		}
		if (!is.numeric(data_event[[value]])) {
			cli::cli_abort(
				"Column {.val {value}} named by {.arg {argument}} must be numeric."
			)
		}
		# Event data may contain several counting-process intervals per subject.
		# The first match is sufficient only when the time declaration is constant
		# within subject, so check that condition explicitly below.
		row_keys <- as.character(data_event[[id_var]])
		mapped <- as.numeric(data_event[[value]][match(id_keys, row_keys)])
		if (anyDuplicated(row_keys)) {
			group_range <- tapply(
				as.numeric(data_event[[value]]),
				row_keys,
				.constant_finite_value
			)
			mapped <- as.numeric(group_range[id_keys])
		}
	} else {
		if (!is.numeric(value) || !length(value)) {
			cli::cli_abort(
				"{.arg {argument}} must be numeric or name a column in {.arg newdataEvent}."
			)
		}
		if (length(value) == 1L) {
			mapped <- rep(as.numeric(value), length(ids))
		} else if (!is.null(names(value))) {
			if (anyDuplicated(names(value))) {
				cli::cli_abort(
					"Names in {.arg {argument}} must be unique subject identifiers."
				)
			}
			mapped <- as.numeric(value[match(id_keys, names(value))])
		} else if (length(value) == length(ids)) {
			mapped <- as.numeric(value)
		} else {
			cli::cli_abort(
				"{.arg {argument}} must be scalar, named by subject, or contain one value per subject."
			)
		}
	}

	if (length(mapped) != length(ids) || any(!is.finite(mapped))) {
		cli::cli_abort(c(
			x = "{.arg {argument}} could not be matched to a finite time for every subject.",
			i = "Check subject names and, for a column specification, use one constant value within each subject."
		))
	}
	stats::setNames(mapped, id_keys)
}

#' Reduce event-process rows to one observed outcome per subject
#'
#' @description
#' Resolve the survival response, order counting-process intervals by their
#' endpoints, and retain the final interval for each subject.  This produces
#' the subject-level failure or censoring outcome required by dynamic
#' discrimination estimators without multiplying subjects that have
#' interval-split covariates.
#'
#' @param formula_event Survival model formula fitted by `joinme()`.
#' @param data_event Event-process evaluation data.
#' @param id_var Character name of the subject identifier column.
#' @param context Character label used in diagnostic messages.
#'
#' @return A data frame with one row per subject and columns `event_time`,
#'   `event_status`, and `event_type`, in addition to the identifier column.
#' @keywords internal
#' @noRd
.subject_event_outcomes <- function(
	formula_event,
	data_event,
	id_var,
	context = "concordance.JoiNMeFit()"
) {
	if (!(id_var %in% names(data_event)) || anyNA(data_event[[id_var]])) {
		cli::cli_abort(
			"{context}: {.arg newdataEvent} must contain a non-missing subject identifier in {.field {id_var}}."
		)
	}
	event_vars <- .resolve_event_model_vars(
		formulaEvent = formula_event,
		dataEvent = data_event,
		context = context
	)
	if (!(event_vars$surv_type %in% c("right", "counting"))) {
		cli::cli_abort(c(
			x = "{context}: survival discrimination currently requires right-censored or counting-process outcomes.",
			i = "Left- and interval-censored outcomes require a censoring-specific discrimination estimator."
		))
	}

	outcome <- data.frame(
		.__id = data_event[[id_var]],
		.__event_start = as.numeric(event_vars$event_start),
		event_time = as.numeric(event_vars$event_time),
		event_status = event_vars$event_status,
		stringsAsFactors = FALSE
	)
	# In counting-process form the terminal interval contains the observed
	# subject outcome.  Ordering before de-duplication makes the result invariant
	# to the incoming row order.
	ord <- order(outcome$.__id, outcome$event_time, outcome$.__event_start)
	outcome <- outcome[ord, , drop = FALSE]
	outcome <- outcome[
		!duplicated(as.character(outcome$.__id), fromLast = TRUE),
		,
		drop = FALSE
	]
	decoded <- .derive_event_outcomes(outcome$event_status, context = context)
	outcome$event_status <- decoded$d_event
	outcome$event_type <- decoded$event_type
	outcome$.__event_start <- NULL
	names(outcome)[names(outcome) == ".__id"] <- id_var
	rownames(outcome) <- NULL
	outcome
}

#' Concordance of conditional survival predictions
#'
#' @description
#' Estimates a single, follow-up-wide concordance index from the conditional
#' survival curves of a fitted joint model.  The estimator follows the
#' survival-curve ordering of Antolini and colleagues: when subject \eqn{i}
#' experiences the event before subject \eqn{j}, their pair is concordant when
#' \eqn{\widehat S_i(r_i) < \widehat S_j(r_i)}, where \eqn{r_i} is the observed
#' residual event time and both curves are evaluated at that same residual
#' time.
#'
#' @param object A fitted object of class `JoiNMeFit`.
#' @param newdataLong Longitudinal evaluation data.  The fitted longitudinal
#'   data are used when this and `newdataEvent` are both omitted.
#' @param newdataEvent Event-process evaluation data.  The fitted event data
#'   are used when this and `newdataLong` are both omitted.
#' @param time_start Optional conditioning time.  A scalar applies to every
#'   subject; a named numeric vector is matched by subject identifier; and a
#'   character scalar names a subject-constant column in `newdataEvent`.  When
#'   omitted, each subject's latest longitudinal measurement strictly before
#'   their observed event or censoring time is used.
#' @param cause Positive integer identifying the event cause of interest.
#'   Other observed causes are treated as censoring at their event times.
#' @param predict_control Named list for prediction configuration.
#' @param seed Integer seed used when posterior draws are subsampled.
#' @param type_weights Character event-time weighting rule.  `"none"` and
#'   `"n"` give Harrell's pair weighting.  `"S"`, `"S/G"`, `"n/G2"` (Uno's
#'   weighting), and `"I"` have the meanings used by
#'   [survival::concordance()].
#' @param ... Additional arguments passed to [predict.JoiNMeFit()].
#'
#' @details
#' Residual follow-up is measured from the conditioning time:
#' \eqn{R_i=T_i-T_{0i}}.  Conditional survival is one at residual time zero.
#' For every observed event at \eqn{R_i}, the method compares subject \eqn{i}
#' with subjects known to have survived beyond \eqn{R_i}.  A subject censored
#' before \eqn{R_i} is not comparable; a subject censored exactly at
#' \eqn{R_i} is known to be event-free at that time and remains comparable.
#' Thus premature censoring is never interpreted as exceptionally long
#' survival.
#'
#' The censoring and event-time weights are obtained from
#' [survival::concordance()].  Under the usual independent-censoring
#' assumption, `"n/G2"` reduces sensitivity to the censoring distribution by
#' applying Uno's inverse-censoring weighting.  No method can recover the
#' ordering of a pair whose event-time order was not observed without further
#' assumptions.
#'
#' The comparison deliberately evaluates both predicted curves at the earlier
#' event time.  Evaluating each subject's curve at their own observed outcome
#' time and then ranking those values would use the response to construct the
#' predictor and can produce optimistically biased concordance.
#'
#' Dynamic predictions contain the complete fitted event model, including
#' baseline hazard, event covariates, latent longitudinal trajectories, and
#' every identity, functional, monotone-spline, or ordered piecewise-linear
#' association.  Consequently no separate piecewise-linear approximation is
#' made by this method.
#'
#' @return A one-row data frame with the concordance estimate, weighted counts
#'   of concordant, discordant, and predictor-tied pairs, total comparison
#'   weight `n_pairs`, and the numbers of analysed subjects and cause-specific
#'   events.
#' @references
#' Antolini L, Boracchi P, Biganzoli E (2005). A time-dependent discrimination
#' index for survival data. *Statistics in Medicine*, 24, 3927--3944.
#' \doi{10.1002/sim.2427}.
#' @seealso [auc.JoiNMeFit()], [predict.JoiNMeFit()], [survival::concordance()]
#' @importFrom survival Surv concordance
#' @export
concordance.JoiNMeFit <- function(
	object,
	newdataLong = NULL,
	newdataEvent = NULL,
	time_start = NULL,
	cause = 1,
	predict_control = list(n_samples = 200, n_pred_draws = 20),
	seed = .Random.seed[[1]],
	type_weights = "none",
	...
) {
  
  .warn_experimental("concordance")
  
	if (
		!is.numeric(cause) ||
			length(cause) != 1L ||
			!is.finite(cause) ||
			cause < 1 ||
			cause != as.integer(cause)
	) {
		cli::cli_abort("{.arg cause} must be a single positive integer.")
	}
	.concordance_time_weight(type_weights)

	data_in <- .resolve_train_data(
		object,
		newdataLong,
		newdataEvent,
		purpose = "conditional-survival concordance"
	)
	curves <- .concordance_survival_curves(
		object = object,
		newdataLong = data_in$newdataLong,
		newdataEvent = data_in$newdataEvent,
		time_start = time_start,
		cause = as.integer(cause),
		predict_control = predict_control,
		seed = seed,
		...
	)
	out <- .concordance_from_survival_curves(
		outcomes = curves$outcomes,
		survival = curves$survival,
		type_weights = type_weights
	)
	class(out) <- c("concordance_JoiNMeFit", "data.frame")
	out
}

#' Determine the final pre-outcome longitudinal measurement
#'
#' @description
#' Finds one conditioning time per event-record subject.  Measurements at the
#' event or censoring time are excluded so the prognostic history precedes the
#' outcome used to assess it.
#'
#' @param ids Subject identifiers in event-outcome order.
#' @param event_time Numeric observed event or censoring times aligned with
#'   `ids`.
#' @param data_long Longitudinal data containing identifier and time columns.
#' @param id_var,time_var Character column names.
#'
#' @return A named numeric vector aligned with `ids`; subjects without a
#'   pre-outcome measurement receive `NA_real_`.
#' @keywords internal
#' @noRd
.last_preoutcome_measurement <- function(
	ids,
	event_time,
	data_long,
	id_var,
	time_var
) {
	id_keys <- as.character(ids)
	event_map <- stats::setNames(as.numeric(event_time), id_keys)
	long_keys <- as.character(data_long[[id_var]])
	long_times <- as.numeric(data_long[[time_var]])
	subject_end <- as.numeric(event_map[long_keys])
	tolerance <- sqrt(.Machine$double.eps) * pmax(1, abs(subject_end))
	keep <- !is.na(long_keys) &
		is.finite(long_times) &
		is.finite(subject_end) &
		long_times < subject_end - tolerance
	if (!any(keep)) {
		return(stats::setNames(rep(NA_real_, length(ids)), id_keys))
	}

	candidates <- data.frame(
		id = long_keys[keep],
		time = long_times[keep],
		stringsAsFactors = FALSE
	)
	candidates <- candidates[
		order(candidates$id, candidates$time),
		,
		drop = FALSE
	]
	candidates <- candidates[
		!duplicated(candidates$id, fromLast = TRUE),
		,
		drop = FALSE
	]
	landmark <- candidates$time[match(id_keys, candidates$id)]
	stats::setNames(as.numeric(landmark), id_keys)
}

#' Build one exact residual-time prediction grid
#'
#' @param landmark Absolute conditioning time for one subject.
#' @param residual_end Observed residual event or censoring time for that
#'   subject.
#' @param event_residuals Distinct residual times of cause-specific events.
#'
#' @return A sorted numeric vector of absolute prediction times.  It contains
#'   every event residual time at which the subject can enter a comparison and
#'   a regular grid for stable curve interpolation.
#' @keywords internal
#' @noRd
.concordance_prediction_grid <- function(
	landmark,
	residual_end,
	event_residuals
) {
	relevant_events <- event_residuals[
		event_residuals >= 0 & event_residuals <= residual_end
	]
	regular_grid <- seq(0, residual_end, length.out = 50L)
	landmark + sort(unique(c(0, regular_grid, relevant_events, residual_end)))
}

#' Predict conditional survival curves required by concordance
#'
#' @description
#' Aligns terminal event outcomes, derives valid subject-specific landmarks,
#' restricts longitudinal histories to information available at those
#' landmarks, and predicts survival at every residual event time needed by a
#' comparable pair.
#'
#' @inheritParams concordance.JoiNMeFit
#'
#' @return A list with `outcomes`, a subject-level analysis data frame, and
#'   `survival`, a matrix whose rows are distinct residual event times and
#'   whose columns follow the subject order in `outcomes`.
#' @keywords internal
#' @noRd
.concordance_survival_curves <- function(
	object,
	newdataLong,
	newdataEvent,
	time_start,
	cause,
	predict_control,
	seed,
	...
) {
	id_var <- eval(object$call$id_var) %||% "id"
	time_var <- eval(object$call$time_var) %||% "time"
	if (
		!(id_var %in% names(newdataLong)) ||
			!(time_var %in% names(newdataLong))
	) {
		cli::cli_abort(
			"{.arg newdataLong} must contain the fitted identifier and time columns."
		)
	}

	outcomes <- .subject_event_outcomes(
		formula_event = object$formulaEvent,
		data_event = newdataEvent,
		id_var = id_var,
		context = "concordance.JoiNMeFit()"
	)
	ids <- outcomes[[id_var]]
	id_keys <- as.character(ids)
	if (is.null(time_start)) {
		landmark <- .last_preoutcome_measurement(
			ids = ids,
			event_time = outcomes$event_time,
			data_long = newdataLong,
			id_var = id_var,
			time_var = time_var
		)
	} else {
		landmark <- .subject_time_map(
			time_start,
			ids,
			newdataEvent,
			id_var,
			"time_start"
		)
	}

	eligible <- is.finite(landmark) &
		is.finite(outcomes$event_time) &
		outcomes$event_time > landmark
	if (sum(eligible) < 2L) {
		cli::cli_abort(c(
			x = "Fewer than two subjects have positive follow-up after a valid conditioning time.",
			i = "Supply earlier {.arg time_start} values or check the longitudinal observation times."
		))
	}
	outcomes <- outcomes[eligible, , drop = FALSE]
	landmark <- landmark[eligible]
	ids <- outcomes[[id_var]]
	id_keys <- as.character(ids)
	names(landmark) <- id_keys
	outcomes$time_start <- as.numeric(landmark)
	outcomes$residual_time <- outcomes$event_time - outcomes$time_start
	outcomes$cause_event <- as.integer(
		outcomes$event_status == 1L & outcomes$event_type == cause
	)
	if (!any(outcomes$cause_event == 1L)) {
		cli::cli_abort("No observed events match the requested {.arg cause}.")
	}

	# Only measurements available by the chosen landmark may inform the
	# conditional random-effects distribution.  This prevents measurements
	# recorded after the prognostic origin from leaking into the survival curve.
	long_keys <- as.character(newdataLong[[id_var]])
	long_landmark <- as.numeric(landmark[long_keys])
	keep_history <- is.finite(long_landmark) &
		is.finite(newdataLong[[time_var]]) &
		newdataLong[[time_var]] <= long_landmark
	history <- newdataLong[keep_history, , drop = FALSE]
	history_ids <- unique(as.character(history[[id_var]]))
	has_history <- id_keys %in% history_ids
	if (!all(has_history)) {
		outcomes <- outcomes[has_history, , drop = FALSE]
		landmark <- landmark[has_history]
		ids <- outcomes[[id_var]]
		id_keys <- as.character(ids)
		history <- history[
			as.character(history[[id_var]]) %in% id_keys,
			,
			drop = FALSE
		]
	}
	if (nrow(outcomes) < 2L || !any(outcomes$cause_event == 1L)) {
		cli::cli_abort(
			"Insufficient longitudinal histories remain for concordance estimation."
		)
	}

	event_residuals <- sort(unique(
		outcomes$residual_time[outcomes$cause_event == 1L]
	))
	times_map <- Map(
		.concordance_prediction_grid,
		as.numeric(outcomes$time_start),
		as.numeric(outcomes$residual_time),
		MoreArgs = list(event_residuals = event_residuals)
	)
	names(times_map) <- id_keys
	landmark <- stats::setNames(as.numeric(outcomes$time_start), id_keys)
	event_keys <- as.character(newdataEvent[[id_var]])
	prediction_event <- newdataEvent[event_keys %in% id_keys, , drop = FALSE]

	prediction <- predict(
		object,
		newdataLong = history,
		newdataEvent = prediction_event,
		process = "event",
		times = times_map,
		time_start = landmark,
		control = predict_control,
		seed = seed,
		...
	)
	survival_df <- prediction$predictions$survival
	if (is.null(survival_df) || nrow(survival_df) == 0L) {
		cli::cli_abort(
			"No conditional survival predictions were returned for concordance."
		)
	}
	survival_matrix <- .concordance_survival_matrix(
		survival_df = survival_df,
		ids = ids,
		landmark = landmark,
		event_residuals = event_residuals
	)
	list(outcomes = outcomes, survival = survival_matrix)
}

#' Align predicted curves on residual event times
#'
#' @description
#' Interpolates each posterior mean survival curve only within its predicted
#' support.  Exact event times were inserted in the prediction grid, so
#' interpolation ordinarily reproduces a stored ordinate; it also protects
#' against harmless floating-point changes during time scaling.
#'
#' @param survival_df Prediction summary with columns `id`, `time`, and
#'   `Survival`.
#' @param ids Subject identifiers in analysis order.
#' @param landmark Named absolute conditioning times.
#' @param event_residuals Sorted distinct cause-specific residual event times.
#'
#' @return Numeric matrix indexed by residual event time and subject.
#' @keywords internal
#' @noRd
.concordance_survival_matrix <- function(
	survival_df,
	ids,
	landmark,
	event_residuals
) {
	required <- c("id", "time", "Survival")
	if (!all(required %in% names(survival_df))) {
		cli::cli_abort("Conditional survival predictions lack required columns.")
	}
	id_keys <- as.character(ids)
	out <- matrix(
		NA_real_,
		nrow = length(event_residuals),
		ncol = length(ids),
		dimnames = list(as.character(event_residuals), id_keys)
	)
	prediction_keys <- as.character(survival_df$id)
	for (j in seq_along(ids)) {
		rows <- prediction_keys == id_keys[[j]] &
			is.finite(survival_df$time) &
			is.finite(survival_df$Survival)
		if (!any(rows)) {
			next
		}
		residual_grid <- as.numeric(survival_df$time[rows]) -
			as.numeric(landmark[[id_keys[[j]]]])
		ord <- order(residual_grid)
		out[, j] <- stats::approx(
			x = residual_grid[ord],
			y = as.numeric(survival_df$Survival[rows])[ord],
			xout = event_residuals,
			method = "linear",
			rule = 1,
			ties = mean
		)$y
	}
	out
}

#' Recover event-specific comparison weights from `survival`
#'
#' @description
#' Calls [survival::concordance()] with an arbitrary fixed ordering solely to
#' obtain its documented event-time weights.  Dividing the returned time
#' weight by the risk-set size gives the weight of each pair formed by that
#' event.  The arbitrary predictor cannot affect these weights.
#'
#' @param outcomes Subject-level data containing `residual_time` and
#'   `cause_event`.
#' @param type_weights JoiNMe concordance weighting declaration.
#'
#' @return Numeric vector aligned with rows of `outcomes`; non-event rows have
#'   zero weight.
#' @keywords internal
#' @noRd
.concordance_event_weights <- function(outcomes, type_weights) {
	timewt <- .concordance_time_weight(type_weights)
	weight_data <- data.frame(
		residual_time = outcomes$residual_time,
		cause_event = outcomes$cause_event,
		working_score = seq_len(nrow(outcomes))
	)
	rownames(weight_data) <- as.character(seq_len(nrow(weight_data)))
	engine <- survival::concordance(
		survival::Surv(residual_time, cause_event) ~ working_score,
		data = weight_data,
		timewt = timewt,
		ranks = TRUE
	)
	weights <- numeric(nrow(outcomes))
	if (is.null(engine$ranks) || nrow(engine$ranks) == 0L) {
		return(weights)
	}
	ranks <- engine$ranks
	# survival currently simplifies a single event's rank table to one column
	# with statistics in the rows.  Restore the ordinary one-row layout before
	# aligning the event.  This also keeps the method stable across survival
	# releases that preserve or drop this dimension.
	required_rank_columns <- c("time", "rank", "timewt", "casewt")
	if (
		!all(required_rank_columns %in% names(ranks)) &&
			all(required_rank_columns %in% rownames(ranks))
	) {
		ranks <- as.data.frame(t(as.matrix(ranks)), stringsAsFactors = FALSE)
	}
	if (!all(c("time", "timewt") %in% names(ranks))) {
		cli::cli_abort(
			"The concordance weighting engine returned an unrecognised rank table."
		)
	}

	event_rows <- suppressWarnings(as.integer(rownames(ranks)))
	valid_rows <- is.finite(event_rows) &
		event_rows >= 1L &
		event_rows <= nrow(outcomes)
	if (!all(valid_rows)) {
		# When row names are unavailable, match only observed event rows.  A
		# censoring observation may share the same time and must not receive the
		# event's comparison weight.
		event_candidates <- which(outcomes$cause_event == 1L)
		event_rows <- event_candidates[
			match(ranks$time, outcomes$residual_time[event_candidates])
		]
		valid_rows <- !is.na(event_rows) & is.finite(event_rows)
	}
	risk_size <- vapply(
		ranks$time,
		function(event_time) sum(outcomes$residual_time >= event_time),
		numeric(1)
	)
	pair_weight <- ranks$timewt / risk_size
	weights[event_rows[valid_rows]] <- pair_weight[valid_rows]
	weights
}

#' Calculate survival-curve concordance over all observable pairs
#'
#' @description
#' At each cause-specific event, compares the event subject's survival
#' probability with every subject whose residual failure time is known to be
#' later.  The implementation stores survival only on distinct event times and
#' visits each observable pair once, avoiding a three-dimensional subject by
#' subject by time array.
#'
#' @param outcomes Subject-level outcome data returned by
#'   `.concordance_survival_curves()`.
#' @param survival Matrix returned by `.concordance_survival_matrix()`.
#' @param type_weights Concordance event-time weighting rule.
#'
#' @return One-row data frame containing the concordance estimate and weighted
#'   comparison counts.
#' @keywords internal
#' @noRd
.concordance_from_survival_curves <- function(
	outcomes,
	survival,
	type_weights = "none"
) {
	required <- c("residual_time", "cause_event")
	if (!all(required %in% names(outcomes))) {
		cli::cli_abort("Concordance outcomes lack residual time or event status.")
	}
	event_weights <- .concordance_event_weights(outcomes, type_weights)
	concordant <- 0
	discordant <- 0
	tied <- 0

	event_rows <- which(outcomes$cause_event == 1L & event_weights > 0)
	for (i in event_rows) {
		event_time <- outcomes$residual_time[[i]]
		time_row <- match(as.character(event_time), rownames(survival))
		if (is.na(time_row) || !is.finite(survival[time_row, i])) {
			next
		}
		# Censoring at the event time implies survival strictly beyond that
		# time.  Failures tied at the same time do not have an observable order.
		comparator <- outcomes$residual_time > event_time |
			(outcomes$residual_time == event_time &
				outcomes$cause_event == 0L)
		comparator[[i]] <- FALSE
		comparator <- comparator & is.finite(survival[time_row, ])
		if (!any(comparator)) {
			next
		}

		difference <- survival[time_row, comparator] - survival[time_row, i]
		weight <- event_weights[[i]]
		concordant <- concordant + weight * sum(difference > 0)
		discordant <- discordant + weight * sum(difference < 0)
		tied <- tied + weight * sum(difference == 0)
	}

	n_pairs <- concordant + discordant + tied
	estimate <- if (n_pairs > 0) {
		(concordant + 0.5 * tied) / n_pairs
	} else {
		NA_real_
	}
	data.frame(
		concordance = estimate,
		concordant = concordant,
		discordant = discordant,
		tied = tied,
		n_pairs = n_pairs,
		n_subjects = nrow(outcomes),
		n_events = sum(outcomes$cause_event == 1L)
	)
}

#' Build a landmark-specific dynamic discrimination risk set
#'
#' @description
#' Obtain posterior mean conditional survival probabilities from the fitted
#' joint model and align them with one terminal event outcome per subject.  The
#' resulting data are the statistical input for cumulative/dynamic AUC.
#'
#' @param object A fitted `JoiNMeFit` object.
#' @param newdataLong Longitudinal histories used for dynamic prediction.
#' @param newdataEvent Event-process outcomes and covariates.
#' @param time_start Scalar or named subject-specific landmark times.
#' @param time_horizon Scalar or named subject-specific horizon times.
#' @param cause Integer event cause treated as the failure of interest.
#' @param n_samples Number of posterior draws used by dynamic prediction.
#' @param seed Integer posterior-subsampling seed.
#' @param ... Additional arguments passed to [predict.JoiNMeFit()].
#'
#' @return A subject-level data frame containing conditional event risk,
#'   follow-up from the landmark, the cause-specific event indicator, and the
#'   observed event or censoring outcome.
#' @keywords internal
#' @noRd
.dynamic_discrimination_risk_set <- function(
	object,
	newdataLong,
	newdataEvent,
	time_start,
	time_horizon,
	cause,
	n_samples,
	seed,
	...
) {
	id_var <- eval(object$call$id_var) %||% "id"
	time_var <- eval(object$call$time_var) %||% "time"
	event_df <- .subject_event_outcomes(
		formula_event = object$formulaEvent,
		data_event = newdataEvent,
		id_var = id_var,
		context = "auc.JoiNMeFit()"
	)
	ids <- event_df[[id_var]]
	id_keys <- as.character(ids)
	time_start_map <- .subject_time_map(
		time_start,
		ids,
		newdataEvent,
		id_var,
		"time_start"
	)
	time_horizon_map <- .subject_time_map(
		time_horizon,
		ids,
		newdataEvent,
		id_var,
		"time_horizon"
	)
	if (any(time_horizon_map <= time_start_map)) {
		cli::cli_abort(
			"{.arg time_horizon} must be greater than {.arg time_start} for every subject."
		)
	}

	# A direct indexed comparison avoids constructing and recombining one data
	# frame per subject.  This is material for repeated landmark calculations,
	# while preserving every longitudinal record observed by its landmark.
	long_keys <- as.character(newdataLong[[id_var]])
	long_landmark <- as.numeric(time_start_map[long_keys])
	keep_history <- is.finite(long_landmark) &
		is.finite(newdataLong[[time_var]]) &
		newdataLong[[time_var]] <= long_landmark
	newdataLong <- newdataLong[keep_history, , drop = FALSE]
	if (nrow(newdataLong) == 0L) {
		cli::cli_abort(
			"No longitudinal history is available at or before the landmark time(s)."
		)
	}

	# Fifty points preserve the established prediction behaviour and include the
	# requested horizon exactly.  The same dense survival grid is therefore used
	# by discrimination summaries and ordinary dynamic prediction plots.
	times_map <- Map(
		.discrimination_time_grid,
		as.numeric(time_start_map[id_keys]),
		as.numeric(time_horizon_map[id_keys])
	)
	names(times_map) <- id_keys

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
	if (is.null(surv_df) || nrow(surv_df) == 0L) {
		cli::cli_abort(
			"No survival predictions were returned for time-varying discrimination."
		)
	}

	# Match each prediction row to its requested subject-specific horizon.  The
	# relative tolerance protects against harmless floating-point rounding in a
	# generated time grid without admitting a neighbouring grid point.
	pred_keys <- as.character(surv_df$id)
	pred_horizon <- as.numeric(time_horizon_map[pred_keys])
	tolerance <- sqrt(.Machine$double.eps) * pmax(1, abs(pred_horizon))
	at_horizon <- is.finite(pred_horizon) &
		is.finite(surv_df$time) &
		abs(surv_df$time - pred_horizon) <= tolerance
	surv_horizon <- surv_df[at_horizon, , drop = FALSE]
	pred_index <- match(id_keys, as.character(surv_horizon$id))
	if (all(is.na(pred_index))) {
		cli::cli_abort(
			"Failed to align survival predictions at the requested horizon."
		)
	}

	risk_set <- event_df
	risk_set$time_start <- as.numeric(time_start_map[id_keys])
	risk_set$time_horizon <- as.numeric(time_horizon_map[id_keys])
	risk_set$risk <- 1 - as.numeric(surv_horizon$Survival[pred_index])
	# A competing event is a censoring observation for a cause-specific C-index.
	# Censoring occurs at its observed time, rather than at the prediction
	# horizon, so it contributes only comparisons that were still observable.
	cause_event <- risk_set$event_status == 1L & risk_set$event_type == cause
	risk_set$event_window <- as.integer(
		cause_event & risk_set$event_time <= risk_set$time_horizon
	)
	risk_set$time_window <- pmin(risk_set$event_time, risk_set$time_horizon) -
		risk_set$time_start
	# Dynamic prediction is conditional on survival beyond the landmark.  Remove
	# subjects whose outcome was already known at or before that time.
	risk_set <- risk_set[
		is.finite(risk_set$event_time) & risk_set$event_time > risk_set$time_start,
		,
		drop = FALSE
	]
	rownames(risk_set) <- NULL
	risk_set
}

# ---- Time-dependent AUC --------------------------------------------------

#' Time-dependent area under the ROC curve
#'
#' @description
#' Estimate the cumulative/dynamic area under the receiver operating
#' characteristic curve for the survival process at one or more landmark
#' times.  Cases experience the requested cause after the landmark and by the
#' horizon; controls remain event-free beyond the horizon.  Subjects censored
#' before the horizon are excluded because their case/control state is unknown.
#'
#' Posterior mean dynamic risks are obtained through `predict.JoiNMeFit()`, so
#' the calculation automatically uses the complete fitted hazard: baseline
#' hazard, event covariates, marker weights, and any identity, functional,
#' monotone-spline, or ordered piecewise-linear association transform.
#'
#' @param object A fitted object.  Methods are currently provided for
#'   `JoiNMeFit`.
#' @param ... Arguments passed to a class-specific method.
#'
#' @return For a `JoiNMeFit`, a data frame with landmark and horizon times,
#'   cumulative/dynamic AUC, and the numbers of cases, controls, and comparable
#'   case-control pairs.
#' @export
auc <- function(object, ...) {
	UseMethod("auc")
}

#' Calculate cumulative/dynamic AUC from an aligned risk set
#'
#' @description
#' Compare posterior mean risks for cases observed by the horizon with risks
#' for controls known to remain event-free beyond the horizon.  The
#' Mann--Whitney rank identity gives half credit to tied risks and avoids
#' allocating the full case-by-control comparison matrix.
#'
#' @param risk_set Subject-level output from
#'   `.dynamic_discrimination_risk_set()`.
#' @param time_start Numeric landmark used to label the result.
#' @param time_horizon Numeric horizon used to label the result.
#'
#' @return A one-row data frame containing AUC and its case, control, and pair
#'   counts.
#' @keywords internal
#' @noRd
.auc_from_risk_set <- function(risk_set, time_start, time_horizon) {
	required <- c("risk", "event_window", "event_time", "time_horizon")
	if (!all(required %in% names(risk_set))) {
		return(data.frame(
			time_start = time_start,
			time_horizon = time_horizon,
			auc = NA_real_,
			n_cases = 0L,
			n_controls = 0L,
			n_pairs = 0L
		))
	}

	case_risk <- as.numeric(risk_set$risk[risk_set$event_window == 1L])
	control_risk <- as.numeric(risk_set$risk[
		risk_set$event_window == 0L & risk_set$event_time > risk_set$time_horizon
	])
	case_risk <- case_risk[is.finite(case_risk)]
	control_risk <- control_risk[is.finite(control_risk)]
	n_cases <- length(case_risk)
	n_controls <- length(control_risk)
	n_pairs <- n_cases * n_controls
	auc_value <- NA_real_
	if (n_pairs > 0L) {
		# The first n_cases ranks belong to cases.  Subtracting their minimum
		# possible rank sum leaves the number of control risks below a case risk,
		# with average ranks contributing one half for a tied case-control pair.
		combined_risk <- c(case_risk, control_risk)
		case_rank_sum <- sum(
			rank(combined_risk, ties.method = "average")[seq_len(n_cases)]
		)
		auc_value <- (case_rank_sum - n_cases * (n_cases + 1) / 2) / n_pairs
	}
	data.frame(
		time_start = time_start,
		time_horizon = time_horizon,
		auc = auc_value,
		n_cases = n_cases,
		n_controls = n_controls,
		n_pairs = n_pairs
	)
}

#' @rdname auc
#' @param newdataLong Longitudinal evaluation data; defaults to the fitted data.
#' @param newdataEvent Event evaluation data; defaults to the fitted data.
#' @param time_start Numeric landmark time or vector of landmark times.  The
#'   default is zero.
#' @param time_horizon Numeric horizon time.  It may be scalar or have the same
#'   length as `time_start`.  When both `time_horizon` and `Dt` are omitted,
#'   the fitted maximum time is used.
#' @param Dt Positive horizon width used when `time_horizon` is omitted.
#' @param cause Integer competing-risk cause, with one denoting the primary
#'   event type.
#' @param n_samples Number of posterior draws used for dynamic prediction.
#' @param seed Integer seed used when posterior draws are subsampled.
#'
#' @details
#' For case risks \eqn{r_i} and control risks \eqn{r_j}, the estimator is the
#' proportion of comparable pairs satisfying \eqn{r_i > r_j}, with half credit
#' for ties.  This is the empirical cumulative/dynamic AUC.  It does not apply
#' inverse-probability-of-censoring weights; early-censored subjects are omitted
#' explicitly and the returned counts make the resulting comparison set clear.
#'
#' @export
auc.JoiNMeFit <- function(
	object,
	newdataLong = NULL,
	newdataEvent = NULL,
	time_start = 0,
	time_horizon = NULL,
	Dt = NULL,
	cause = 1,
	n_samples = 200,
	seed = 123,
	...
) {
	.warn_experimental("auc")
	if (
		!is.numeric(time_start) ||
			!length(time_start) ||
			any(!is.finite(time_start))
	) {
		cli::cli_abort(
			"{.arg time_start} must be a finite numeric value or vector."
		)
	}
	if (
		!is.numeric(cause) ||
			length(cause) != 1L ||
			!is.finite(cause) ||
			cause < 1 ||
			cause != as.integer(cause)
	) {
		cli::cli_abort("{.arg cause} must be a single positive integer.")
	}
	data_in <- .resolve_train_data(
		object,
		newdataLong,
		newdataEvent,
		purpose = "time-dependent AUC"
	)
	if (is.null(time_horizon)) {
		if (is.null(Dt)) {
			time_horizon <- rep(
				.default_discrimination_horizon(object, data_in$newdataEvent),
				length(time_start)
			)
		} else {
			if (!is.numeric(Dt) || length(Dt) != 1L || !is.finite(Dt) || Dt <= 0) {
				cli::cli_abort("{.arg Dt} must be a positive finite scalar.")
			}
			time_horizon <- time_start + Dt
		}
	} else {
		if (
			!is.numeric(time_horizon) ||
				any(!is.finite(time_horizon)) ||
				!(length(time_horizon) %in% c(1L, length(time_start)))
		) {
			cli::cli_abort(
				"{.arg time_horizon} must be finite and have length one or length(time_start)."
			)
		}
		if (length(time_horizon) == 1L) {
			time_horizon <- rep(time_horizon, length(time_start))
		}
	}
	if (any(time_horizon <= time_start)) {
		cli::cli_abort(
			"Every {.arg time_horizon} must be greater than its landmark time."
		)
	}

	rows <- lapply(seq_along(time_start), function(i) {
		risk_set <- .dynamic_discrimination_risk_set(
			object = object,
			newdataLong = data_in$newdataLong,
			newdataEvent = data_in$newdataEvent,
			time_start = time_start[i],
			time_horizon = time_horizon[i],
			cause = as.integer(cause),
			n_samples = n_samples,
			seed = seed,
			...
		)
		.auc_from_risk_set(
			risk_set = risk_set,
			time_start = time_start[i],
			time_horizon = time_horizon[i]
		)
	})
	out <- do.call(rbind, rows)
	class(out) <- c("tvAUC_JoiNMeFit", "data.frame")
	out
}

#' @rdname auc
#' @export
AUC <- auc

# ---- Stan diagnostics ----------------------------------------------------

#' Stan diagnostics for JoiNMe models
#'
#' @name stan_diagnostics.JoiNMeFit
#' @rdname stan_diagnostics.JoiNMeFit
#' @param object A fitted object of class `JoiNMeFit`.
#' @param pars Optional character vector of parameter names to include.
#' @param regex_pars Optional regular expression for parameter selection.
#' @param draws Optional number of posterior draws to subset.
#' @param seed Random seed for draw subsetting.
#' @param type Diagnostic type (for ESS and MCSE).
#' @param ... Unused.
#'
#' @return A data frame of diagnostics by parameter.
#' @export
stan_rhat.JoiNMeFit <- function(
	object,
	pars = NULL,
	regex_pars = NULL,
	draws = NULL,
	seed = 1,
	...
) {
	
	draws_obj <- draws(object, format = "draws_array")
	vars <- posterior::variables(draws_obj)
	vars <- .filter_diag_vars(vars, pars, regex_pars)
	if (length(vars) == 0) {
		return(data.frame(variable = character(0), rhat = numeric(0)))
	}
	draws_obj <- draws(
		object,
		variables = vars,
		draws = draws,
		seed = seed,
		format = "draws_array"
	)
	.diag_summary_df(draws_obj, metric = "rhat", vars = vars)
}

#' @rdname stan_diagnostics.JoiNMeFit
#' @export
stan_ess.JoiNMeFit <- function(
	object,
	pars = NULL,
	regex_pars = NULL,
	draws = NULL,
	seed = 1,
	type = c("bulk", "tail"),
	...
) {
	
	type <- match.arg(type)
	draws_obj <- draws(object, format = "draws_array")
	vars <- posterior::variables(draws_obj)
	vars <- .filter_diag_vars(vars, pars, regex_pars)
	if (length(vars) == 0) {
		return(data.frame(variable = character(0), ess = numeric(0)))
	}
	draws_obj <- draws(
		object,
		variables = vars,
		draws = draws,
		seed = seed,
		format = "draws_array"
	)
	metric_name <- if (type == "bulk") "ess_bulk" else "ess_tail"
	metric_df <- .diag_summary_df(draws_obj, metric = metric_name, vars = vars)
	names(metric_df)[2] <- "ess"
	metric_df
}

#' @rdname stan_diagnostics.JoiNMeFit
#' @export
stan_mcse.JoiNMeFit <- function(
	object,
	pars = NULL,
	regex_pars = NULL,
	draws = NULL,
	seed = 1,
	type = c("mean", "sd", "median"),
	...
) {
	
	type <- match.arg(type)
	draws_obj <- draws(object, format = "draws_array")
	vars <- posterior::variables(draws_obj)
	vars <- .filter_diag_vars(vars, pars, regex_pars)
	if (length(vars) == 0) {
		return(data.frame(variable = character(0), mcse = numeric(0)))
	}
	draws_obj <- draws(
		object,
		variables = vars,
		draws = draws,
		seed = seed,
		format = "draws_array"
	)
	metric_name <- switch(
		type,
		mean = "mcse_mean",
		sd = "mcse_sd",
		median = "mcse_median"
	)
	metric_df <- .diag_summary_df(draws_obj, metric = metric_name, vars = vars)
	names(metric_df)[2] <- "mcse"
	metric_df
}

#' @keywords internal
#' @noRd

.diag_summary_df <- function(draws_obj, metric, vars) {
	metric_df <- suppressWarnings(tryCatch(
		posterior::summarise_draws(draws_obj, metric),
		error = function(e) NULL
	))
	if (is.null(metric_df) || !all(c("variable", metric) %in% names(metric_df))) {
		cli::cli_abort("Failed to summarise posterior diagnostics for plotting.")
	}
	metric_df <- as.data.frame(
		metric_df[, c("variable", metric), drop = FALSE],
		row.names = NULL,
		stringsAsFactors = FALSE
	)
	if (nrow(metric_df) == 0L) {
		out <- data.frame(
			variable = character(0),
			value = numeric(0),
			row.names = NULL
		)
		names(out)[2] <- metric
		return(out)
	}
	if (!is.null(vars)) {
		metric_df <- metric_df[match(vars, metric_df$variable), , drop = FALSE]
		if (nrow(metric_df) != length(vars) || anyNA(metric_df$variable)) {
			cli::cli_abort(
				"Posterior diagnostic output does not align with selected parameters."
			)
		}
	}
	out <- data.frame(
		variable = metric_df$variable,
		value = metric_df[[metric]],
		row.names = NULL
	)
	names(out)[2] <- metric
	out
}

#' @keywords internal
#' @noRd
.filter_diag_vars <- function(vars, pars = NULL, regex_pars = NULL) {
	if (!is.null(pars)) {
		vars <- intersect(vars, pars)
	}
	if (!is.null(regex_pars)) {
		vars <- vars[grepl(regex_pars, vars)]
	}
	vars
}
