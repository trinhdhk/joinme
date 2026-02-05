#' @name joinme_diagnosis
#' @title Model Diagnosis
#'
#' @description
#' Posterior predictive checks and information criteria.
#'
#' @keywords internal
NULL

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
pp_check.JoinMeFit <- function(object, newdataLong, newdataEvent, ci_level = 0.95, n_samples = 200, seed = 123, plot = FALSE,...) {
	assertthat::assert_that(inherits(object, "JoinMeFit"), msg = "Object must be a JoinMeFit instance.")
	assertthat::assert_that(!missing(newdataLong), msg = "newdataLong is required.")
	assertthat::assert_that(!missing(newdataEvent), msg = "newdataEvent is required.")

	ci_level <- .validate_ci_levels(ci_level)
	probs <- .quantile_probs_from_ci(ci_level)

	pred <- posterior_epred(object,
				newdataLong = newdataLong,
				newdataEvent = newdataEvent,
				n_samples = n_samples,
				seed = seed,
				...)

	draws_fit <- pred$draws$longitudinal_fitted
	if (is.null(draws_fit) || length(draws_fit) == 0) {
		cli::cli_abort("No fitted draws found. Ensure longitudinal data is provided.")
	}

	id_var <- pred$metadata$id_var %||% (eval(object$call$id_var) %||% "id")
	time_var <- pred$metadata$time_var %||% (eval(object$call$time_var) %||% "time")
	marker_var <- pred$metadata$marker_var %||% (eval(object$call$marker_var) %||% "marker")
	y_var <- pred$metadata$response_var %||% (eval(object$call$y_var) %||% "y")

	out_list <- list()
	for (id in names(draws_fit)) {
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

# ---- Concordance ----------------------------------------------------------

#' Concordance index for survival component
#'
#' @rdname concordance.JoinMeFit
#' @param object A fitted object of class `JoinMeFit`.
#' @param cause Integer; cause index for competing risks (default 1).
#' @param draws Optional number of posterior draws to subset.
#' @param seed Random seed for draw subsetting.
#' @param ... Additional arguments passed to `survival::concordance()`.
#'
#' @importFrom survival Surv concordance
#' @return A `concordance` object.
#' @export
concordance.JoinMeFit <- function(object, cause = 1, draws = NULL, seed = 1, ...) {
	sd <- object$stan_data
	if (is.null(sd$W) || sd$p_w < 1) {
		cli::cli_abort("No baseline survival covariates found for concordance calculation.")
	}

	if (!is.numeric(cause) || length(cause) != 1 || cause < 1 || cause > (sd$K_event %||% 1)) {
		cli::cli_abort("cause must be a single integer between 1 and K_event.")
	}

	vars <- paste0("gamma_w[", cause, ",", seq_len(sd$p_w), "]")
	ddf <- .get_draws_df(object$fit, variables = vars, draws = draws, seed = seed)
	gamma_mean <- colMeans(as.matrix(ddf[, vars, drop = FALSE]))

	concord_data <- data.frame(
		time = sd$S_event,
		status = if (!is.null(sd$event_type)) {
			ifelse(sd$d_event == 1 & sd$event_type == cause, 1, 0)
		} else {
			sd$d_event
		},
		risk = as.numeric(sd$W %*% gamma_mean)
	)

	survival::concordance(survival::Surv(time, status) ~ risk, data = concord_data, ...)
}
