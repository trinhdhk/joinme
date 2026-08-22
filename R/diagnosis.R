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
		"Mean",
		"Median",
		"Est.Error",
		"SD",
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
			out$family <- if (length(path) >= 1L) path[[1L]] else NA_character_
			# out$subsection <- if (length(path) > 1L) {
			# 	paste(path[-1L], collapse = " / ")
			# } else {
			# 	NA_character_
			# }
			out$parameter <- .diagnosis_parameter_label(out)
			keep_id_cols <- setdiff(
				id_cols,
				c("family", 
				# "subsection", 
				"parameter")
			)
			out <- out[,
				c(
					"family",
					# "subsection",
					"parameter",
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
		"_seed=",
		seed,
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

#' Build posterior features used to interpret Hamiltonian energy movement
#'
#' @description
#' Reduces large random-effect and scale blocks to one value per iteration and
#' chain.  Correlation with energy is an exploratory diagnostic: it identifies
#' posterior directions along which energy changes, but it does not establish
#' that any one parameter caused a low E-BFMI value.
#'
#' @param object A fitted latent-progress mixture.
#'
#' @return A named list of iteration-by-chain matrices.
#' @keywords internal
#' @noRd
.mixture_energy_features <- function(object) {
	all_variables <- tryCatch(
		posterior::variables(.get_draws_obj(object$fit)),
		error = function(error) character(0)
	) # complete fitted parameter names used to locate each interpretable block
	if (length(all_variables) == 0L) {
		return(list())
	}

	feature_from_pattern <- function(pattern, summary_function) {
		selected_variables <- grep(
			pattern,
			all_variables,
			value = TRUE
		) # fitted variables belonging to one scale or latent-effect block
		if (length(selected_variables) == 0L) {
			return(NULL)
		}
		draw_array <- as.array(.get_draws_array(
			object$fit,
			variables = selected_variables
		)) # iteration-by-chain draws for the selected posterior block
		summary_function(draw_array)
	}

	rms_summary <- function(draw_array) {
		sqrt(apply(draw_array^2, c(1L, 2L), mean, na.rm = TRUE))
	}
	minimum_summary <- function(draw_array) {
		apply(draw_array, c(1L, 2L), min, na.rm = TRUE)
	}
	geometric_mean_summary <- function(draw_array) {
		exp(apply(
			log(pmax(abs(draw_array), sqrt(.Machine$double.eps))),
			c(1L, 2L),
			mean,
			na.rm = TRUE
		))
	}

	features <- list(
		"within-class scale RMS" = feature_from_pattern(
			"^mix_scale\\[",
			rms_summary
		),
		"smallest within-class scale" = feature_from_pattern(
			"^mix_scale\\[",
			minimum_summary
		),
		"class-location RMS" = feature_from_pattern(
			"^mix_location\\[",
			rms_summary
		),
		"subject scale geometric mean" = feature_from_pattern(
			"^tau_u\\[",
			geometric_mean_summary
		),
		"marker scale geometric mean" = feature_from_pattern(
			"^tau_v\\[",
			geometric_mean_summary
		),
		"covariance loading geometric mean" = feature_from_pattern(
			"^lambda_L\\[",
			geometric_mean_summary
		),
		"subject standardised-effect RMS" = feature_from_pattern(
			"^z_u\\[",
			rms_summary
		),
		"marker standardised-effect RMS" = feature_from_pattern(
			"^z_v\\[",
			rms_summary
		),
		"covariance standardised-effect RMS" = feature_from_pattern(
			"^z_L\\[",
			rms_summary
		),
		"class-regression coefficient RMS" = feature_from_pattern(
			"^mix_class_coefficient_",
			rms_summary
		)
	) # scientifically named low-dimensional summaries of the main mixture geometry
	Filter(Negate(is.null), features)
}

#' Relate posterior block summaries to Hamiltonian energy
#'
#' @param object A fitted latent-progress mixture.
#'
#' @return A long data frame with one feature-by-chain correlation per row.
#' @keywords internal
#' @noRd
.mixture_energy_association <- function(object) {
	sampler_draws <- .sampler_diagnostic_draws(
		object$fit
	) # sampler states in the same iteration and chain order as posterior draws
	if (is.null(sampler_draws)) {
		return(data.frame())
	}
	sampler_array <- as.array(
		sampler_draws
	) # numeric sampler-state array used to retrieve chain-wise energy
	if (!"energy__" %in% (dimnames(sampler_array)[[3L]] %||% character(0))) {
		return(data.frame())
	}
	features <- .mixture_energy_features(
		object
	) # posterior block summaries to compare with the sampled energy sequence
	if (length(features) == 0L) {
		return(data.frame())
	}

	rows <- list() # feature-by-chain diagnostic rows accumulated below
	row_position <- 1L # next free position in the diagnostic row list
	for (feature_name in names(features)) {
		feature_values <- features[[feature_name]] # iteration-by-chain feature matrix
		number_chains <- min(
			ncol(feature_values),
			dim(sampler_array)[2L]
		) # chains jointly available in model and sampler-state draws
		for (chain_index in seq_len(number_chains)) {
			number_iterations <- min(
				nrow(feature_values),
				dim(sampler_array)[1L]
			) # aligned retained iterations for this feature and energy
			feature_chain <- as.numeric(
				feature_values[seq_len(number_iterations), chain_index]
			) # posterior feature sequence for one chain
			energy_chain <- as.numeric(
				sampler_array[
					seq_len(number_iterations),
					chain_index,
					"energy__"
				]
			) # Hamiltonian energy sequence for the same chain and iterations
			finite_values <- is.finite(feature_chain) &
				is.finite(energy_chain) # paired iterations eligible for correlation
			correlation <- if (
				sum(finite_values) > 2L &&
					stats::sd(feature_chain[finite_values]) > 0 &&
					stats::sd(energy_chain[finite_values]) > 0
			) {
				stats::cor(
					feature_chain[finite_values],
					energy_chain[finite_values]
				)
			} else {
				NA_real_
			} # within-chain linear association, used only as an exploratory clue
			rows[[row_position]] <- data.frame(
				feature = feature_name,
				chain = chain_index,
				correlation_with_energy = correlation,
				# absolute_correlation = abs(correlation),
				stringsAsFactors = FALSE
			) # one interpretable energy association for this feature and chain
			row_position <- row_position + 1L
		}
	}
	output <- do.call(
		rbind,
		rows
	) # complete long-form energy association table
	# output[order(-output$correlation_with_energy), , drop = FALSE]
	dplyr::arrange(output, feature, chain)
}

#' Diagnose the likelihood scale trade-off in a mixture block
#'
#' @description
#' The selected standardised component scale is subsequently multiplied by an
#' ordinary random-effect scale.  This helper reports the within-chain
#' correlation between their logarithmic geometric means.  A strong negative
#' value is direct posterior evidence of the scale ridge discussed in the Stan
#' parameter declaration; it is not evidence of a coordinate-indexing error.
#'
#' @param object A fitted latent-progress mixture.
#'
#' @return A data frame with one row per active class type and chain.
#' @keywords internal
#' @noRd
.mixture_scale_tradeoff <- function(object) {
	stan_data <- object$stan_data # exact selected coordinates supplied to Stan
	number_classes <- as.integer(
		stan_data$n_classes %||% 0L
	) # shared number of mixture components
	if (number_classes < 2L) {
		return(data.frame())
	}

	block_specifications <- list(
		subject = list(
			dimension = as.integer(stan_data$mix_dim_subject %||% 0L),
			start = as.integer(stan_data$mix_start_subject %||% 0L),
			indices = as.integer(stan_data$mix_idx_subject %||% integer(0)),
			ordinary_prefix = "tau_u",
			likelihood_combination = "subject Cholesky scale x class coordinate"
		),
		marker = list(
			dimension = as.integer(stan_data$mix_dim_marker %||% 0L),
			start = as.integer(stan_data$mix_start_marker %||% 0L),
			indices = as.integer(stan_data$mix_idx_marker %||% integer(0)),
			ordinary_prefix = "tau_v",
			likelihood_combination = "marker Cholesky scale x class coordinate"
		),
		covariance = list(
			dimension = as.integer(stan_data$mix_dim_covariance %||% 0L),
			start = as.integer(stan_data$mix_start_covariance %||% 0L),
			indices = as.integer(stan_data$mix_idx_covariance %||% integer(0)),
			ordinary_prefix = "lambda_L",
			likelihood_combination = "covariance loading x class coordinate"
		)
	) # mapping from each active class block to its ordinary posterior scale

	rows <- list() # block-by-chain scale correlations accumulated below
	row_position <- 1L # next free position in the result list
	for (block_name in names(block_specifications)) {
		block <- block_specifications[[block_name]] # selected coordinate contract for one block
		if (
			block$dimension < 1L ||
				block$start < 1L ||
				length(block$indices) != block$dimension
		) {
			next
		}
		ordinary_variables <- paste0(
			block$ordinary_prefix,
			"[",
			block$indices,
			"]"
		) # ordinary scale variables paired with the selected latent coordinates
		packed_coordinates <- block$start +
			seq_len(block$dimension) - 1L # packed mixture columns belonging to this block
		component_variables <- as.vector(outer(
			seq_len(number_classes),
			packed_coordinates,
			function(class_index, coordinate_index) {
				paste0(
					"mix_scale[",
					class_index,
					",",
					coordinate_index,
					"]"
				)
			}
		)) # within-class scale variables for all classes in the current block
		ordinary_draws <- tryCatch(
			as.array(.get_draws_array(
				object$fit,
				variables = ordinary_variables
			)),
			error = function(error) NULL
		) # posterior draws of the ordinary scale factor
		component_draws <- tryCatch(
			as.array(.get_draws_array(
				object$fit,
				variables = component_variables
			)),
			error = function(error) NULL
		) # posterior draws of the class-conditional standardised scales
		if (is.null(ordinary_draws) || is.null(component_draws)) {
			next
		}
		log_ordinary_scale <- apply(
			log(pmax(ordinary_draws, sqrt(.Machine$double.eps))),
			c(1L, 2L),
			mean
		) # logarithmic geometric mean of ordinary scales per draw and chain
		log_component_scale <- apply(
			log(pmax(component_draws, sqrt(.Machine$double.eps))),
			c(1L, 2L),
			mean
		) # logarithmic geometric mean of component scales per draw and chain
		number_chains <- min(
			ncol(log_ordinary_scale),
			ncol(log_component_scale)
		) # jointly represented Markov chains
		for (chain_index in seq_len(number_chains)) {
			correlation <- stats::cor(
				log_ordinary_scale[, chain_index],
				log_component_scale[, chain_index],
				use = "complete.obs"
			) # posterior correlation revealing compensation between scale factors
			rows[[row_position]] <- data.frame(
				class_type = if (identical(block_name, "covariance")) {
					.mixture_covariance_class_type(
						object$mixture %||% object$config$mixture
					) %||% "covariance"
				} else {
					block_name
				},
				chain = chain_index,
				ordinary_scale = block$ordinary_prefix,
				likelihood_combination = block$likelihood_combination,
				correlation = correlation,
				# absolute_correlation = abs(correlation),
				# interpretation = if (
				# 	is.finite(correlation) && correlation <= -0.3
				# ) {
				# 	"material compensating scale movement"
				# } else {
				# 	"no strong negative linear trade-off detected"
				# },
				stringsAsFactors = FALSE
			) # transparent class-block scale assessment for one chain
			row_position <- row_position + 1L
		}
	}
	if (length(rows) == 0L) {
		return(data.frame())
	}
	do.call(rbind, rows)
}

#' Diagnostic summary for a latent-progress mixture
#'
#' @description
#' Extends the inherited diagnostics with chain-specific Hamiltonian energy
#' behaviour, exploratory associations between energy and model blocks, and a
#' direct assessment of the ordinary-scale/component-scale trade-off.
#'
#' @inheritParams diagnosis
#'
#' @return A `JoiNMeMix_diagnosis` object inheriting from
#'   `JoiNMe_diagnosis`.
#' @export
diagnosis.JoiNMeMixFit <- function(
	object,
	draws = NULL,
	seed = 1,
	digits = 3,
	include_corr = TRUE,
	...
) {
	cache_key <- paste0(
		"mixture_diagnosis_",
		draws %||% "default",
		"_seed=",
		seed,
		"_digits=",
		digits,
		"_corr=",
		as.integer(include_corr)
	) # cache identity for the expanded mixture diagnostic report
	cached <- object$cache_get(
		cache_key
	) # previously computed report for the same posterior subset and presentation
	if (!is.null(cached)) {
		return(cached)
	}

	result <- diagnosis.JoiNMeFit(
		object,
		draws = draws,
		seed = seed,
		digits = digits,
		include_corr = include_corr,
		...
	) # inherited sampler and parameter-level diagnostics
	backend_diagnostics <- result$sampler %||%
		list() # backend-wide values already collected by the inherited summary
	if (
		!is.data.frame(backend_diagnostics$energy_by_chain) ||
			nrow(backend_diagnostics$energy_by_chain) == 0L
	) {
		backend_diagnostics <- .joinme_sampler_diagnostics(
			object$fit
		) # fresh sampler states for an object whose earlier cache predates this schema
	}
	backend_table <- .diagnostics_table_from_sampler(
		backend_diagnostics
	) # common metric-value layout used by every JoiNMe summary
	if (is.data.frame(result$summary) && nrow(result$summary) > 0L) {
		for (metric_name in intersect(
			result$summary$metric,
			backend_table$metric
		)) {
			backend_value <- backend_table$value[
				match(metric_name, backend_table$metric)
			] # current backend-wide value for this established metric
			if (is.finite(backend_value)) {
				result$summary$value[
					match(metric_name, result$summary$metric)
				] <- backend_value
			}
		}
	}
	result$sampler <- backend_diagnostics
	result$energy_by_chain <- backend_diagnostics$energy_by_chain %||%
		data.frame()
	result$energy_association <- .mixture_energy_association(
		object
	) # ranked exploratory clues about which posterior blocks move with energy
	result$scale_tradeoff <- .mixture_scale_tradeoff(
		object
	) # direct posterior evidence for or against compensating scale movement
	result$metadata$mixture_geometry <- paste(
		"Class scales act before the ordinary random-effect transformation;",
		"the likelihood can therefore identify their combination more readily",
		"than its separate factors."
	) # concise interpretation retained with the diagnostic object
	class(result) <- c(
		"JoiNMeMix_diagnosis",
		"JoiNMe_diagnosis"
	) # specialised class which retains all inherited diagnosis behaviour
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
	if (
		is.data.frame(x$energy_by_chain) &&
			nrow(x$energy_by_chain) > 0L
	) {
		.cli_print_table_section(
			"Hamiltonian energy by chain",
			x$energy_by_chain,
			level = 2L
		)
		if (any(x$energy_by_chain$E_BFMI < 0.3, na.rm = TRUE)) {
			.cli_print_bullets(c(
				"E-BFMI below 0.3 indicates that at least one chain did not move efficiently through the posterior energy distribution.",
				"High energy persistence, long trajectories or a small step size support a posterior-geometry explanation; they do not by themselves identify a coding error."
			))
		}
	}
	if (
		is.data.frame(x$scale_tradeoff) &&
			nrow(x$scale_tradeoff) > 0L
	) {
		.cli_print_table_section(
			"Ordinary-scale and within-class-scale trade-off",
			x$scale_tradeoff,
			level = 2L
		)
	}
	if (
		is.data.frame(x$energy_association) &&
			nrow(x$energy_association) > 0L
	) {
		energy_table <- utils::head(
			x$energy_association,
			max_rows
		) # strongest absolute within-chain energy associations for concise printing
		# .cli_print_bullets(
		# 	"These correlations are exploratory indicators of posterior geometry, not causal tests or model-selection statistics."
		# )
		.cli_print_table_section(
			"Posterior features associated with energy",
			energy_table,
			level = 2L
		)
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
#'   every ordered piecewise-linear transform; earlier `y` values are not used.
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
		.get_response_var(
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
.get_train_data <- function(
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
	event_vars <- .get_event_model_vars(
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
#' conversion in one helper ensures that concordance and time-dependent ROC
#' analysis apply identical
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
	event_vars <- .get_event_model_vars(
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
#' @seealso [tvROC.JoiNMeFit()], [tvAUC.JoiNMeFit()],
#'   [predict.JoiNMeFit()], [survival::concordance()]
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
	.require_fitted_survival_process(
		object,
		"Concordance assessment"
	)
  
  .warn_experimental("concordance")

	# Concordance is evaluated over all observable event--comparator pairs and
	# consequently has no fixed prediction horizon.  Detect the former AUC
	# arguments here so that they cannot silently alter dynamic prediction.
	dot_names <- names(list(...))
	if (any(c("time_horizon", "Dt") %in% dot_names)) {
		cli::cli_abort(c(
			x = "{.fn concordance} does not use a prediction horizon.",
			i = "Use {.code tvROC(object, ...)} or {.code tvAUC(object, ...)} for landmark--horizon discrimination."
		))
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
	.concordance_time_weight(type_weights)

	data_in <- .get_train_data(
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

#' Extract horizon-specific survival draws for discrimination
#'
#' @description
#' Aligns the raw conditional-survival draws retained by
#' [predict.JoiNMeFit()] with one requested horizon for each analysed subject.
#'
#' @param prediction A `JoiNMeDynPred` object containing raw survival draws.
#' @param ids Subject identifiers in the required row order.
#' @param horizons Named numeric horizon vector indexed by subject identifier.
#'
#' @return A numeric subject-by-draw matrix of conditional event risks. Returns
#'   a zero-column matrix when raw survival draws are unavailable; the ROC
#'   calculation then uses posterior mean risks only.
#' @keywords internal
#' @noRd
.discrimination_horizon_risk_draws <- function(prediction, ids, horizons) {
	id_keys <- as.character(ids)
	draw_lists <- prediction$draws$survival %||% list()
	n_draws_by_subject <- integer(length(id_keys))
	for (subject_index in seq_along(id_keys)) {
		entry <- draw_lists[[id_keys[[subject_index]]]]
		if (!is.null(entry$matrix) && is.matrix(entry$matrix)) {
			n_draws_by_subject[[subject_index]] <- nrow(entry$matrix)
		}
	}
	available <- n_draws_by_subject > 0L
	if (!any(available)) {
		return(matrix(numeric(0), nrow = length(ids), ncol = 0L))
	}

	# Prediction ordinarily retains the same number of draws for every subject.
	# Taking the common minimum gives a rectangular matrix if an older object
	# contains unequal draw counts, without recycling or inventing draws.
	n_draws <- min(n_draws_by_subject[available])
	out <- matrix(
		NA_real_,
		nrow = length(ids),
		ncol = n_draws,
		dimnames = list(id_keys, paste0("draw_", seq_len(n_draws)))
	)
	for (subject_index in seq_along(id_keys)) {
		entry <- draw_lists[[id_keys[[subject_index]]]]
		if (is.null(entry$matrix) || !is.matrix(entry$matrix) || nrow(entry$matrix) < n_draws) {
			next
		}
		time_grid <- as.numeric(entry$time)
		horizon <- as.numeric(horizons[[id_keys[[subject_index]]]])
		tolerance <- sqrt(.Machine$double.eps) * max(1, abs(horizon))
		horizon_column <- which(
			is.finite(time_grid) & abs(time_grid - horizon) <= tolerance
		)[1L]
		if (is.na(horizon_column)) {
			next
		}
		out[subject_index, ] <- 1 - as.numeric(
			entry$matrix[seq_len(n_draws), horizon_column]
		)
	}
	out
}

#' Build a landmark-specific dynamic discrimination risk set
#'
#' @description
#' Obtain posterior mean conditional survival probabilities from the fitted
#' joint model and align them with one terminal event outcome per subject.  The
#' resulting data are the statistical input for cumulative/dynamic ROC and AUC
#' estimation.
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
#'   observed event or censoring outcome. Attribute `"risk_draws"` contains the
#'   aligned subject-by-draw event-risk matrix used for posterior ROC curves.
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
		context = "tvROC.JoiNMeFit()"
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
	keep_subject <- is.finite(risk_set$event_time) &
		risk_set$event_time > risk_set$time_start
	risk_draws <- .discrimination_horizon_risk_draws(
		prediction = pred,
		ids = ids,
		horizons = time_horizon_map
	)
	risk_set <- risk_set[keep_subject, , drop = FALSE]
	if (ncol(risk_draws) > 0L) {
		risk_draws <- risk_draws[keep_subject, , drop = FALSE]
	}
	rownames(risk_set) <- NULL
	attr(risk_set, "risk_draws") <- risk_draws
	risk_set
}

# ---- Time-dependent ROC and AUC ------------------------------------------

#' Time-dependent ROC curves and areas
#'
#' @description
#' Re-exports the `JMbayes2` generics [JMbayes2::tvROC()] and
#' [JMbayes2::tvAUC()]. JoiNMe supplies methods for fitted `JoiNMeFit` objects,
#' while ROC objects retain class `"tvROC"` and can therefore use the plotting,
#' printing, and AUC methods supplied by `JMbayes2`.
#'
#' `tvROC.JoiNMeFit()` estimates a cumulative/dynamic ROC curve using
#' longitudinal information available through a landmark time. Cases
#' experience the requested cause in the interval from the landmark to the
#' horizon; controls remain event-free beyond the horizon. Censoring is handled
#' by model-based expected status or Kaplan--Meier inverse-probability weights.
#'
#' Dynamic risks are calculated by [predict.JoiNMeFit()]. They therefore include
#' the fitted baseline hazard, survival covariates, longitudinal trajectories,
#' marker weights, and every identity, functional, monotone-spline, or ordered
#' piecewise-linear association.
#'
#' @param object A fitted `JoiNMeFit` object.
#' @param newdataLong Longitudinal evaluation data; defaults to the fitted data.
#' @param newdataEvent Event-process evaluation data; defaults to the fitted
#'   data.
#' @param time_start Numeric landmark time or vector of landmark times. The
#'   default is zero.
#' @param time_horizon Numeric horizon. It may be scalar or aligned with
#'   `time_start`. When this and `Dt` are omitted, the fitted maximum follow-up
#'   time is used.
#' @param Dt Positive prediction-window width used when `time_horizon` is
#'   omitted.
#' @param cause Positive integer identifying the event cause of interest.
#'   Other causes are treated as censoring for cause-specific discrimination.
#' @param n_samples Positive integer number of posterior draws used for dynamic
#'   prediction.
#' @param seed Integer seed used when posterior draws are subsampled.
#' @param type_weights Censoring treatment: `"model-based"` imputes the
#'   expected case/control status of subjects censored before the horizon;
#'   `"IPCW"` applies inverse Kaplan--Meier censoring weights to observable
#'   cases and controls.
#' @param ... Additional arguments passed to [predict.JoiNMeFit()].
#'
#' @details
#' For a threshold \eqn{c}, a subject is classified as a case when their
#' predicted conditional survival probability is below \eqn{c}. Sensitivity
#' and one minus specificity are evaluated at 101 thresholds from zero to one,
#' following `JMbayes2::tvROC()`. Posterior-draw curves are retained in `tp` and
#' `fp`; `TP` and `FP` are calculated from posterior mean risks.
#'
#' With model-based weighting, an observed case has case weight one, a subject
#' known to be event-free beyond the horizon has case weight zero, and a
#' subject censored within the prediction window contributes their model-based
#' conditional event probability. IPCW uses the Kaplan--Meier estimator of the
#' censoring survival distribution conditional on remaining under observation
#' at the landmark.
#'
#' A scalar landmark returns a standard `"tvROC"` object compatible with
#' `JMbayes2::tvAUC()` and `plot()`. Multiple landmarks return a
#' `"tvROC_JoiNMeFit_list"` containing one compatible curve per
#' landmark--horizon pair.
#'
#' @return `tvROC.JoiNMeFit()` returns a `"tvROC"` object, or a
#'   `"tvROC_JoiNMeFit_list"` for multiple landmarks. `tvAUC.JoiNMeFit()`
#'   returns the standard `"tvAUC"` object for one landmark and a data frame of
#'   areas for multiple landmarks.
#' @references
#' Heagerty PJ, Zheng Y (2005). Survival model predictive accuracy and ROC
#' curves. *Biometrics*, 61, 92--105.
#'
#' Rizopoulos D (2011). Dynamic predictions and prospective accuracy in joint
#' models. *Biometrics*, 67, 819--829.
#' @seealso [JMbayes2::tvROC()], [JMbayes2::tvAUC()],
#'   [predict.JoiNMeFit()], [concordance.JoiNMeFit()]
#' @name tvROC
#' @importFrom JMbayes2 tvROC tvAUC
#' @export
JMbayes2::tvROC

#' @rdname tvROC
#' @export
JMbayes2::tvAUC

#' Landmark and horizon vectors for time-dependent ROC analysis
#'
#' @param object A fitted `JoiNMeFit` object.
#' @param data_event Event-process evaluation data.
#' @param time_start Finite numeric landmark vector.
#' @param time_horizon Optional finite numeric horizon vector.
#' @param Dt Optional positive prediction-window width.
#'
#' @return A list containing aligned numeric `time_start` and `time_horizon`
#'   vectors.
#' @keywords internal
#' @noRd
.get_tvroc_times <- function(
	object,
	data_event,
	time_start,
	time_horizon,
	Dt
) {
	if (!is.numeric(time_start) || !length(time_start) || any(!is.finite(time_start))) {
		cli::cli_abort(
			"{.arg time_start} must be a finite numeric value or vector."
		)
	}
	time_start <- as.numeric(time_start)
	if (is.null(time_horizon)) {
		if (is.null(Dt)) {
			time_horizon <- rep(
				.default_discrimination_horizon(object, data_event),
				length(time_start)
			)
		} else {
			if (!is.numeric(Dt) || length(Dt) != 1L || !is.finite(Dt) || Dt <= 0) {
				cli::cli_abort("{.arg Dt} must be a positive finite scalar.")
			}
			time_horizon <- time_start + as.numeric(Dt)
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
		time_horizon <- rep(as.numeric(time_horizon), length.out = length(time_start))
	}
	if (any(time_horizon <= time_start)) {
		cli::cli_abort(
			"Every {.arg time_horizon} must be greater than its landmark time."
		)
	}
	list(time_start = time_start, time_horizon = time_horizon)
}

#' Kaplan-Meier curve
#'
#' @param fit A `survfit` object for the censoring distribution.
#' @param times Finite numeric evaluation times.
#' @param before Logical; evaluate immediately before each requested time.
#'
#' @return Numeric censoring-survival probabilities aligned with `times`.
#' @keywords internal
#' @noRd
.censoring_survival_at <- function(fit, times, before = FALSE) {
	times <- as.numeric(times)
	if (isTRUE(before)) {
		times <- times - sqrt(.Machine$double.eps) * pmax(1, abs(times))
	}
	index <- findInterval(times, as.numeric(fit$time))
	curve <- c(1, as.numeric(fit$surv))
	probability <- curve[index + 1L]
	probability[!is.finite(probability) | probability <= 0] <- NA_real_
	probability
}

#' Construct case and control weights for a dynamic ROC curve
#'
#' @param risk_set Subject-level dynamic risk set.
#' @param type_weights Character censoring method.
#' @param cause Positive integer identifying the event cause of interest.
#'
#' @return A list containing non-negative `case` and `control` weights.
#' @keywords internal
#' @noRd
.tvroc_status_weights <- function(risk_set, type_weights, cause) {
	type_weights <- match.arg(type_weights, c("model-based", "IPCW"))
	n_subjects <- nrow(risk_set)
	case_weight <- numeric(n_subjects)
	control_weight <- numeric(n_subjects)
	observed_case <- risk_set$event_window == 1L
	beyond_horizon <- risk_set$event_time > risk_set$time_horizon

	if (identical(type_weights, "model-based")) {
		unknown_status <- !observed_case & !beyond_horizon
		case_weight[observed_case] <- 1
		case_weight[unknown_status] <- pmin(
			pmax(as.numeric(risk_set$risk[unknown_status]), 0),
			1
		)
		control_weight <- 1 - case_weight
	} else {
		# For cause-specific ROC analysis, administrative censoring and competing
		# events terminate observation of the requested cause.  Requested-cause
		# events after the horizon would fail, rather than censoring events,
		# when estimating the censoring distribution.
		cause_failure <- risk_set$event_status == 1L &
			risk_set$event_type == cause
		censoring_event <- as.integer(!cause_failure)
		censoring_fit <- survival::survfit(
			survival::Surv(risk_set$event_time, censoring_event) ~ 1
		)
		if (any(observed_case)) {
			g_case <- .censoring_survival_at(
				censoring_fit,
				risk_set$event_time[observed_case],
				before = TRUE
			)
			case_weight[observed_case] <- ifelse(
				is.finite(g_case),
				1 / g_case,
				0
			)
		}
		if (any(beyond_horizon)) {
			g_control <- .censoring_survival_at(
				censoring_fit,
				risk_set$time_horizon[beyond_horizon]
			)
			control_weight[beyond_horizon] <- ifelse(
				is.finite(g_control),
				1 / g_control,
				0
			)
		}
	}
	list(case = case_weight, control = control_weight)
}

#' Select the central optimal ROC threshold
#'
#' @description
#' Identifies all thresholds attaining the largest finite criterion value and
#' returns their median. This reproduces the threshold convention used by
#' `JMbayes2::tvROC()` while remaining defined when the first threshold has an
#' indeterminate criterion, as occurs for the F1 score with no positive
#' classifications.
#'
#' @param thresholds Numeric vector of candidate survival-probability
#'   thresholds.
#' @param criterion Numeric criterion evaluated at the candidate thresholds.
#'
#' @return A numeric scalar, or `NA_real_` when no criterion value is finite.
#' @keywords internal
#' @noRd
.central_optimal_threshold <- function(thresholds, criterion) {
	finite <- is.finite(criterion)
	if (!any(finite)) {
		return(NA_real_)
	}
	best <- max(criterion[finite])
	stats::median(thresholds[finite & criterion == best])
}

#' Calculate one JMbayes2-compatible dynamic ROC curve
#'
#' @param risk_set Subject-level output from
#'   `.dynamic_discrimination_risk_set()`.
#' @param time_start Scalar landmark.
#' @param time_horizon Scalar prediction horizon.
#' @param type_weights Character censoring method.
#' @param cause Positive integer event cause.
#' @param object_name Character expression used to identify the fitted object.
#'
#' @return An object of class `"tvROC"` compatible with
#'   `JMbayes2::tvAUC()` and the JMbayes2 plotting method.
#' @keywords internal
#' @noRd
.tvroc_from_risk_set <- function(
	risk_set,
	time_start,
	time_horizon,
	type_weights,
	cause,
	object_name
) {
	required <- c(
		"risk", "event_window", "event_time", "event_status", "event_type",
		"time_horizon"
	)
	if (!all(required %in% names(risk_set)) || nrow(risk_set) == 0L) {
		cli::cli_abort("The dynamic ROC risk set contains no analysable subjects.")
	}

	# A missing posterior mean risk cannot define a threshold classification.
	# Remove such subjects once, before forming the case and control weights,
	# and retain the corresponding posterior-draw rows when they are available.
	risk_draws <- attr(risk_set, "risk_draws", exact = TRUE)
	complete_subject <- is.finite(risk_set$risk) &
		is.finite(risk_set$event_time) &
		is.finite(risk_set$time_horizon)
	if (!all(complete_subject)) {
		risk_set <- risk_set[complete_subject, , drop = FALSE]
		if (!is.null(risk_draws) && nrow(risk_draws) == length(complete_subject)) {
			risk_draws <- risk_draws[complete_subject, , drop = FALSE]
		}
	}
	if (nrow(risk_set) == 0L) {
		cli::cli_abort("The dynamic ROC risk set contains no finite predicted risks.")
	}
	if (!any(risk_set$event_window == 1L)) {
		cli::cli_abort(
			"No requested-cause event occurred between the landmark and horizon."
		)
	}

	thresholds <- seq(0, 1, length.out = 101L)
	weights <- .tvroc_status_weights(risk_set, type_weights, cause)
	case_total <- sum(weights$case)
	control_total <- sum(weights$control)
	if (case_total <= 0 || control_total <= 0) {
		cli::cli_abort(
			"Time-dependent ROC estimation requires positive case and control weight."
		)
	}

	# JMbayes2 formulates the threshold comparison on conditional survival:
	# smaller survival indicates greater event risk. Matrix multiplication
	# evaluates every threshold without constructing case-control pairs.
	mean_survival <- 1 - as.numeric(risk_set$risk)
	mean_positive <- outer(mean_survival, thresholds, "<")
	nTP <- as.numeric(crossprod(weights$case, mean_positive))
	nFP <- as.numeric(crossprod(weights$control, mean_positive))
	nFN <- case_total - nTP
	nTN <- control_total - nFP
	TP <- nTP / case_total
	FP <- nFP / control_total

	if (is.null(risk_draws) || ncol(risk_draws) == 0L) {
		risk_draws <- matrix(as.numeric(risk_set$risk), ncol = 1L)
	} else {
		# A draw containing a missing risk would make its complete ROC curve
		# undefined. Discard only those columns and retain the posterior mean
		# curve as a deterministic fallback if none remain.
		finite_draw <- colSums(!is.finite(risk_draws)) == 0L
		risk_draws <- risk_draws[, finite_draw, drop = FALSE]
		if (ncol(risk_draws) == 0L) {
			risk_draws <- matrix(as.numeric(risk_set$risk), ncol = 1L)
		}
	}
	tp <- matrix(NA_real_, nrow = length(thresholds), ncol = ncol(risk_draws))
	fp <- matrix(NA_real_, nrow = length(thresholds), ncol = ncol(risk_draws))
	for (draw_index in seq_len(ncol(risk_draws))) {
		draw_positive <- outer(1 - risk_draws[, draw_index], thresholds, "<")
		tp[, draw_index] <- as.numeric(crossprod(weights$case, draw_positive)) /
			case_total
		fp[, draw_index] <- as.numeric(crossprod(weights$control, draw_positive)) /
			control_total
	}

	f1 <- 2 * nTP / (2 * nTP + nFN + nFP)
	youden <- TP - FP
	f1_cutoff <- .central_optimal_threshold(thresholds, f1)
	youden_cutoff <- .central_optimal_threshold(thresholds, youden)
	out <- list(
		TP = TP,
		FP = FP,
		nTP = nTP,
		nFN = nFN,
		nFP = nFP,
		nTN = nTN,
		tp = tp,
		fp = fp,
		thrs = thresholds,
		thr = thresholds,
		F1score = f1_cutoff,
		Youden = youden_cutoff,
		Tstart = as.numeric(time_start),
		Thoriz = as.numeric(time_horizon),
		nr = nrow(risk_set),
		classObject = "JoiNMeFit",
		type_weights = type_weights,
		nameObject = object_name,
		cause = as.integer(cause)
	)
	class(out) <- "tvROC"
	out
}

#' @rdname tvROC
#' @export
tvROC.JoiNMeFit <- function(
	object,
	newdataLong = NULL,
	newdataEvent = NULL,
	time_start = 0,
	time_horizon = NULL,
	Dt = NULL,
	cause = 1,
	n_samples = 200,
	seed = 123,
	type_weights = c("model-based", "IPCW"),
	...
) {
	.require_fitted_survival_process(
		object,
		"Time-dependent ROC assessment"
	)
	.warn_experimental("tvROC")
	if (
		!is.numeric(cause) ||
			length(cause) != 1L ||
			!is.finite(cause) ||
			cause < 1 ||
			cause != as.integer(cause)
	) {
		cli::cli_abort("{.arg cause} must be a single positive integer.")
	}
	if (
		!is.numeric(n_samples) ||
			length(n_samples) != 1L ||
			!is.finite(n_samples) ||
			n_samples < 1 ||
			n_samples != as.integer(n_samples)
	) {
		cli::cli_abort("{.arg n_samples} must be a positive integer.")
	}
	type_weights <- match.arg(type_weights)
	data_in <- .get_train_data(
		object,
		newdataLong,
		newdataEvent,
		purpose = "time-dependent ROC"
	)
	times <- .get_tvroc_times(
		object = object,
		data_event = data_in$newdataEvent,
		time_start = time_start,
		time_horizon = time_horizon,
		Dt = Dt
	)

	curves <- vector("list", length(times$time_start))
	for (time_index in seq_along(curves)) {
		risk_set <- .dynamic_discrimination_risk_set(
			object = object,
			newdataLong = data_in$newdataLong,
			newdataEvent = data_in$newdataEvent,
			time_start = times$time_start[[time_index]],
			time_horizon = times$time_horizon[[time_index]],
			cause = as.integer(cause),
			n_samples = as.integer(n_samples),
			seed = seed,
			...
		)
		curves[[time_index]] <- .tvroc_from_risk_set(
			risk_set = risk_set,
			time_start = times$time_start[[time_index]],
			time_horizon = times$time_horizon[[time_index]],
			type_weights = type_weights,
			cause = as.integer(cause),
			object_name = deparse(substitute(object))
		)
	}
	if (length(curves) == 1L) {
		return(curves[[1L]])
	}
	names(curves) <- paste0(
		"t", format(times$time_start, trim = TRUE),
		"_h", format(times$time_horizon, trim = TRUE)
	)
	structure(
		list(
			curves = curves,
			Tstart = times$time_start,
			Thoriz = times$time_horizon,
			type_weights = type_weights,
			cause = as.integer(cause)
		),
		class = "tvROC_JoiNMeFit_list"
	)
}

#' @rdname tvROC
#' @export
tvAUC.JoiNMeFit <- function(
	object,
	newdataLong = NULL,
	newdataEvent = NULL,
	time_start = 0,
	time_horizon = NULL,
	Dt = NULL,
	cause = 1,
	n_samples = 200,
	seed = 123,
	type_weights = c("model-based", "IPCW"),
	...
) {
	roc <- tvROC(
		object = object,
		newdataLong = newdataLong,
		newdataEvent = newdataEvent,
		time_start = time_start,
		time_horizon = time_horizon,
		Dt = Dt,
		cause = cause,
		n_samples = n_samples,
		seed = seed,
		type_weights = type_weights,
		...
	)
	if (inherits(roc, "tvROC")) {
		return(JMbayes2::tvAUC(roc))
	}

	areas <- vector("list", length(roc$curves))
	time_start_out <- numeric(length(roc$curves))
	time_horizon_out <- numeric(length(roc$curves))
	auc_out <- numeric(length(roc$curves))
	n_subjects_out <- integer(length(roc$curves))
	type_weights_out <- character(length(roc$curves))
	for (curve_index in seq_along(roc$curves)) {
		areas[[curve_index]] <- JMbayes2::tvAUC(roc$curves[[curve_index]])
		time_start_out[[curve_index]] <- areas[[curve_index]]$Tstart
		time_horizon_out[[curve_index]] <- areas[[curve_index]]$Thoriz
		auc_out[[curve_index]] <- areas[[curve_index]]$auc
		n_subjects_out[[curve_index]] <- areas[[curve_index]]$nr
		type_weights_out[[curve_index]] <- areas[[curve_index]]$type_weights
	}
	out <- data.frame(
		time_start = time_start_out,
		time_horizon = time_horizon_out,
		auc = auc_out,
		n_subjects = n_subjects_out,
		type_weights = type_weights_out,
		stringsAsFactors = FALSE
	)
	class(out) <- c("tvAUC_JoiNMeFit", "data.frame")
	out
}

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
	
	draws_obj <- posterior_draws(object, format = "draws_array")
	vars <- posterior::variables(draws_obj)
	vars <- .filter_diag_vars(vars, pars, regex_pars)
	if (length(vars) == 0) {
		return(data.frame(variable = character(0), rhat = numeric(0)))
	}
	draws_obj <- posterior_draws(
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
	draws_obj <- posterior_draws(object, format = "draws_array")
	vars <- posterior::variables(draws_obj)
	vars <- .filter_diag_vars(vars, pars, regex_pars)
	if (length(vars) == 0) {
		return(data.frame(variable = character(0), ess = numeric(0)))
	}
	draws_obj <- posterior_draws(
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
	draws_obj <- posterior_draws(object, format = "draws_array")
	vars <- posterior::variables(draws_obj)
	vars <- .filter_diag_vars(vars, pars, regex_pars)
	if (length(vars) == 0) {
		return(data.frame(variable = character(0), mcse = numeric(0)))
	}
	draws_obj <- posterior_draws(
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
