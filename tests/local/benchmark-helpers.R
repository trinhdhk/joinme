benchmark_config <- function(
  seed = 20260717L,
  n_rep = 20L,
  n_id = 40L,
  n_markers = 40L,
  iter_warmup_joinme = 500L,
  iter_sampling_joinme = 700L,
  n_iter_jm = 18000L,
  n_burnin_jm = 2000L,
  n_chains_jm = 2L,
  n_thin_jm = 10L,
  output_dir = ".artifacts",
  use_mirai = TRUE,
  n_workers = 10L,
  tmax = 5,
  benchmark_stamp = format(Sys.time(), "%Y%m%d-%H%M%S")
) {
  list(
    seed = as.integer(seed),
    n_rep = as.integer(n_rep),
    n_id = as.integer(n_id),
    n_markers = as.integer(n_markers),
    iter_warmup_joinme = as.integer(iter_warmup_joinme),
    iter_sampling_joinme = as.integer(iter_sampling_joinme),
    n_iter_jm = as.integer(n_iter_jm),
    n_burnin_jm = as.integer(n_burnin_jm),
    n_chains_jm = as.integer(n_chains_jm),
    n_thin_jm = as.integer(n_thin_jm),
    output_dir = as.character(output_dir),
    use_mirai = isTRUE(use_mirai),
    n_workers = as.integer(n_workers),
    tmax = as.numeric(tmax),
    benchmark_stamp = as.character(benchmark_stamp)
  )
}

cfg <- benchmark_config(
  seed = 20260330L,
  n_rep = 500L,
  n_id = 50L,
  n_markers = 40L,
  iter_warmup_joinme = 500L,
  iter_sampling_joinme = 1200L,
  n_iter_jm = 25000L,
  n_burnin_jm = 2000L,
  n_chains_jm = 2L,
  n_thin_jm = 10L,
  output_dir = ".artifacts",
  use_mirai = TRUE,
  n_workers = 4L,
  benchmark_stamp = 'test-jul'
)

if (!dir.exists(cfg$output_dir)) {
  dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
}

# benchmark_stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
benchmark_stamp <- cfg$benchmark_stamp

base_name <- paste0(
  "joinme-vs-jmbayes2-",
  cfg$n_markers,
  "markers-",
  benchmark_stamp
)

output_rds <- file.path(
  cfg$output_dir,
  paste0(base_name, '.rds')
)

# Keep the longitudinal basis identical across simulation, joinme, and the
# per-marker JMbayes2 submodels so term-level recovery compares like with like.
# benchmark_long_basis <- paste0("splines::bs(time, knots = 1, Boundary.knots = c(0, ", cfg$tmax, "), degree = 2)")
benchmark_long_basis <-
  paste0(
    'splines::ns(time, knots = 1, Boundary.knots = c(0, ', cfg$tmax, '))'
  )
benchmark_formula_long <- stats::as.formula(
  paste0(
    "y ~ 1 + ",
    benchmark_long_basis,
    " + (1 | id) + (1 + ",
    benchmark_long_basis,
    " + (1 | id) | marker)"
  )
)
benchmark_formula_jm_long <- function(outcome_var) {
  stats::as.formula(paste0(outcome_var, " ~ 1 + ", benchmark_long_basis))
}
benchmark_beta_long <- stats::setNames(
  # c(1.0, 0.5, 0.6, 0.7),
  c(1.0, 0.5, 0.2),
  c("(Intercept)", paste0(benchmark_long_basis, 1:2))
)
benchmark_marker_re_params <- list(
  sd = c(1.0, 0.25, 0.2), #0.15),
  corr = diag(3)
)

safe_cor <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 2L) {
    return(NA_real_)
  }
  stats::cor(x[keep], y[keep])
}

null_or <- function(x, default) {
  if (is.null(x)) default else x
}

time_try <- function(expr) {
  value <- NULL
  err <- NULL
  elapsed <- system.time({
    value <- tryCatch(
      eval.parent(substitute(expr)),
      error = function(e) {
        err <<- conditionMessage(e)
        NULL
      }
    )
  })[["elapsed"]]
  list(value = value, error = err, elapsed = unname(elapsed))
}

normalize_marker_draws <- function(marker_draws, marker_levels, marker_terms) {
  if (is.data.frame(marker_draws)) {
    marker_draws <- as.matrix(marker_draws)
  }
  if (is.null(dim(marker_draws))) {
    marker_draws <- matrix(marker_draws, ncol = length(marker_terms))
  }
  if (!is.matrix(marker_draws)) {
    stop("sim$truth$re_draws$marker must be a matrix-like object.")
  }
  if (
    nrow(marker_draws) != length(marker_levels) &&
      ncol(marker_draws) == length(marker_levels)
  ) {
    marker_draws <- t(marker_draws)
  }
  if (nrow(marker_draws) != length(marker_levels)) {
    stop(
      "Could not align sim$truth$re_draws$marker with marker levels: ",
      nrow(marker_draws),
      " rows for ",
      length(marker_levels),
      " markers."
    )
  }
  if (ncol(marker_draws) != length(marker_terms)) {
    stop(
      "Could not align sim$truth$re_draws$marker with marker terms: ",
      ncol(marker_draws),
      " columns for ",
      length(marker_terms),
      " terms."
    )
  }
  rownames(marker_draws) <- marker_levels
  colnames(marker_draws) <- marker_terms
  marker_draws
}

infer_marker_terms <- function(sim, marker_levels) {
  marker_draws <- sim$truth$re_draws$marker
  if (is.data.frame(marker_draws)) {
    marker_draws <- as.matrix(marker_draws)
  }
  if (is.null(dim(marker_draws))) {
    term_count <- 1L
    term_names <- names(marker_draws)
  } else if (nrow(marker_draws) == length(marker_levels)) {
    term_count <- ncol(marker_draws)
    term_names <- colnames(marker_draws)
  } else if (ncol(marker_draws) == length(marker_levels)) {
    term_count <- nrow(marker_draws)
    term_names <- rownames(marker_draws)
  } else {
    term_count <- ncol(marker_draws)
    term_names <- colnames(marker_draws)
  }

  if (is.null(term_names) || any(!nzchar(term_names))) {
    beta_names <- names(sim$truth$beta_long)
    if (length(beta_names) < term_count) {
      stop(
        "Could not infer marker-level longitudinal terms from simulation truth."
      )
    }
    term_names <- beta_names[seq_len(term_count)]
  }
  as.character(term_names)
}

build_truth_long <- function(sim, marker_levels) {
  marker_terms <- infer_marker_terms(sim, marker_levels)
  marker_draws <- normalize_marker_draws(
    sim$truth$re_draws$marker,
    marker_levels,
    marker_terms
  )
  fixed_truth <- as.numeric(sim$truth$beta_long[marker_terms])
  fixed_truth[is.na(fixed_truth)] <- 0
  do.call(
    rbind,
    lapply(marker_levels, function(marker_label) {
      data.frame(
        marker = marker_label,
        term = marker_terms,
        truth = fixed_truth +
          as.numeric(marker_draws[marker_label, marker_terms, drop = TRUE]),
        stringsAsFactors = FALSE
      )
    })
  )
}

infer_survival_terms <- function(sim) {
  beta_event <- sim$truth$beta_event
  term_names <- names(beta_event)
  if (
    !is.null(term_names) && length(term_names) > 0L && all(nzchar(term_names))
  ) {
    return(as.character(term_names))
  }

  formula_event <- sim$truth$formulaEvent
  event_cols <- colnames(stats::model.matrix(formula_event, sim$dataEvent))
  as.character(event_cols %||% character(0))
}

build_truth_survival <- function(sim) {
  survival_terms <- infer_survival_terms(sim)
  truth_vals <- rep(0, length(survival_terms))
  names(truth_vals) <- survival_terms

  beta_event <- sim$truth$beta_event
  if (!is.null(beta_event) && length(beta_event) > 0L) {
    beta_names <- names(beta_event)
    if (!is.null(beta_names) && any(nzchar(beta_names))) {
      keep <- intersect(survival_terms, beta_names)
      if (length(keep) > 0L) {
        truth_vals[keep] <- as.numeric(beta_event[keep])
      }
    } else {
      n_copy <- min(length(survival_terms), length(beta_event))
      if (n_copy > 0L) {
        truth_vals[seq_len(n_copy)] <- as.numeric(beta_event)[seq_len(n_copy)]
      }
    }
  }

  data.frame(
    term = survival_terms,
    truth = as.numeric(truth_vals),
    stringsAsFactors = FALSE
  )
}

build_truth_assoc <- function(sim, marker_levels) {
  assoc_coefs <- sim$truth$assoc_coefs
  marker_weights_by_term <- null_or(
    sim$truth$marker_weights_by_term,
    list()
  )
  default_weights <- as.numeric(null_or(
    sim$truth$marker_weights,
    rep(1, length(marker_levels))
  ))
  cv_marker_weights <- as.numeric(null_or(
    marker_weights_by_term$cv_marker,
    default_weights
  ))
  if (length(cv_marker_weights) != length(marker_levels)) {
    stop("cv_marker weights do not align with marker levels.")
  }
  cv_mean_coef <- as.numeric(null_or(assoc_coefs[["cv_mean"]], 0))
  cv_marker_coef <- as.numeric(null_or(assoc_coefs[["cv_marker"]], 0))
  n_markers <- length(marker_levels) # normalising divisor used by the fitted weighted marker aggregate
  data.frame(
    marker = marker_levels,
    truth = (cv_mean_coef + cv_marker_coef * cv_marker_weights) / n_markers,
    stringsAsFactors = FALSE
  )
}

summarise_draw_matrix <- function(draw_matrix, id_col) {
  if (
    is.null(draw_matrix) || !is.matrix(draw_matrix) || ncol(draw_matrix) == 0L
  ) {
    return(data.frame())
  }
  out <- data.frame(
    id = colnames(draw_matrix),
    estimate = colMeans(draw_matrix),
    std_error = apply(draw_matrix, 2, stats::sd),
    conf_low = apply(draw_matrix, 2, stats::quantile, probs = 0.025),
    conf_high = apply(draw_matrix, 2, stats::quantile, probs = 0.975),
    stringsAsFactors = FALSE
  )
  names(out)[names(out) == "id"] <- id_col
  out
}

summarise_draw_array_with_rhat <- function(draw_array, id_col, labels) {
  if (
    is.null(draw_array) ||
      length(dim(draw_array)) != 3L ||
      dim(draw_array)[3] == 0L
  ) {
    return(data.frame())
  }
  dimnames(draw_array)[[3L]] <- labels
  summary_tbl <- posterior::summarise_draws(
    posterior::as_draws_array(draw_array),
    "mean",
    "sd",
    ~posterior::quantile2(.x, probs = c(0.025, 0.975)),
    "rhat"
  )
  out <- data.frame(
    id = labels,
    estimate = summary_tbl$mean,
    std_error = summary_tbl$sd,
    conf_low = summary_tbl$q2.5,
    conf_high = summary_tbl$q97.5,
    rhat = summary_tbl$rhat,
    stringsAsFactors = FALSE
  )
  names(out)[names(out) == "id"] <- id_col
  out
}

bind_rows_or_empty <- function(x) {
  x <- Filter(function(obj) !is.null(obj) && nrow(obj) > 0L, x)
  if (!length(x)) {
    return(data.frame())
  }
  do.call(rbind, x)
}

draw_array_slice <- function(draw_array, index = 1L) {
  out <- draw_array[,, index, drop = FALSE]
  dim(out) <- dim(draw_array)[1:2]
  out
}

build_joinme_assoc_coef_draws <- function(
  fit_joinme,
  marker_levels,
  draws = NULL,
  seed = 1
) {
  association_raw <- tryCatch(
    extract(
      fit_joinme,
      what = "raw",
      variable = c("alpha_cv_mean", "alpha_cv_marker"),
      draws = draws,
      seed = seed,
      keep_chains = TRUE
    )$posterior_draws,
    error = function(error) NULL
  )
  available_variables <- if (is.null(association_raw)) character(0) else dimnames(association_raw)[[3L]]
  cv_mean_var <- if ("alpha_cv_mean" %in% available_variables) "alpha_cv_mean" else NULL
  cv_marker_var <- if ("alpha_cv_marker" %in% available_variables) "alpha_cv_marker" else NULL

  cv_mean_arr <- if (!is.null(cv_mean_var)) {
    association_raw[, , cv_mean_var, drop = FALSE]
  } else {
    NULL
  }
  cv_marker_arr <- if (!is.null(cv_marker_var)) {
    association_raw[, , cv_marker_var, drop = FALSE]
  } else {
    NULL
  }
  weight_arr <- if (!is.null(cv_marker_var)) {
    extract(
      fit_joinme,
      what = "marker_weights",
      term = "cv_marker",
      draws = draws,
      seed = seed,
      keep_chains = TRUE
    )$posterior_draws
  } else {
    NULL
  }

  if (is.null(cv_mean_arr) && is.null(cv_marker_arr)) {
    return(NULL)
  }

  template_arr <- null_or(cv_mean_arr, cv_marker_arr)
  n_markers <- length(marker_levels)
  out_arr <- array(
    0,
    dim = c(dim(template_arr)[1], dim(template_arr)[2], n_markers),
    dimnames = list(
      iteration = dimnames(template_arr)[[1]],
      chain = dimnames(template_arr)[[2]],
      variable = marker_levels
    )
  )

  if (!is.null(cv_mean_arr)) {
    cv_mean_mat <- draw_array_slice(cv_mean_arr, 1L)
    for (marker_index in seq_len(n_markers)) {
      out_arr[,, marker_index] <- out_arr[,, marker_index] + cv_mean_mat
    }
  }
  if (!is.null(cv_marker_arr)) {
    if (is.null(weight_arr)) {
      stop(
        "Could not recover marker weights for cv_marker association coefficient reconstruction."
      )
    }
    cv_marker_mat <- draw_array_slice(cv_marker_arr, 1L)
    for (marker_index in seq_len(n_markers)) {
      weight_mat <- draw_array_slice(weight_arr, marker_index)
      out_arr[,, marker_index] <- out_arr[,, marker_index] +
        cv_marker_mat * weight_mat
    }
  }

  out_arr <- out_arr / n_markers # convert the mean and marker channels to coefficients of the separate marker trajectories used by JMbayes2

  posterior::as_draws_matrix(posterior::as_draws_array(out_arr))
}

build_joinme_assoc_coef_array <- function(
  fit_joinme,
  marker_levels,
  draws = NULL,
  seed = 1
) {
  association_raw <- tryCatch(
    extract(
      fit_joinme,
      what = "raw",
      variable = c("alpha_cv_mean", "alpha_cv_marker"),
      draws = draws,
      seed = seed,
      keep_chains = TRUE
    )$posterior_draws,
    error = function(error) NULL
  )
  available_variables <- if (is.null(association_raw)) character(0) else dimnames(association_raw)[[3L]]
  cv_mean_var <- if ("alpha_cv_mean" %in% available_variables) "alpha_cv_mean" else NULL
  cv_marker_var <- if ("alpha_cv_marker" %in% available_variables) "alpha_cv_marker" else NULL

  cv_mean_arr <- if (!is.null(cv_mean_var)) {
    association_raw[, , cv_mean_var, drop = FALSE]
  } else {
    NULL
  }
  cv_marker_arr <- if (!is.null(cv_marker_var)) {
    association_raw[, , cv_marker_var, drop = FALSE]
  } else {
    NULL
  }
  weight_arr <- if (!is.null(cv_marker_var)) {
    extract(
      fit_joinme,
      what = "marker_weights",
      term = "cv_marker",
      draws = draws,
      seed = seed,
      keep_chains = TRUE
    )$posterior_draws
  } else {
    NULL
  }

  if (is.null(cv_mean_arr) && is.null(cv_marker_arr)) {
    return(NULL)
  }

  template_arr <- null_or(cv_mean_arr, cv_marker_arr)
  n_markers <- length(marker_levels)
  out_arr <- array(
    0,
    dim = c(dim(template_arr)[1], dim(template_arr)[2], n_markers),
    dimnames = list(
      iteration = dimnames(template_arr)[[1]],
      chain = dimnames(template_arr)[[2]],
      variable = marker_levels
    )
  )

  if (!is.null(cv_mean_arr)) {
    cv_mean_mat <- draw_array_slice(cv_mean_arr, 1L)
    for (marker_index in seq_len(n_markers)) {
      out_arr[,, marker_index] <- out_arr[,, marker_index] + cv_mean_mat
    }
  }
  if (!is.null(cv_marker_arr)) {
    if (is.null(weight_arr)) {
      stop(
        "Could not recover marker weights for cv_marker association coefficient reconstruction."
      )
    }
    cv_marker_mat <- draw_array_slice(cv_marker_arr, 1L)
    for (marker_index in seq_len(n_markers)) {
      weight_mat <- draw_array_slice(weight_arr, marker_index)
      out_arr[,, marker_index] <- out_arr[,, marker_index] +
        cv_marker_mat * weight_mat
    }
  }

  out_arr <- out_arr / n_markers # match the marker-count normalisation in the fitted hazard exactly

  out_arr
}

#' Check the marker-weight data-generating parameterisation used by a benchmark
#'
#' The comparison is meaningful only when simulation and fitting use the same
#' latent marker coordinates. This check verifies that a Normal family draw is
#' inserted directly into the effective weight, with no unrecorded multiplier,
#' and that the reported effective weights equal offset plus common location
#' plus that departure. It also guards the marker ordering used by both fitted
#' models before a computationally expensive replication begins.
#'
#' @param simulation Result returned by `simulate_joinme()`.
#' @param marker_levels Marker labels in the benchmark fitting order.
#'
#' @return The unchanged simulation result, suitable as the final expression
#'   of `simulate_replication()`.
validate_benchmark_marker_weight_dgp <- function(simulation, marker_levels) {
  truth <- simulation$truth # complete data-generating record returned by simulate_joinme()
  direct_departures <- truth$z_marker_weight_sets # set-by-marker departures entering effective weights directly
  stored_departures <- truth$marker_weight_standardised_sets # same family-transformed coordinates retained for recovery checks
  if (!is.matrix(direct_departures) || ncol(direct_departures) != length(marker_levels)) {
    stop("Simulated marker-weight departures do not align with the benchmark markers.")
  }
  if (!isTRUE(all.equal(direct_departures, stored_departures, tolerance = 0))) {
    stop("The simulator has transformed marker-weight departures after drawing them.")
  }
  if (!isTRUE(all.equal(as.numeric(t(direct_departures)), truth$marker_weight_prior_raw, tolerance = 0))) {
    stop("The Normal marker-weight draw no longer matches the coordinate used by the fitted model.")
  }
  if (!any(abs(direct_departures) > sqrt(.Machine$double.eps))) {
    stop("The benchmark requires randomly drawn marker-weight departures, but all departures are zero.")
  }
  effective_weights <- as.numeric(truth$marker_weights_by_term$cv_marker) # marker weights that generated the event process
  declared_offsets <- as.numeric(truth$marker_weights_offset_by_term$cv_marker) # fixed marker-specific reference values
  common_location <- as.numeric(truth$marker_weight_mean_by_term$cv_marker) # population location drawn once because marker_weight_mean is NULL
  if (length(common_location) != 1L || !is.finite(common_location)) {
    stop("The benchmark requires one finite randomly generated marker-weight mean.")
  }
  expected_weights <- declared_offsets + common_location + as.numeric(direct_departures[1L, ]) # direct offset-plus-location-plus-departure parameterisation
  if (!isTRUE(all.equal(effective_weights, expected_weights, tolerance = 1e-12))) {
    stop("The benchmark marker weights do not follow offset + common location + direct departure.")
  }
  if ("marker_weight_scale" %in% names(truth) || "marker_weight_scale" %in% names(truth$stan_fit)) {
    stop("The simulated truth unexpectedly contains a separate marker-weight scale.")
  }
  simulation
}

#' Summarise the finite-marker realisation of the marker-weight hierarchy
#'
#' The population location, average realised departure and average effective
#' weight are different statistical quantities when only finitely many markers
#' are sampled. This compact table prevents a recovery benchmark from treating
#' the realised average weight as if it had to equal its generating population
#' location. Increasing the number of subjects improves information about the
#' realised weights, but does not increase the number of marker-level draws
#' informing their population location.
#'
#' @param simulation Result returned by `simulate_joinme()`.
#' @param term Weighted association term represented in the benchmark.
#'
#' @return One-row data frame containing the population location and the two
#'   finite-marker averages.
benchmark_marker_weight_truth <- function(simulation, term = "cv_marker") {
  truth <- simulation$truth # complete data-generating record for this replication
  population_location <- as.numeric(truth$marker_weight_mean_by_term[[term]]) # location held fixed while markers and subjects are generated
  realised_departures <- as.numeric(truth$marker_weights_latent_by_term[[term]]) # finite collection of marker-specific random departures
  effective_weights <- as.numeric(truth$marker_weights_by_term[[term]]) # offset-plus-location-plus-departure values entering the event process
  data.frame(
    assoc_term = term,
    n_markers = length(effective_weights),
    population_mean_weight = population_location,
    realised_departure_mean = mean(realised_departures),
    realised_effective_weight_mean = mean(effective_weights),
    stringsAsFactors = FALSE
  )
}

parse_assoc_marker <- function(x, fallback = NULL) {
  if (is.null(x)) {
    return(fallback)
  }
  out <- as.character(x)
  out <- gsub("^association\\.?", "", out)
  out <- gsub("^value\\(", "", out)
  out <- gsub("\\)$", "", out)
  out <- gsub("^y_", "", out)
  out <- gsub(".*y_", "", out)
  out[nchar(out) == 0L] <- fallback[nchar(out) == 0L]
  out
}

combine_truth_and_estimates <- function(
  truth_df,
  est_df,
  by_cols,
  replication,
  engine
) {
  out <- merge(truth_df, est_df, by = by_cols, all.x = TRUE, sort = FALSE)
  out$replication <- replication
  out$engine <- engine
  out$error <- out$estimate - out$truth
  out
}

extract_joinme_long <- function(fit_joinme, truth_long, replication) {
  coef_summary <- coef(fit_joinme, summary = TRUE)
  est_df <- coef_summary$formulaLong$marker[, c(
    "marker",
    "term",
    "Estimate",
    "Est.Error",
    "Q2.5",
    "Q97.5",
    "Rhat"
  )]
  names(est_df) <- c(
    "marker",
    "term",
    "estimate",
    "std_error",
    "conf_low",
    "conf_high",
    "rhat"
  )
  combine_truth_and_estimates(
    truth_long,
    est_df,
    c("marker", "term"),
    replication,
    "joinme"
  )
}

extract_joinme_survival <- function(fit_joinme, truth_survival, replication) {
  coef_summary <- coef(fit_joinme, summary = TRUE)
  est_df <- coef_summary$formulaEvent[, c(
    "term",
    "Estimate",
    "Est.Error",
    "Q2.5",
    "Q97.5",
    "Rhat"
  )]
  names(est_df) <- c(
    "term",
    "estimate",
    "std_error",
    "conf_low",
    "conf_high",
    "rhat"
  )
  combine_truth_and_estimates(
    truth_survival,
    est_df,
    "term",
    replication,
    "joinme"
  )
}

extract_joinme_assoc <- function(
  fit_joinme,
  truth_assoc,
  marker_levels,
  replication,
  seed
) {
  draw_arr <- build_joinme_assoc_coef_array(
    fit_joinme,
    marker_levels,
    seed = seed
  )
  est_df <- summarise_draw_array_with_rhat(draw_arr, "marker", marker_levels)
  combine_truth_and_estimates(
    truth_assoc,
    est_df,
    "marker",
    replication,
    "joinme"
  )
}

extract_jm_long <- function(
  jm_summary,
  truth_long,
  marker_levels,
  replication
) {
  truth_terms <- unique(as.character(truth_long$term))
  outcome_tabs <- grep("^Outcome", names(jm_summary), value = TRUE)
  est_df <- do.call(
    rbind,
    lapply(seq_along(outcome_tabs), function(idx) {
      tab <- jm_summary[[outcome_tabs[idx]]]
      if (is.null(tab) || !nrow(tab)) {
        return(NULL)
      }
      out <- data.frame(
        marker = marker_levels[idx],
        term = rownames(tab),
        estimate = as.numeric(tab[, "Mean"]),
        std_error = as.numeric(tab[, "StDev"]),
        conf_low = as.numeric(tab[, "2.5%"]),
        conf_high = as.numeric(tab[, "97.5%"]),
        rhat = if ("Rhat" %in% names(tab)) {
          as.numeric(tab[, "Rhat"])
        } else {
          NA_real_
        },
        stringsAsFactors = FALSE
      )
      out[out$term %in% truth_terms, , drop = FALSE]
    })
  )
  combine_truth_and_estimates(
    truth_long,
    est_df,
    c("marker", "term"),
    replication,
    "JMbayes2"
  )
}

extract_jm_survival <- function(jm_summary, truth_survival, replication) {
  tab <- jm_summary$Survival
  if (is.null(tab) || !nrow(tab)) {
    return(data.frame())
  }
  truth_terms <- unique(as.character(truth_survival$term))
  est_df <- data.frame(
    term = rownames(tab),
    estimate = as.numeric(tab[, "Mean"]),
    std_error = as.numeric(tab[, "StDev"]),
    conf_low = as.numeric(tab[, "2.5%"]),
    conf_high = as.numeric(tab[, "97.5%"]),
    rhat = if ("Rhat" %in% names(tab)) as.numeric(tab[, "Rhat"]) else NA_real_,
    stringsAsFactors = FALSE
  )
  est_df <- est_df[est_df$term %in% truth_terms, , drop = FALSE]
  combine_truth_and_estimates(
    truth_survival,
    est_df,
    "term",
    replication,
    "JMbayes2"
  )
}

extract_jm_assoc <- function(
  jm_summary,
  truth_assoc,
  truth_survival,
  marker_levels,
  replication
) {
  tab <- jm_summary$Survival
  if (is.null(tab) || !nrow(tab)) {
    return(data.frame())
  }
  surv_terms <- unique(as.character(truth_survival$term))
  keep <- !(rownames(tab) %in% surv_terms)
  if (!any(keep)) {
    return(data.frame())
  }
  est_df <- data.frame(
    marker = parse_assoc_marker(
      rownames(tab)[keep],
      fallback = marker_levels[seq_len(sum(keep))]
    ),
    estimate = as.numeric(tab[keep, "Mean"]),
    std_error = as.numeric(tab[keep, "StDev"]),
    conf_low = as.numeric(tab[keep, "2.5%"]),
    conf_high = as.numeric(tab[keep, "97.5%"]),
    rhat = if ("Rhat" %in% names(tab)) {
      as.numeric(tab[keep, "Rhat"])
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
  combine_truth_and_estimates(
    truth_assoc,
    est_df,
    "marker",
    replication,
    "JMbayes2"
  )
}

summarise_runtime_metrics <- function(runtime_records, total_reps) {
  if (!nrow(runtime_records)) {
    return(data.frame())
  }
  split_idx <- split(seq_len(nrow(runtime_records)), runtime_records$engine)
  out <- lapply(split_idx, function(idx) {
    block <- runtime_records[idx, , drop = FALSE]
    data.frame(
      engine = block$engine[1],
      n_rep = total_reps,
      mean_runtime = mean(block$elapsed_seconds),
      sd_runtime = stats::sd(block$elapsed_seconds),
      median_runtime = stats::median(block$elapsed_seconds),
      success_rate = 100 * mean(block$success),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

summarise_replication_metrics <- function(records, group_cols, total_reps) {
  if (!nrow(records)) {
    return(data.frame())
  }
  split_idx <- split(
    seq_len(nrow(records)),
    interaction(records[c("engine", group_cols)], drop = TRUE, lex.order = TRUE)
  )
  out <- lapply(split_idx, function(idx) {
    block <- records[idx, , drop = FALSE]
    available <- is.finite(block$estimate)
    covered <- available &
      is.finite(block$conf_low) &
      is.finite(block$conf_high) &
      block$conf_low <= block$truth &
      block$truth <= block$conf_high

    est_vals <- block$estimate[available]
    err_vals <- block$error[available]
    se_vals <- block$std_error[available]
    rhat_vals <- block$rhat[available & is.finite(block$rhat)]
    bad_rhat <- available & is.finite(block$rhat) & block$rhat > 1.01

    data.frame(
      block[1, c("engine", group_cols), drop = FALSE],
      n_rep = total_reps,
      n_available = sum(available),
      success_rate = 100 * sum(available) / total_reps,
      truth_mean = mean(block$truth),
      estimate_mean = if (length(est_vals)) mean(est_vals) else NA_real_,
      bias = if (length(err_vals)) mean(err_vals) else NA_real_,
      average_std_error = if (length(se_vals)) mean(se_vals) else NA_real_,
      asymptotic_std_error = if (length(est_vals) > 1L) {
        stats::sd(est_vals)
      } else {
        NA_real_
      },
      average_rhat = if (length(rhat_vals)) mean(rhat_vals) else NA_real_,
      max_rhat = if (length(rhat_vals)) max(rhat_vals) else NA_real_,
      bad_rhat_rate = 100 * sum(bad_rhat) / total_reps,
      coverage_rate = 100 * sum(covered) / total_reps,
      coverage_rate_available = if (sum(available) > 0L) {
        100 * mean(covered[available])
      } else {
        NA_real_
      },
      cor = if (length(est_vals) > 1L) {
        safe_cor(block$truth[available], est_vals)
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

build_rhat_comparison <- function(summary_df, key_cols) {
  if (!nrow(summary_df)) {
    return(data.frame())
  }
  keep_cols <- c(
    key_cols,
    "engine",
    "average_rhat",
    "max_rhat",
    "bad_rhat_rate"
  )
  wide <- summary_df[, keep_cols, drop = FALSE]
  split_df <- split(wide, wide$engine)
  joinme_df <- split_df$joinme
  jm_df <- split_df$JMbayes2
  if (is.null(joinme_df) || is.null(jm_df)) {
    return(data.frame())
  }
  names(joinme_df)[
    names(joinme_df) %in% c("average_rhat", "max_rhat", "bad_rhat_rate")
  ] <- paste0(c("average_rhat", "max_rhat", "bad_rhat_rate"), "_joinme")
  names(jm_df)[
    names(jm_df) %in% c("average_rhat", "max_rhat", "bad_rhat_rate")
  ] <- paste0(c("average_rhat", "max_rhat", "bad_rhat_rate"), "_JMbayes2")
  joinme_df$engine <- NULL
  jm_df$engine <- NULL
  out <- merge(joinme_df, jm_df, by = key_cols, all = TRUE, sort = FALSE)
  if (all(c("average_rhat_joinme", "average_rhat_JMbayes2") %in% names(out))) {
    out$average_rhat_diff <- out$average_rhat_joinme - out$average_rhat_JMbayes2
  }
  if (all(c("max_rhat_joinme", "max_rhat_JMbayes2") %in% names(out))) {
    out$max_rhat_diff <- out$max_rhat_joinme - out$max_rhat_JMbayes2
  }
  out
}

simulate_replication <- function(seed) {
  simulation <- simulate_joinme(
    formulaLong = benchmark_formula_long,
    formulaEvent = survival::Surv(time, event) ~ x1 + x2,
    n_id = cfg$n_id,
    families = rep("gaussian", cfg$n_markers),
    times_obs = seq(0, cfg$tmax, length.out = 5),
    assoc = c("cv_mean", "cv_marker"),
    truth = jm_truth(
      longitudinal = list(
        intercept = benchmark_beta_long["(Intercept)"],
        slope = benchmark_beta_long[setdiff(names(benchmark_beta_long), "(Intercept)")]
      ),
      survival = list(slope = c(x1 = -0.5, x2 = 0.25)),
      assoc_coef = list(slope = c(cv_mean = 1, cv_marker = 1)),
      vcov = list(
        sd = list(intercept = -0.5, latent = 1),
        corr = list(intercept = numeric(0), latent = numeric(0))
      ),
      marker_weights = list(
        # Centred unit-scale family for the random departures used directly
        # around the offsets; the association slope supplies hazard magnitude.
        intercept = 0.5,
        family = "normal"
      ),
      basehaz = list(type = "weibull", shape = 0.8, scale = 12),
      re_params = list(
        id = list(sd = c(1)),
        marker = benchmark_marker_re_params
      )
    ),
    use_mirai = cfg$use_mirai,
    n_workers = cfg$n_workers,
    seed = seed
  )
  validate_benchmark_marker_weight_dgp(
    simulation,
    paste0("m", seq_len(cfg$n_markers))
  ) # default marker labels generated by simulate_joinme()
}
