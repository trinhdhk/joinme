options(error = function() {
  traceback(2)
  quit(status = 1)
})

pkgload::load_all(".", helpers = FALSE)
library(nlme)
library(survival)
library(JMbayes2)
library(splines)

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
    benchmark_stamp = as.character(benchmark_stamp)
  )
}

cfg <- benchmark_config(seed = 20260330L,
  n_rep = 100L,
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

if (!dir.exists(cfg$output_dir)) 
  dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)

# benchmark_stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
benchmark_stamp <- cfg$benchmark_stamp

base_name <- paste0(
  "joinme-vs-jmbayes2-", 
  cfg$n_markers, "markers-", 
  benchmark_stamp
)

output_rds <- file.path(
  cfg$output_dir,
  paste0(base_name, '.rds')
)

# Keep the longitudinal basis identical across simulation, joinme, and the
# per-marker JMbayes2 submodels so term-level recovery compares like with like.
benchmark_long_basis <- "splines::bs(time, knots = 1, degree = 2)"
benchmark_formula_long <- stats::as.formula(
  paste0(
    "y ~ 1 + ", benchmark_long_basis,
    " + (1 | id) + (1 + ", benchmark_long_basis,
    " + (1 | id) | marker)"
  )
)
benchmark_formula_jm_long <- function(outcome_var) {
  stats::as.formula(paste0(outcome_var, " ~ 1 + ", benchmark_long_basis))
}
benchmark_beta_long <- stats::setNames(
  c(1.0, 0.5, 0.6, 0.7),
  c("(Intercept)", paste0(benchmark_long_basis, 1:3))
)
benchmark_marker_re_params <- list(
  sd = c(1.0, 0.25, 0.2, 0.15),
  corr = diag(4)
)

safe_cor <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 2L) return(NA_real_)
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
    stop("sim$true_params$re_draws$marker must be a matrix-like object.")
  }
  if (nrow(marker_draws) != length(marker_levels) && ncol(marker_draws) == length(marker_levels)) {
    marker_draws <- t(marker_draws)
  }
  if (nrow(marker_draws) != length(marker_levels)) {
    stop(
      "Could not align sim$true_params$re_draws$marker with marker levels: ",
      nrow(marker_draws), " rows for ", length(marker_levels), " markers."
    )
  }
  if (ncol(marker_draws) != length(marker_terms)) {
    stop(
      "Could not align sim$true_params$re_draws$marker with marker terms: ",
      ncol(marker_draws), " columns for ", length(marker_terms), " terms."
    )
  }
  rownames(marker_draws) <- marker_levels
  colnames(marker_draws) <- marker_terms
  marker_draws
}

infer_marker_terms <- function(sim, marker_levels) {
  marker_draws <- sim$true_params$re_draws$marker
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
    beta_names <- names(sim$true_params$beta_long)
    if (length(beta_names) < term_count) {
      stop("Could not infer marker-level longitudinal terms from simulation truth.")
    }
    term_names <- beta_names[seq_len(term_count)]
  }
  as.character(term_names)
}

build_truth_long <- function(sim, marker_levels) {
  marker_terms <- infer_marker_terms(sim, marker_levels)
  marker_draws <- normalize_marker_draws(sim$true_params$re_draws$marker, marker_levels, marker_terms)
  fixed_truth <- as.numeric(sim$true_params$beta_long[marker_terms])
  fixed_truth[is.na(fixed_truth)] <- 0
  do.call(rbind, lapply(marker_levels, function(marker_label) {
    data.frame(
      marker = marker_label,
      term = marker_terms,
      truth = fixed_truth + as.numeric(marker_draws[marker_label, marker_terms, drop = TRUE]),
      stringsAsFactors = FALSE
    )
  }))
}

infer_survival_terms <- function(sim) {
  beta_event <- sim$true_params$beta_event
  term_names <- names(beta_event)
  if (!is.null(term_names) && length(term_names) > 0L && all(nzchar(term_names))) {
    return(as.character(term_names))
  }

  formula_event <- sim$true_params$formulaEvent %||% sim$truth$formulaEvent
  event_cols <- colnames(joinme:::.mm_event(formula_event, sim$dataEvent))
  as.character(event_cols %||% character(0))
}

build_truth_survival <- function(sim) {
  survival_terms <- infer_survival_terms(sim)
  truth_vals <- rep(0, length(survival_terms))
  names(truth_vals) <- survival_terms

  beta_event <- sim$true_params$beta_event
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
  assoc_coefs <- sim$true_params$assoc_coefs
  marker_weights_by_term <- null_or(sim$true_params$marker_weights_by_term, list())
  default_weights <- as.numeric(null_or(sim$true_params$marker_weights, rep(1, length(marker_levels))))
  cv_marker_weights <- as.numeric(null_or(marker_weights_by_term$cv_marker, default_weights))
  if (length(cv_marker_weights) != length(marker_levels)) {
    stop("cv_marker weights do not align with marker levels.")
  }
  cv_mean_coef <- as.numeric(null_or(assoc_coefs[["cv_mean"]], 0))
  cv_marker_coef <- as.numeric(null_or(assoc_coefs[["cv_marker"]], 0))
  data.frame(
    marker = marker_levels,
    truth = cv_mean_coef + cv_marker_coef * cv_marker_weights,
    stringsAsFactors = FALSE
  )
}

summarise_draw_matrix <- function(draw_matrix, id_col) {
  if (is.null(draw_matrix) || !is.matrix(draw_matrix) || ncol(draw_matrix) == 0L) {
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
  if (is.null(draw_array) || length(dim(draw_array)) != 3L || dim(draw_array)[3] == 0L) {
    return(data.frame())
  }
  summary_tbl <- joinme:::.assoc_summary_from_draw_array(draw_array, term_labels = labels, digits = 8)
  out <- data.frame(
    id = labels,
    estimate = summary_tbl$Estimate,
    std_error = summary_tbl$Est.Error,
    conf_low = summary_tbl$Q2.5,
    conf_high = summary_tbl$Q97.5,
    rhat = summary_tbl$Rhat,
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
  out <- draw_array[, , index, drop = FALSE]
  dim(out) <- dim(draw_array)[1:2]
  out
}

build_joinme_assoc_coef_draws <- function(fit_joinme, marker_levels, draws = NULL, seed = 1) {
  fit <- fit_joinme$fit
  all_vars <- tryCatch(posterior::variables(joinme:::.get_draws_obj(fit)), error = function(e) character(0))
  cv_mean_var <- if ("alpha_cv_mean" %in% all_vars) "alpha_cv_mean" else NULL
  cv_marker_var <- if ("alpha_cv_marker" %in% all_vars) "alpha_cv_marker" else NULL

  cv_mean_arr <- if (!is.null(cv_mean_var)) {
    joinme:::.get_draws_array(fit, variables = cv_mean_var, draws = draws, seed = seed)
  } else {
    NULL
  }
  cv_marker_arr <- if (!is.null(cv_marker_var)) {
    joinme:::.get_draws_array(fit, variables = cv_marker_var, draws = draws, seed = seed)
  } else {
    NULL
  }
  weight_arr <- if (!is.null(cv_marker_var)) {
    joinme:::.association_marker_weight_array(
      fit_joinme,
      term_key = "cv_marker",
      draws = draws,
      seed = seed,
      all_vars = all_vars
    )
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
      out_arr[, , marker_index] <- out_arr[, , marker_index] + cv_mean_mat
    }
  }
  if (!is.null(cv_marker_arr)) {
    if (is.null(weight_arr)) {
      stop("Could not recover marker weights for cv_marker association coefficient reconstruction.")
    }
    cv_marker_mat <- draw_array_slice(cv_marker_arr, 1L)
    for (marker_index in seq_len(n_markers)) {
      weight_mat <- draw_array_slice(weight_arr, marker_index)
      out_arr[, , marker_index] <- out_arr[, , marker_index] + cv_marker_mat * weight_mat
    }
  }

  joinme:::.assoc_matrix_from_draw_array(out_arr, marker_levels)
}

build_joinme_assoc_coef_array <- function(fit_joinme, marker_levels, draws = NULL, seed = 1) {
  fit <- fit_joinme$fit
  all_vars <- tryCatch(posterior::variables(joinme:::.get_draws_obj(fit)), error = function(e) character(0))
  cv_mean_var <- if ("alpha_cv_mean" %in% all_vars) "alpha_cv_mean" else NULL
  cv_marker_var <- if ("alpha_cv_marker" %in% all_vars) "alpha_cv_marker" else NULL

  cv_mean_arr <- if (!is.null(cv_mean_var)) {
    joinme:::.get_draws_array(fit, variables = cv_mean_var, draws = draws, seed = seed)
  } else {
    NULL
  }
  cv_marker_arr <- if (!is.null(cv_marker_var)) {
    joinme:::.get_draws_array(fit, variables = cv_marker_var, draws = draws, seed = seed)
  } else {
    NULL
  }
  weight_arr <- if (!is.null(cv_marker_var)) {
    joinme:::.association_marker_weight_array(
      fit_joinme,
      term_key = "cv_marker",
      draws = draws,
      seed = seed,
      all_vars = all_vars
    )
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
      out_arr[, , marker_index] <- out_arr[, , marker_index] + cv_mean_mat
    }
  }
  if (!is.null(cv_marker_arr)) {
    if (is.null(weight_arr)) {
      stop("Could not recover marker weights for cv_marker association coefficient reconstruction.")
    }
    cv_marker_mat <- draw_array_slice(cv_marker_arr, 1L)
    for (marker_index in seq_len(n_markers)) {
      weight_mat <- draw_array_slice(weight_arr, marker_index)
      out_arr[, , marker_index] <- out_arr[, , marker_index] + cv_marker_mat * weight_mat
    }
  }

  out_arr
}

parse_assoc_marker <- function(x, fallback = NULL) {
  if (is.null(x)) return(fallback)
  out <- as.character(x)
  out <- gsub("^association\\.?", "", out)
  out <- gsub("^value\\(", "", out)
  out <- gsub("\\)$", "", out)
  out <- gsub("^y_", "", out)
  out <- gsub(".*y_", "", out)
  out[nchar(out) == 0L] <- fallback[nchar(out) == 0L]
  out
}

combine_truth_and_estimates <- function(truth_df, est_df, by_cols, replication, engine) {
  out <- merge(truth_df, est_df, by = by_cols, all.x = TRUE, sort = FALSE)
  out$replication <- replication
  out$engine <- engine
  out$error <- out$estimate - out$truth
  out
}

extract_joinme_long <- function(fit_joinme, truth_long, replication) {
  coef_summary <- coef(fit_joinme, summary = TRUE)
  est_df <- coef_summary$formulaLong$marker[, c("marker", "term", "Estimate", "Est.Error", "Q2.5", "Q97.5", "Rhat")]
  names(est_df) <- c("marker", "term", "estimate", "std_error", "conf_low", "conf_high", "rhat")
  combine_truth_and_estimates(truth_long, est_df, c("marker", "term"), replication, "joinme")
}

extract_joinme_survival <- function(fit_joinme, truth_survival, replication) {
  coef_summary <- coef(fit_joinme, summary = TRUE)
  est_df <- coef_summary$formulaEvent[, c("term", "Estimate", "Est.Error", "Q2.5", "Q97.5", "Rhat")]
  names(est_df) <- c("term", "estimate", "std_error", "conf_low", "conf_high", "rhat")
  combine_truth_and_estimates(truth_survival, est_df, "term", replication, "joinme")
}

extract_joinme_assoc <- function(fit_joinme, truth_assoc, marker_levels, replication, seed) {
  draw_arr <- build_joinme_assoc_coef_array(fit_joinme, marker_levels, seed = seed)
  est_df <- summarise_draw_array_with_rhat(draw_arr, "marker", marker_levels)
  combine_truth_and_estimates(truth_assoc, est_df, "marker", replication, "joinme")
}

extract_jm_long <- function(jm_summary, truth_long, marker_levels, replication) {
  truth_terms <- unique(as.character(truth_long$term))
  outcome_tabs <- grep("^Outcome", names(jm_summary), value = TRUE)
  est_df <- do.call(rbind, lapply(seq_along(outcome_tabs), function(idx) {
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
      rhat = if ("Rhat" %in% names(tab)) as.numeric(tab[, "Rhat"]) else NA_real_,
      stringsAsFactors = FALSE
    )
    out[out$term %in% truth_terms, , drop = FALSE]
  }))
  combine_truth_and_estimates(truth_long, est_df, c("marker", "term"), replication, "JMbayes2")
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
  combine_truth_and_estimates(truth_survival, est_df, "term", replication, "JMbayes2")
}

extract_jm_assoc <- function(jm_summary, truth_assoc, truth_survival, marker_levels, replication) {
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
    marker = parse_assoc_marker(rownames(tab)[keep], fallback = marker_levels[seq_len(sum(keep))]),
    estimate = as.numeric(tab[keep, "Mean"]),
    std_error = as.numeric(tab[keep, "StDev"]),
    conf_low = as.numeric(tab[keep, "2.5%"]),
    conf_high = as.numeric(tab[keep, "97.5%"]),
    rhat = if ("Rhat" %in% names(tab)) as.numeric(tab[keep, "Rhat"]) else NA_real_,
    stringsAsFactors = FALSE
  )
  combine_truth_and_estimates(truth_assoc, est_df, "marker", replication, "JMbayes2")
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
      asymptotic_std_error = if (length(est_vals) > 1L) stats::sd(est_vals) else NA_real_,
      average_rhat = if (length(rhat_vals)) mean(rhat_vals) else NA_real_,
      max_rhat = if (length(rhat_vals)) max(rhat_vals) else NA_real_,
      bad_rhat_rate = 100 * sum(bad_rhat) / total_reps,
      coverage_rate = 100 * sum(covered) / total_reps,
      coverage_rate_available = if (sum(available) > 0L) 100 * mean(covered[available]) else NA_real_,
      cor = if (length(est_vals) > 1L) safe_cor(block$truth[available], est_vals) else NA_real_,
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
  keep_cols <- c(key_cols, "engine", "average_rhat", "max_rhat", "bad_rhat_rate")
  wide <- summary_df[, keep_cols, drop = FALSE]
  split_df <- split(wide, wide$engine)
  joinme_df <- split_df$joinme
  jm_df <- split_df$JMbayes2
  if (is.null(joinme_df) || is.null(jm_df)) {
    return(data.frame())
  }
  names(joinme_df)[names(joinme_df) %in% c("average_rhat", "max_rhat", "bad_rhat_rate")] <- paste0(c("average_rhat", "max_rhat", "bad_rhat_rate"), "_joinme")
  names(jm_df)[names(jm_df) %in% c("average_rhat", "max_rhat", "bad_rhat_rate")] <- paste0(c("average_rhat", "max_rhat", "bad_rhat_rate"), "_JMbayes2")
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
  simulate_joinme(
    formulaLong = benchmark_formula_long,
    formulaEvent = survival::Surv(time, event) ~ x1 + x2,
    n_id = cfg$n_id,
    families = rep("gaussian", cfg$n_markers),
    times_obs = seq(0, 5, length.out = 5),
    n_obs_per_marker_per_id = 5,
    assoc = c("cv_mean", "cv_marker"),
    assoc_coefs = c(cv_mean = 1, cv_marker = 1),
    marker_weights = rep(1, cfg$n_markers),
    fixed_marker_weights = TRUE,
    shared_marker_weights = TRUE,
    beta_long = benchmark_beta_long,
    beta_event = c(-0.5, 0.25),
    baseline_hazard = list(type = "weibull", shape = 0.8, scale = 12),
    re_params = list(
      id = list(sd = c(1)),
      marker = benchmark_marker_re_params,
      id_marker_cov = list(
        latent = list(sd = 1),
        alpha = c(-0.5),
        lambda = c(1)
      )
    ),
    use_mirai = cfg$use_mirai,
    n_workers = cfg$n_workers,
    seed = seed
  )
}

runtime_records <- list()
longitudinal_records <- list()
survival_records <- list()
association_records <- list()
error_records <- list()

cat("Running ", cfg$n_rep, " benchmark replications with ", cfg$n_markers, " biomarkers...\n", sep = "")

# create tmp directory
dir.create('sim', showWarnings = FALSE)
for (replication in seq_len(cfg$n_rep)) {
  rep_seed <- cfg$seed + replication - 1L
  cat("[", replication, "/", cfg$n_rep, "] seed=", rep_seed, "\n", sep = "")

  # Check for existing simulation result for this replication
  if (file.exists(
    file.path("sim", 
      paste0(base_name, '-rep', replication, '.rds'))
  )) {
    cat("Found existing simulation for replication ", replication, "\n", sep = "")
    fit_all <- readRDS(
      file.path("sim", 
        paste0(base_name, '-rep', replication, '.rds'))
    )
    joinme_result <- fit_all$joinme_result

    runtime_records[[length(runtime_records) + 1L]] <- data.frame(
      replication = replication,
      engine = "joinme",
      elapsed_seconds = joinme_result$elapsed,
      success = is.null(joinme_result$error),
      stringsAsFactors = FALSE
    )

    fit_joinme <- joinme_result$value
    longitudinal_records[[length(longitudinal_records) + 1L]] <- extract_joinme_long(fit_joinme, truth_long, replication)
    survival_records[[length(survival_records) + 1L]] <- extract_joinme_survival(fit_joinme, truth_survival, replication)
    association_records[[length(association_records) + 1L]] <- extract_joinme_assoc(fit_joinme, truth_assoc, marker_levels, replication, seed = rep_seed)

    jm_lme_result <- fit_all$jm_lme_result
    runtime_records[[length(runtime_records) + 1L]] <- data.frame(
      replication = replication,
      engine = "JMbayes2_lme",
      elapsed_seconds = jm_lme_result$elapsed,
      success = is.null(jm_lme_result$error),
      stringsAsFactors = FALSE
    )

    jm_surv_result <- fit_all$jm_surv_result
    runtime_records[[length(runtime_records) + 1L]] <- data.frame(
      replication = replication,
      engine = "JMbayes2_survival",
      elapsed_seconds = jm_surv_result$elapsed,
      success = is.null(jm_surv_result$error),
      stringsAsFactors = FALSE
    )

    jm_joint_result <- fit_all$jm_joint_result
    runtime_records[[length(runtime_records) + 1L]] <- data.frame(
      replication = replication,
      engine = "JMbayes2_joint",
      elapsed_seconds = jm_joint_result$elapsed,
      success = is.null(jm_joint_result$error),
      stringsAsFactors = FALSE
    )

    runtime_records[[length(runtime_records) + 1L]] <- data.frame(
      replication = replication,
      engine = "JMbayes2_total",
      elapsed_seconds = jm_lme_result$elapsed + jm_surv_result$elapsed + jm_joint_result$elapsed,
      success = is.null(jm_joint_result$error),
      stringsAsFactors = FALSE
    )

    jm_summary <- summary(jm_joint_result$value)
    longitudinal_records[[length(longitudinal_records) + 1L]] <- extract_jm_long(jm_summary, truth_long, marker_levels, replication)
    survival_records[[length(survival_records) + 1L]] <- extract_jm_survival(jm_summary, truth_survival, replication)
    association_records[[length(association_records) + 1L]] <- extract_jm_assoc(jm_summary, truth_assoc, truth_survival, marker_levels, replication)

    next
  }


  sim <- simulate_replication(rep_seed)
  marker_levels <- as.character(sim$marker_info$names)
  truth_long <- build_truth_long(sim, marker_levels)
  truth_survival <- build_truth_survival(sim)
  truth_assoc <- build_truth_assoc(sim, marker_levels)

  joinme_result <- time_try(
    joinme(
      formulaLong = benchmark_formula_long,
      dataLong = sim$dataLong,
      formulaEvent = survival::Surv(time, event) ~ x1 + x2,
      dataEvent = sim$dataEvent,
      assoc = c("cv_mean", "cv_marker"),
      families = rep("gaussian", cfg$n_markers),
      marker_weights = rep(1, cfg$n_markers),
      fixed_marker_weights = TRUE,
      shared_marker_weights = TRUE,
      control = list(
        engine = "cmdstanr",
        chains = 2,
        parallel_chains = 2,
        threads_per_chain = 6,
        iter_warmup = cfg$iter_warmup_joinme,
        iter_sampling = cfg$iter_sampling_joinme,
        refresh = 0,
        show_messages = FALSE,
        init = 1,
        seed = rep_seed,
        adapt_delta = 0.72,
        max_treedepth = 11
      )
    )
  )

  runtime_records[[length(runtime_records) + 1L]] <- data.frame(
    replication = replication,
    engine = "joinme",
    elapsed_seconds = joinme_result$elapsed,
    success = is.null(joinme_result$error),
    stringsAsFactors = FALSE
  )

  if (is.null(joinme_result$error)) {
    fit_joinme <- joinme_result$value
    longitudinal_records[[length(longitudinal_records) + 1L]] <- extract_joinme_long(fit_joinme, truth_long, replication)
    survival_records[[length(survival_records) + 1L]] <- extract_joinme_survival(fit_joinme, truth_survival, replication)
    association_records[[length(association_records) + 1L]] <- extract_joinme_assoc(fit_joinme, truth_assoc, marker_levels, replication, seed = rep_seed)
  } else {
    error_records[[length(error_records) + 1L]] <- data.frame(
      replication = replication,
      engine = "joinme",
      stage = "fit",
      error = joinme_result$error,
      stringsAsFactors = FALSE
    )
  }

  wide_long <- reshape(
    sim$dataLong[, c("id", "time", "marker", "y")],
    idvar = c("id", "time"),
    timevar = "marker",
    direction = "wide"
  )
  wide_long <- wide_long[order(wide_long$id, wide_long$time), ]
  names(wide_long) <- sub("^y\\.", "y_", names(wide_long))
  outcome_vars <- paste0("y_", marker_levels)

  jm_lme_result <- time_try({
    fits <- lapply(outcome_vars, function(outcome_var) {
      nlme::lme(
        benchmark_formula_jm_long(outcome_var),
        random = ~ 1 | id,
        data = wide_long,
        na.action = na.exclude,
        control = nlme::lmeControl(msMaxIter = 100, msMaxEval = 200, pnlsMaxIter = 25, returnObject = TRUE)
      )
    })
    names(fits) <- outcome_vars
    fits
  })

  runtime_records[[length(runtime_records) + 1L]] <- data.frame(
    replication = replication,
    engine = "JMbayes2_lme",
    elapsed_seconds = jm_lme_result$elapsed,
    success = is.null(jm_lme_result$error),
    stringsAsFactors = FALSE
  )

  if (!is.null(jm_lme_result$error)) {
    error_records[[length(error_records) + 1L]] <- data.frame(
      replication = replication,
      engine = "JMbayes2",
      stage = "lme",
      error = jm_lme_result$error,
      stringsAsFactors = FALSE
    )
  }

  jm_surv_result <- time_try(
    survival::coxph(survival::Surv(time, event) ~ x1 + x2, data = sim$dataEvent, x = TRUE)
  )

  runtime_records[[length(runtime_records) + 1L]] <- data.frame(
    replication = replication,
    engine = "JMbayes2_survival",
    elapsed_seconds = jm_surv_result$elapsed,
    success = is.null(jm_surv_result$error),
    stringsAsFactors = FALSE
  )

  if (!is.null(jm_surv_result$error)) {
    error_records[[length(error_records) + 1L]] <- data.frame(
      replication = replication,
      engine = "JMbayes2",
      stage = "survival",
      error = jm_surv_result$error,
      stringsAsFactors = FALSE
    )
  }

  jm_joint_result <- if (is.null(jm_lme_result$error) && is.null(jm_surv_result$error)) {
    time_try(
      JMbayes2::jm(
        jm_surv_result$value,
        jm_lme_result$value,
        time_var = "time",
        n_iter = cfg$n_iter_jm,
        n_burnin = cfg$n_burnin_jm,
        n_chains = cfg$n_chains_jm,
        n_thin = cfg$n_thin_jm,
        cores = cfg$n_chains_jm
      )
    )
  } else {
    list(value = NULL, error = "Skipped because JMbayes2 preprocessing failed.", elapsed = 0)
  }

  runtime_records[[length(runtime_records) + 1L]] <- data.frame(
    replication = replication,
    engine = "JMbayes2_joint",
    elapsed_seconds = jm_joint_result$elapsed,
    success = is.null(jm_joint_result$error),
    stringsAsFactors = FALSE
  )
  runtime_records[[length(runtime_records) + 1L]] <- data.frame(
    replication = replication,
    engine = "JMbayes2_total",
    elapsed_seconds = jm_lme_result$elapsed + jm_surv_result$elapsed + jm_joint_result$elapsed,
    success = is.null(jm_joint_result$error),
    stringsAsFactors = FALSE
  )

  if (!is.null(jm_joint_result$error)) {
    error_records[[length(error_records) + 1L]] <- data.frame(
      replication = replication,
      engine = "JMbayes2",
      stage = "joint",
      error = jm_joint_result$error,
      stringsAsFactors = FALSE
    )
  } else {
    jm_summary <- summary(jm_joint_result$value)
    longitudinal_records[[length(longitudinal_records) + 1L]] <- extract_jm_long(jm_summary, truth_long, marker_levels, replication)
    survival_records[[length(survival_records) + 1L]] <- extract_jm_survival(jm_summary, truth_survival, replication)
    association_records[[length(association_records) + 1L]] <- extract_jm_assoc(jm_summary, truth_assoc, truth_survival, marker_levels, replication)
  }

  # Save results to rds
  saveRDS(list(
    joinme_result = joinme_result, 
    jm_lme_result = jm_lme_result, 
    jm_surv_result = jm_surv_result,
    jm_joint_result = jm_joint_result),
    file.path("sim", paste0(base_name, '-rep', replication, '.rds'))
  )
}

runtime_records <- do.call(rbind, runtime_records)
longitudinal_records <- bind_rows_or_empty(longitudinal_records)
survival_records <- bind_rows_or_empty(survival_records)
association_records <- bind_rows_or_empty(association_records)
error_records <- bind_rows_or_empty(error_records)

runtime_summary <- summarise_runtime_metrics(runtime_records, total_reps = cfg$n_rep)
longitudinal_summary <- summarise_replication_metrics(longitudinal_records, c("marker", "term"), total_reps = cfg$n_rep)
survival_summary <- summarise_replication_metrics(survival_records, c("term"), total_reps = cfg$n_rep)
association_summary <- summarise_replication_metrics(association_records, c("marker"), total_reps = cfg$n_rep)

longitudinal_rhat_comparison <- build_rhat_comparison(longitudinal_summary, c("marker", "term"))
survival_rhat_comparison <- build_rhat_comparison(survival_summary, c("term"))
association_rhat_comparison <- build_rhat_comparison(association_summary, c("marker"))

results <- list(
  config = cfg,
  runtime_replications = runtime_records,
  runtime_summary = runtime_summary,
  longitudinal_replications = longitudinal_records,
  longitudinal_summary = longitudinal_summary,
  survival_replications = survival_records,
  survival_summary = survival_summary,
  association_replications = association_records,
  association_summary = association_summary,
  rhat_comparison = list(
    longitudinal = longitudinal_rhat_comparison,
    survival = survival_rhat_comparison,
    association = association_rhat_comparison
  ),
  errors = error_records
)

saveRDS(results, output_rds)

cat("Benchmark replications complete.\n")
cat("Results saved to: ", output_rds, "\n", sep = "")
cat("\nRuntime summary\n")
print(runtime_summary, row.names = FALSE)
cat("\nLongitudinal summary (first 12 rows)\n")
print(utils::head(longitudinal_summary, 12), row.names = FALSE)
cat("\nSurvival summary\n")
print(survival_summary, row.names = FALSE)
cat("\nAssociation summary (first 12 rows)\n")
print(utils::head(association_summary, 12), row.names = FALSE)
if (nrow(longitudinal_rhat_comparison)) {
  cat("\nLongitudinal Rhat comparison (first 12 rows)\n")
  print(utils::head(longitudinal_rhat_comparison, 12), row.names = FALSE)
}
if (nrow(survival_rhat_comparison)) {
  cat("\nSurvival Rhat comparison\n")
  print(survival_rhat_comparison, row.names = FALSE)
}
if (nrow(association_rhat_comparison)) {
  cat("\nAssociation Rhat comparison (first 12 rows)\n")
  print(utils::head(association_rhat_comparison, 12), row.names = FALSE)
}
if (nrow(error_records)) {
  cat("\nErrors (first 12 rows)\n")
  print(utils::head(error_records, 12), row.names = FALSE)
}
