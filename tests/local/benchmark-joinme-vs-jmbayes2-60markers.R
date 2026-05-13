options(error = function() {
  traceback(2)
  quit(status = 1)
})

pkgload::load_all(".", helpers = FALSE)
library(nlme)
library(survival)
library(JMbayes2)

benchmark_config <- function(
  seed = 20260329L,
  n_id = 50L,
  n_markers = 60L,
  iter_warmup_joinme = 800L,
  iter_sampling_joinme = 1000L,
  n_iter_jm = 6000L,
  n_burnin_jm = 2000L,
  n_chains_jm = 2L,
  output_dir = ".artifacts"
) {
  list(
    seed = as.integer(seed),
    n_id = as.integer(n_id),
    n_markers = as.integer(n_markers),
    iter_warmup_joinme = as.integer(iter_warmup_joinme),
    iter_sampling_joinme = as.integer(iter_sampling_joinme),
    n_iter_jm = as.integer(n_iter_jm),
    n_burnin_jm = as.integer(n_burnin_jm),
    n_chains_jm = as.integer(n_chains_jm),
    output_dir = as.character(output_dir)
  )
}

cfg <- benchmark_config(
  seed = 20260329L,
  n_id = 50L,
  n_markers = 60L,
  iter_warmup_joinme = 500L,
  iter_sampling_joinme = 1000L,
  n_iter_jm = 4000L,
  n_burnin_jm = 1000L,
  n_chains_jm = 2L,
  output_dir = ".artifacts"
)

if (!dir.exists(cfg$output_dir)) dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)

benchmark_stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
output_rds <- file.path(cfg$output_dir, paste0("joinme-vs-jmbayes2-", cfg$n_markers, "markers-", benchmark_stamp, ".rds"))

safe_cor <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 2L) return(NA_real_)
  stats::cor(x[keep], y[keep])
}

null_or <- function(x, default) {
  if (is.null(x)) default else x
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

summarise_recovery <- function(df, group_cols) {
  if (!nrow(df)) return(data.frame())
  split_idx <- if (length(group_cols) == 0L) {
    list(all = seq_len(nrow(df)))
  } else {
    split(seq_len(nrow(df)), interaction(df[group_cols], drop = TRUE, lex.order = TRUE))
  }
  out <- lapply(split_idx, function(idx) {
    block <- df[idx, , drop = FALSE]
    data.frame(
      if (length(group_cols)) block[1, group_cols, drop = FALSE] else data.frame(scope = "all", stringsAsFactors = FALSE),
      n = nrow(block),
      truth_mean = mean(block$truth),
      estimate_mean = mean(block$estimate),
      bias = mean(block$error),
      mae = mean(abs(block$error)),
      rmse = sqrt(mean(block$error^2)),
      cor = safe_cor(block$truth, block$estimate),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
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
  weighted_divisor <- if (length(marker_levels) > 0L) length(marker_levels) else 1L

  data.frame(
    marker = marker_levels,
    truth = cv_mean_coef + cv_marker_coef * cv_marker_weights / weighted_divisor,
    stringsAsFactors = FALSE
  )
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

build_joinme_assoc_draws <- function(fit_joinme, marker_levels) {
  assoc_draws <- posterior_assoc(fit_joinme, summary = FALSE)
  n_markers <- length(marker_levels)

  cv_mean_draws <- assoc_draws[["cv_mean"]]
  cv_marker_draws <- assoc_draws[["cv_marker"]]

  n_draws <- NULL
  if (!is.null(cv_mean_draws)) {
    n_draws <- nrow(cv_mean_draws)
  }
  if (is.null(n_draws) && !is.null(cv_marker_draws)) {
    n_draws <- nrow(cv_marker_draws)
  }
  if (is.null(n_draws)) {
    return(NULL)
  }

  out <- matrix(0, nrow = n_draws, ncol = n_markers)
  colnames(out) <- marker_levels

  if (!is.null(cv_mean_draws)) {
    out <- out + matrix(as.numeric(cv_mean_draws[, 1]), nrow = n_draws, ncol = n_markers)
  }
  if (!is.null(cv_marker_draws)) {
    out <- out + as.matrix(cv_marker_draws[, marker_levels, drop = FALSE])
  }

  out
}

cat("[1/6] Simulating dataset with ", cfg$n_markers, " biomarkers...\n", sep = "")
set.seed(cfg$seed)
sim <- simulate_joinme(
  formulaLong = y ~ 1 + time + (1 | id) + (1 + time + (1 | id) | marker),
  formulaEvent = survival::Surv(time, event) ~ x1 + x2,
  n_id = cfg$n_id,
  families = rep("gaussian", cfg$n_markers),
  times_obs = seq(0, 4, length.out = 5),
  n_obs_per_marker_per_id = 5,
  assoc = c("cv_mean", "cv_marker"),
  assoc_coefs = c(cv_mean = 1, cv_marker = 0.5),
  marker_weights = rep(1, cfg$n_markers),
  shared_marker_weights = TRUE,
  beta_long = c(-1, 0.5),
  beta_event = c(-0.5, 0.25),
  re_params = list(
    id = list(sd = c(1)),
    marker = list(sd = c(1, 0.25), corr = matrix(c(1, -0.25, -0.25, 1), ncol=2)),
    id_marker_cov = list(
        latent = list(
            sd = 1
        ),
        alpha = c(-0.5),
        lambda = c(1)
    )
  ),
  use_mirai = TRUE,
  n_workers = 10,
  seed = cfg$seed
)

marker_levels <- as.character(sim$marker_info$names)
truth_survival <- build_truth_survival(sim)
truth_assoc <- build_truth_assoc(sim, marker_levels)

cat("[2/6] Fitting joinme...\n")
joinme_elapsed <- system.time({
  fit_joinme <- joinme(
    formulaLong = y ~ 1 + time + (1 | id) + (1 + time + (1 | id) | marker),
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
      init = 1,
      seed = cfg$seed,
      adapt_delta = 0.7,
      max_treedepth = 12
    )
  )
})

joinme_coef <- coef(fit_joinme, summary = TRUE)
marker_terms <- unique(as.character(joinme_coef$formulaLong$marker$term))
marker_draws <- normalize_marker_draws(
  marker_draws = sim$true_params$re_draws$marker,
  marker_levels = marker_levels,
  marker_terms = marker_terms
)
fixed_truth <- as.numeric(sim$true_params$beta_long[marker_terms])
fixed_truth[is.na(fixed_truth)] <- 0
truth_long <- do.call(rbind, lapply(marker_levels, function(marker_label) {
  data.frame(
    marker = marker_label,
    term = marker_terms,
    truth = fixed_truth + as.numeric(marker_draws[marker_label, marker_terms, drop = TRUE]),
    stringsAsFactors = FALSE
  )
}))
joinme_long_est <- joinme_coef$formulaLong$marker[, c("marker", "term", "Estimate")]
names(joinme_long_est)[names(joinme_long_est) == "Estimate"] <- "estimate"
joinme_long_cmp <- merge(truth_long, joinme_long_est, by = c("marker", "term"), all.x = TRUE)
joinme_long_cmp$error <- joinme_long_cmp$estimate - joinme_long_cmp$truth

joinme_surv_est <- joinme_coef$formulaEvent[, c("term", "Estimate")]
names(joinme_surv_est)[names(joinme_surv_est) == "Estimate"] <- "estimate"
joinme_surv_cmp <- merge(truth_survival, joinme_surv_est, by = "term", all.x = TRUE)
joinme_surv_cmp$error <- joinme_surv_cmp$estimate - joinme_surv_cmp$truth

joinme_assoc_draws <- build_joinme_assoc_draws(fit_joinme, marker_levels)
joinme_assoc_est <- data.frame(
  marker = colnames(joinme_assoc_draws),
  estimate = colMeans(joinme_assoc_draws),
  stringsAsFactors = FALSE
)
joinme_assoc_cmp <- merge(truth_assoc, joinme_assoc_est, by = "marker", all.x = TRUE)
joinme_assoc_cmp$error <- joinme_assoc_cmp$estimate - joinme_assoc_cmp$truth

cat("[3/6] Reshaping data for JMbayes2...\n")
wide_long <- reshape(
  sim$dataLong[, c("id", "time", "x1", "marker", "y")],
  idvar = c("id", "time", "x1"),
  timevar = "marker",
  direction = "wide"
)
wide_long <- wide_long[order(wide_long$id, wide_long$time), ]
names(wide_long) <- sub("^y\\.", "y_", names(wide_long))
outcome_vars <- paste0("y_", marker_levels)

cat("[4/6] Fitting JMbayes2 longitudinal submodels...\n")
jm_lme_elapsed <- system.time({
  jm_mixed <- lapply(outcome_vars, function(outcome_var) {
    nlme::lme(
      stats::as.formula(paste0(outcome_var, " ~ time + x1")),
      random = ~ 1 | id,
      data = wide_long,
      na.action = na.exclude,
      control = nlme::lmeControl(msMaxIter = 100, msMaxEval = 200, pnlsMaxIter = 25, returnObject = TRUE)
    )
  })
  names(jm_mixed) <- outcome_vars
})

cat("[5/6] Fitting JMbayes2 survival and joint model...\n")
jm_fit_error <- NULL
jm_surv_elapsed <- system.time({
  jm_surv <- survival::coxph(survival::Surv(time, event) ~ x1 + x2, data = sim$dataEvent, x = TRUE)
})
jm_elapsed <- system.time({
  jm_fit <- tryCatch(
    JMbayes2::jm(
      jm_surv,
      jm_mixed,
      time_var = "time",
    #   which_independent = "all",
    #   priors = list(penalty_alphas = "ridge"),
      n_iter = cfg$n_iter_jm,
      n_burnin = cfg$n_burnin_jm,
      n_chains = cfg$n_chains_jm,
      cores = 2
    ),
    error = function(e) {
      jm_fit_error <<- conditionMessage(e)
      NULL
    }
  )
})

jm_long_cmp <- data.frame()
jm_surv_cmp <- data.frame()
jm_assoc_cmp <- data.frame()

if (is.null(jm_fit)) {
  cat("JMbayes2 fit failed: ", jm_fit_error, "\n", sep = "")
} else {
  jm_fixef <- fixef(jm_fit)
  jm_long_est <- do.call(rbind, lapply(names(jm_fixef), function(outcome_name) {
    coeffs <- jm_fixef[[outcome_name]]
    data.frame(
      marker = sub("^y_", "", outcome_name),
      term = names(coeffs),
      estimate = as.numeric(coeffs),
      stringsAsFactors = FALSE
    )
  }))
  jm_long_est <- jm_long_est[jm_long_est$term %in% marker_terms, , drop = FALSE]
  jm_long_cmp <- merge(truth_long, jm_long_est, by = c("marker", "term"), all.x = TRUE)
  jm_long_cmp$error <- jm_long_cmp$estimate - jm_long_cmp$truth

  jm_coef <- coef(jm_fit)
  jm_gamma <- jm_coef$gammas %||% numeric(0)
  jm_assoc <- jm_coef$association %||% numeric(0)

  if (length(jm_gamma) > 0L) {
    jm_surv_est <- data.frame(term = names(jm_gamma), estimate = as.numeric(jm_gamma), stringsAsFactors = FALSE)
    jm_surv_cmp <- merge(truth_survival, jm_surv_est, by = "term", all.x = TRUE)
    jm_surv_cmp$error <- jm_surv_cmp$estimate - jm_surv_cmp$truth
  }

  if (length(jm_assoc) > 0L) {
    assoc_markers <- parse_assoc_marker(names(jm_assoc), fallback = marker_levels[seq_along(jm_assoc)])
    jm_assoc_est <- data.frame(marker = assoc_markers, estimate = as.numeric(jm_assoc), stringsAsFactors = FALSE)
    jm_assoc_cmp <- merge(truth_assoc, jm_assoc_est, by = "marker", all.x = TRUE)
    jm_assoc_cmp$error <- jm_assoc_cmp$estimate - jm_assoc_cmp$truth
  }
}

runtime <- data.frame(
  engine = c("joinme", "JMbayes2_lme", "JMbayes2_survival", "JMbayes2_joint", "JMbayes2_total"),
  elapsed_seconds = c(
    unname(joinme_elapsed[["elapsed"]]),
    unname(jm_lme_elapsed[["elapsed"]]),
    unname(jm_surv_elapsed[["elapsed"]]),
    unname(jm_elapsed[["elapsed"]]),
    unname(jm_lme_elapsed[["elapsed"]] + jm_surv_elapsed[["elapsed"]] + jm_elapsed[["elapsed"]])
  ),
  stringsAsFactors = FALSE
)

recovery <- list(
  longitudinal = do.call(rbind, Filter(Negate(is.null), list(
    transform(summarise_recovery(joinme_long_cmp, c("term")), engine = "joinme"),
    if (nrow(jm_long_cmp)) transform(summarise_recovery(jm_long_cmp, c("term")), engine = "JMbayes2") else NULL
  ))),
  survival = do.call(rbind, Filter(Negate(is.null), list(
    transform(summarise_recovery(joinme_surv_cmp, c("term")), engine = "joinme"),
    if (nrow(jm_surv_cmp)) transform(summarise_recovery(jm_surv_cmp, c("term")), engine = "JMbayes2") else NULL
  ))),
  association = do.call(rbind, Filter(Negate(is.null), list(
    transform(summarise_recovery(joinme_assoc_cmp, character(0)), engine = "joinme"),
    if (nrow(jm_assoc_cmp)) transform(summarise_recovery(jm_assoc_cmp, character(0)), engine = "JMbayes2") else NULL
  )))
)

results <- list(
  config = cfg,
  runtime = runtime,
  recovery = recovery,
  truth = list(longitudinal = truth_long, survival = truth_survival, association = truth_assoc),
  estimates = list(
    joinme = list(longitudinal = joinme_long_cmp, survival = joinme_surv_cmp, association = joinme_assoc_cmp),
    JMbayes2 = list(longitudinal = jm_long_cmp, survival = jm_surv_cmp, association = jm_assoc_cmp)
  ),
  jm_fit_error = jm_fit_error
)

saveRDS(results, output_rds)

cat("[6/6] Benchmark complete.\n")
cat("Results saved to: ", output_rds, "\n", sep = "")
cat("\nRuntime (seconds)\n")
print(runtime, row.names = FALSE)
cat("\nLongitudinal recovery\n")
print(recovery$longitudinal, row.names = FALSE)
cat("\nSurvival recovery\n")
print(recovery$survival, row.names = FALSE)
cat("\nAssociation recovery\n")
print(recovery$association, row.names = FALSE)
if (!is.null(jm_fit_error)) {
  cat("\nJMbayes2 error\n")
  cat(jm_fit_error, "\n")
}
