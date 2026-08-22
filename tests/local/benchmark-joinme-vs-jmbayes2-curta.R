options(error = function() {
  traceback(2)
  quit(status = 1)
})

pkgload::load_all(".", helpers = FALSE)
library(nlme)
library(survival)
library(JMbayes2)
library(splines)

source('benchmark-helpers.R')
runtime_records <- list()
longitudinal_records <- list()
survival_records <- list()
association_records <- list()
error_records <- list()

cat(
  "Running ",
  cfg$n_rep,
  " benchmark replications with ",
  cfg$n_markers,
  " biomarkers...\n",
  sep = ""
)

# create tmp directory
dir.create(".artifacts/records", showWarnings = FALSE)
replication <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = 1L))

rep_seed <- cfg$seed + replication - 1L
cat("[", replication, "/", cfg$n_rep, "] seed=", rep_seed, "\n", sep = "")

# Check for existing simulation result for this replication
if (
  file.exists(
    file.path(".artifacts/records", paste0(base_name, '-rep', replication, '.rds'))
  )
) {
  cat("Found existing simulation for replication ", replication, "\n", sep = "")
  records <- readRDS(
    file.path(".artifacts/records", paste0(base_name, '-rep', replication, '.rds'))
  )

  # runtime_records <- append(runtime_records, records$runtime_records)
  # longitudinal_records <- append(
  #   longitudinal_records,
  #   records$longitudinal_records
  # )
  # survival_records <- append(survival_records, records$survival_records)
  # association_records <- append(
  #   association_records,
  #   records$association_records
  # )
  # error_records <- append(error_records, records$error_records)

  # next
}

cfg$n_id = 300
sim <- simulate_replication(rep_seed)
marker_levels <- as.character(sim$marker_info$names)
marker_weight_truth <- benchmark_marker_weight_truth(sim, term = "cv_marker")
print(marker_weight_truth, row.names = FALSE)
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
    priors = jm_prior(
      intercept = prior_student_t(df=3, mu = 0, scale = 3),
      slope = prior_student_t(df=3, mu = 0, scale = 3),

      marker_weights = list(
        # This prior is deliberately independent of the randomly generated
        # population marker-weight mean retained in sim$truth.
        intercept = prior_student_t(df = 3, mu = 0, scale = 3),
        offset = rep(0, cfg$n_markers),
        # Estimate marker-weight departures on a unit normal scale.
        family = "normal"
      )
    ),
    control = list(
      engine = "cmdstanr",
      chains = 2,
      parallel_chains = 2,
      threads_per_chain = 4,
      iter_warmup = cfg$iter_warmup_joinme,
      iter_sampling = cfg$iter_sampling_joinme,
      refresh = 200,
      show_messages = TRUE,
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
  longitudinal_records[[
    length(longitudinal_records) + 1L
  ]] <- extract_joinme_long(fit_joinme, truth_long, replication)
  survival_records[[length(survival_records) + 1L]] <- extract_joinme_survival(
    fit_joinme,
    truth_survival,
    replication
  )
  association_records[[
    length(association_records) + 1L
  ]] <- extract_joinme_assoc(
    fit_joinme,
    truth_assoc,
    marker_levels,
    replication,
    seed = rep_seed
  )
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
      control = nlme::lmeControl(
        msMaxIter = 100,
        msMaxEval = 200,
        pnlsMaxIter = 25,
        returnObject = TRUE
      )
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
  survival::coxph(
    survival::Surv(time, event) ~ x1 + x2,
    data = sim$dataEvent,
    x = TRUE
  )
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

jm_joint_result <- if (
  is.null(jm_lme_result$error) && is.null(jm_surv_result$error)
) {
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
  list(
    value = NULL,
    error = "Skipped because JMbayes2 preprocessing failed.",
    elapsed = 0
  )
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
  elapsed_seconds = jm_lme_result$elapsed +
    jm_surv_result$elapsed +
    jm_joint_result$elapsed,
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
  longitudinal_records[[length(longitudinal_records) + 1L]] <- extract_jm_long(
    jm_summary,
    truth_long,
    marker_levels,
    replication
  )
  survival_records[[length(survival_records) + 1L]] <- extract_jm_survival(
    jm_summary,
    truth_survival,
    replication
  )
  association_records[[length(association_records) + 1L]] <- extract_jm_assoc(
    jm_summary,
    truth_assoc,
    truth_survival,
    marker_levels,
    replication
  )
}

# Save results to rds
saveRDS(
  list(
    runtime_records = runtime_records,
    longitudinal_records = longitudinal_records,
    survival_records = survival_records,
    association_records = association_records,
    error_records = error_records
  ),
  file.path("records", paste0(base_name, '-rep', replication, '.rds'))
)

# runtime_records <- do.call(rbind, runtime_records)
# longitudinal_records <- bind_rows_or_empty(longitudinal_records)
# survival_records <- bind_rows_or_empty(survival_records)
# association_records <- bind_rows_or_empty(association_records)
# error_records <- bind_rows_or_empty(error_records)

# runtime_summary <- summarise_runtime_metrics(runtime_records, total_reps = cfg$n_rep)
# longitudinal_summary <- summarise_replication_metrics(longitudinal_records, c("marker", "term"), total_reps = cfg$n_rep)
# survival_summary <- summarise_replication_metrics(survival_records, c("term"), total_reps = cfg$n_rep)
# association_summary <- summarise_replication_metrics(association_records, c("marker"), total_reps = cfg$n_rep)

# longitudinal_rhat_comparison <- build_rhat_comparison(longitudinal_summary, c("marker", "term"))
# survival_rhat_comparison <- build_rhat_comparison(survival_summary, c("term"))
# association_rhat_comparison <- build_rhat_comparison(association_summary, c("marker"))

# results <- list(
#   config = cfg,
#   runtime_replications = runtime_records,
#   runtime_summary = runtime_summary,
#   longitudinal_replications = longitudinal_records,
#   longitudinal_summary = longitudinal_summary,
#   survival_replications = survival_records,
#   survival_summary = survival_summary,
#   association_replications = association_records,
#   association_summary = association_summary,
#   rhat_comparison = list(
#     longitudinal = longitudinal_rhat_comparison,
#     survival = survival_rhat_comparison,
#     association = association_rhat_comparison
#   ),
#   errors = error_records
# )

# saveRDS(results, output_rds)

# cat("Benchmark replications complete.\n")
# cat("Results saved to: ", output_rds, "\n", sep = "")
# cat("\nRuntime summary\n")
# print(runtime_summary, row.names = FALSE)
# cat("\nLongitudinal summary (first 12 rows)\n")
# print(utils::head(longitudinal_summary, 12), row.names = FALSE)
# cat("\nSurvival summary\n")
# print(survival_summary, row.names = FALSE)
# cat("\nAssociation summary (first 12 rows)\n")
# print(utils::head(association_summary, 12), row.names = FALSE)
# if (nrow(longitudinal_rhat_comparison)) {
#   cat("\nLongitudinal Rhat comparison (first 12 rows)\n")
#   print(utils::head(longitudinal_rhat_comparison, 12), row.names = FALSE)
# }
# if (nrow(survival_rhat_comparison)) {
#   cat("\nSurvival Rhat comparison\n")
#   print(survival_rhat_comparison, row.names = FALSE)
# }
# if (nrow(association_rhat_comparison)) {
#   cat("\nAssociation Rhat comparison (first 12 rows)\n")
#   print(utils::head(association_rhat_comparison, 12), row.names = FALSE)
# }
# if (nrow(error_records)) {
#   cat("\nErrors (first 12 rows)\n")
#   print(utils::head(error_records, 12), row.names = FALSE)
# }
