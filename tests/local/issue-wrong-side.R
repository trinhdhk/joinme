devtools::load_all(quiet = TRUE)
library(dplyr)

formula_long <- y ~ 1 + time + (1 + time | id) + (1 + time + (1 + time | id) | marker)
formula_event <- survival::Surv(time, event) ~ 1
seed <- 2401

sim <- simulate_joinme(
  formulaLong = formula_long,
  formulaEvent = formula_event,
  formulaVCov = ~ 1,
  families = rep('gaussian', 20),
  # n_markers = 6,
  n_id = 80,
  times_obs = seq(0, 10, length.out = 8),
  obs_time_noise_sd = 0.5,
  time_cens = 10,
  assoc = "vcov",
  truth = jm_truth(
    longitudinal = c(-0.5, 0.5),
    survival = list(slope = 0.5),
    assoc_coef = list(slope = c("vcov[1]" = 2, "vcov[2]" = -1, "vcov[3]" = 0)),
    vcov = list(
      sd = list(intercept = c(0.5, -0.2), latent = c(1, 0.25)),
      corr = list(intercept = -0.2, latent = 0.5)
    ),
    basehaz = list(type = "weibull", shape = 0.8, scale = 10),
    re_params = list(
      id = list(sd = c(1, 0.5), corr = matrix(c(1, -0.5, -0.5, 1), ncol = 2))
    )
  ),
  transforms = joinme_tf(vcov = ~ softplus(x)),
  seed = seed,
  use_mirai = TRUE,
  n_workers = 10
)

sim$dataLong |>
  tidytable::filter(marker == "m1") |>
  tidytable::summarize(n = tidytable::n(), .by = id)
sim$dataEvent$time
sd_check <- joinme_standata(
  formulaLong = formula_long,
  dataLong = sim$dataLong,
  formulaEvent = formula_event,
  dataEvent = sim$dataEvent,
  formulaVCov = ~ 1,
  families = 'gaussian',
  assoc = "vcov",
  transforms = joinme_tf(vcov = ~ log(1 + exp(x)))
)

k_cov <- ncol(stats::model.matrix(sim$truth$formulaVCov, sim$dataEvent))
if ("(Intercept)" %in% colnames(stats::model.matrix(sim$truth$formulaVCov, sim$dataEvent))) {
  k_cov <- k_cov - 1L
}

cat("== Simulation truth ==\n")
cat("Q_idm:", sd_check$Q_idm, "\n")
cat("M_vcov_tf:", sd_check$M_vcov_tf, "\n")
cat("zidm_cols:", paste(sd_check$zidm_cols, collapse = ", "), "\n")
cat("K_cov:", k_cov, "\n")
cat("formulaVCov:", paste(deparse(sim$truth$formulaVCov), collapse = " "), "\n")
cat("requested alpha:", -0.5, "\n")
cat("alpha truth:", paste(format(sim$truth$id_marker_cov_effective$alpha, digits = 6), collapse = ", "), "\n")
cat("lambda truth:", paste(format(sim$truth$id_marker_cov_effective$lambda, digits = 6), collapse = ", "), "\n")
cat("assoc truth:", paste(names(sim$truth$assoc_coefs), format(sim$truth$assoc_coefs, digits = 6), collapse = ", "), "\n")
cat("latent z scale truth: 1\n\n")

vcov_tf_by_id <- t(vapply(
  seq_len(nrow(sim$dataEvent)),
  function(i) sim$helpers$assoc_components_raw(i, 1.2)$vcov_vals_tf,
  numeric(sd_check$M_vcov_tf)
))

cat("== Simulated vcov transform correlation ==\n")
print(round(stats::cor(vcov_tf_by_id), 4))
cat("vcov tf sd:", paste(format(apply(vcov_tf_by_id, 2, stats::sd), digits = 4), collapse = ", "), "\n\n")

control <- list(
  engine = "cmdstanr",
  chains = 2,
  parallel_chains = 2,
  iter_warmup = 800,
  iter_sampling = 800,
  threads_per_chain = 6,
  adapt_delta = 0.68,
  max_treedepth = 13,
  init = 1,
  force_recompile = FALSE,
  refresh = 100,
  seed = seed
)

fit <- joinme(
  formulaLong = formula_long,
  formulaEvent = formula_event,
  formulaVCov = ~ 1,
  families = 'gaussian',
  dataLong = sim$dataLong,
  dataEvent = sim$dataEvent,
  assoc = "vcov",
  transforms = joinme_tf(vcov = ~ softplus(x)),
  control = control
)

draw_vars <- c(
  paste0("alpha_L[", seq_len(sd_check$M_vcov_tf), "]"),
  paste0("lambda_L[", seq_len(sd_check$M_vcov_tf), "]"),
  paste0("alpha_vcov[", seq_len(sd_check$M_vcov_tf), "]"),
  paste0("z_alpha_vcov[", seq_len(sd_check$M_vcov_tf), "]"),
  "s_vcov"
)

draws <- extract(
  fit,
  what = "raw",
  variable = draw_vars,
  seed = seed,
  keep_chains = FALSE
)$posterior_draws

draws_df <- posterior_draws(fit, format = "draws_df")
tracked_vars <- c(
  paste0("alpha_L[", seq_len(sd_check$M_vcov_tf), "]"),
  paste0("lambda_L[", seq_len(sd_check$M_vcov_tf), "]"),
  paste0("alpha_vcov[", seq_len(sd_check$M_vcov_tf), "]"),
  paste0("z_alpha_vcov[", seq_len(sd_check$M_vcov_tf), "]"),
  "s_vcov"
)

chain_tbl <- data.frame(
  chain = draws_df$.chain,
  draws_df[, tracked_vars, drop = FALSE],
  check.names = FALSE,
  stringsAsFactors = FALSE
)

path_tbl <- data.frame(
  draw = seq_len(nrow(draws)),
  draws[, tracked_vars, drop = FALSE],
  check.names = FALSE,
  stringsAsFactors = FALSE
)

summarise_path <- function(x) {
  q <- stats::quantile(x, c(0.05, 0.25, 0.75, 0.95), names = FALSE)
  c(
    mean = mean(x),
    median = stats::median(x),
    q05 = q[1],
    q25 = q[2],
    q75 = q[3],
    q95 = q[4]
  )
}

truth_map <- c(
  stats::setNames(sim$truth$id_marker_cov_effective$alpha, paste0("alpha_L[", seq_len(sd_check$M_vcov_tf), "]")),
  stats::setNames(sim$truth$id_marker_cov_effective$lambda, paste0("lambda_L[", seq_len(sd_check$M_vcov_tf), "]")),
  sim$truth$assoc_coefs[paste0("vcov[", seq_len(sd_check$M_vcov_tf), "]")]
)
names(truth_map)[grepl("^vcov\\[", names(truth_map))] <- paste0("alpha_vcov[", seq_len(sd_check$M_vcov_tf), "]")
truth_map[["s_vcov"]] <- NA_real_

raw_assoc_truth <- data.frame(
  parameter = paste0("z_alpha_vcov[", seq_len(sd_check$M_vcov_tf), "]"),
  note = "standardized latent; compare alpha_vcov to hazard-scale truth",
  row.names = NULL,
  check.names = FALSE
)

recovery_tbl <- do.call(rbind, lapply(tracked_vars, function(var_name) {
  data.frame(
    parameter = var_name,
    truth = unname(truth_map[[var_name]] %||% NA_real_),
    t(summarise_path(path_tbl[[var_name]])),
    row.names = NULL,
    check.names = FALSE
  )
}))
recovery_tbl <- tidytable::left_join(recovery_tbl, raw_assoc_truth, by = "parameter")

prior_scale <- 1
set.seed(seed)
n_prior <- 200000
prior_lambda <- prior_scale * stats::rt(n_prior, df = 6)
prior_ref <- data.frame(
  quantity = "lambda_L",
  q02.5 = stats::quantile(prior_lambda, 0.025),
  median = stats::quantile(prior_lambda, 0.5),
  q97.5 = stats::quantile(prior_lambda, 0.975),
  row.names = NULL,
  check.names = FALSE
)

chain_summary <- stats::aggregate(chain_tbl[, tracked_vars, drop = FALSE], by = list(chain = chain_tbl$chain), FUN = mean)

cat("== Recovery summary ==\n")
print(recovery_tbl, row.names = FALSE, digits = 4)
cat("\nNote: alpha_vcov[...] is the hazard-scale association; z_alpha_vcov[...] is the standardized latent before multiplying by s_vcov.\n")
cat("\n== Chain means ==\n")
print(chain_summary, row.names = FALSE, digits = 4)
cat("\n== Prior reference (current Stan priors) ==\n")
print(prior_ref, row.names = FALSE, digits = 4)
cat("\n== Early draw path ==\n")
print(utils::head(path_tbl[, c("draw", tracked_vars), drop = FALSE], 12), row.names = FALSE, digits = 4)
cat("\n== Late draw path ==\n")
print(utils::tail(path_tbl[, c("draw", tracked_vars), drop = FALSE], 12), row.names = FALSE, digits = 4)

fit_summary <- suppressWarnings(summary(fit, seed = seed))

cat("\n== Summary table: association ==\n")
print(subset(fit_summary$tables$assoc, grepl("^vcov\\[", term)), row.names = FALSE)

cat("\n== Summary table: covariance regression ==\n")
print(fit_summary$tables$id_marker_cov$regression, row.names = FALSE)
