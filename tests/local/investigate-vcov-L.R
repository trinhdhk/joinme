devtools::load_all(quiet = TRUE)

formula_long <- y ~ 1 + time + (1 + time | id) + (0 + (1 + time | id) | marker)
formula_event <- survival::Surv(time, event) ~ 1
seed <- 2401
t_ref <- 1.2

softplus <- function(x) {
  ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))
}

summarise_vec <- function(x) {
  q <- stats::quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95), names = FALSE, na.rm = TRUE)
  c(mean = mean(x, na.rm = TRUE), q05 = q[1], q25 = q[2], q50 = q[3], q75 = q[4], q95 = q[5])
}

component_labels <- c("L11", "L21", "L22")

sim <- simulate_joinme(
  formulaLong = formula_long,
  formulaEvent = formula_event,
  formulaVCov = ~ 1,
  families = c("gaussian", "gaussian", "student_t"),
  n_id = 200,
  times_obs = seq(0, 8, length.out = 5),
  assoc = "vcov",
  truth = jm_truth(
    longitudinal = c(-1, 0.5),
    survival = list(slope = 0.5),
    assoc_coef = list(slope = c("vcov[1]" = 2, "vcov[2]" = -1, "vcov[3]" = -1)),
    vcov = list(
      sd = list(intercept = c(-0.5, 1), latent = c(0.1, 0.2)),
      corr = list(intercept = 0.5, latent = -0.15)
    ),
    re_params = list(id = list(sd = c(0.5, 0.25)))
  ),
  transforms = joinme_tf(vcov = ~ log(1 + exp(x))),
  seed = seed,
  use_mirai = TRUE,
  n_workers = 10
)

sd_check <- joinme_standata(
  formulaLong = formula_long,
  dataLong = sim$dataLong,
  formulaEvent = formula_event,
  dataEvent = sim$dataEvent,
  formulaVCov = ~ 1,
  families = c("gaussian", "gaussian", "student_t"),
  assoc = "vcov",
  transforms = joinme_tf(vcov = ~ log(1 + exp(x)))
)

tmax <- sd_check$tmax
n_id <- sd_check$n_id
n_marker <- sd_check$D

truth_L <- sim$truth$L_i
truth_w <- sim$truth$re_draws$id_marker_cov_scaled

truth_raw_vcov <- t(vapply(
  seq_len(n_id),
  function(i) sim$helpers$assoc_components_raw(i, t_ref)$vcov_vals,
  numeric(sd_check$M_vcov_tf)
))
truth_tf_vcov <- t(vapply(
  seq_len(n_id),
  function(i) sim$helpers$assoc_components_raw(i, t_ref)$vcov_vals_tf,
  numeric(sd_check$M_vcov_tf)
))

control <- list(
  engine = "cmdstanr",
  chains = 2,
  parallel_chains = 2,
  iter_warmup = 200,
  iter_sampling = 200,
  threads_per_chain = 6,
  adapt_delta = 0.7,
  max_treedepth = 10,
  force_recompile = FALSE,
  refresh = 100,
  seed = seed
)

fit <- joinme(
  formulaLong = formula_long,
  formulaEvent = formula_event,
  formulaVCov = ~ 1,
  families = c("gaussian", "gaussian", "student_t"),
  dataLong = sim$dataLong,
  dataEvent = sim$dataEvent,
  assoc = "vcov",
  transforms = joinme_tf(vcov = ~ softplus(x)),
  control = control
)

draw_vars <- c(
  paste0("alpha_L[", seq_len(sd_check$M_vcov_tf), "]"),
  paste0("lambda_L[", seq_len(sd_check$M_vcov_tf), "]"),
  as.vector(outer(seq_len(n_id), seq_len(sd_check$M_vcov_tf), function(i, m) paste0("z_L[", i, ",", m, "]"))),
  as.vector(outer(rep(seq_len(n_id), each = n_marker), rep(seq_len(n_marker), times = n_id), function(i, d) paste0("w_idm[", i, ",", d, ",2]")))
)

draws <- extract(
  fit,
  what = "raw",
  variable = draw_vars,
  seed = seed,
  keep_chains = FALSE
)$posterior_draws
draw_names <- colnames(draws)

reconstruct_L_draws <- function(draw_mat, n_subject) {
  n_draw <- nrow(draw_mat)
  out <- array(NA_real_, dim = c(n_draw, n_subject, 3L), dimnames = list(NULL, NULL, component_labels))
  alpha <- draw_mat[, paste0("alpha_L[", 1:3, "]"), drop = FALSE]
  lambda <- draw_mat[, paste0("lambda_L[", 1:3, "]"), drop = FALSE]
  for (i in seq_len(n_subject)) {
    z_mat <- draw_mat[, paste0("z_L[", i, ",", 1:3, "]"), drop = FALSE]
    lp <- alpha + lambda * z_mat
    out[, i, 1] <- softplus(lp[, 1])
    out[, i, 2] <- lp[, 2]
    out[, i, 3] <- softplus(lp[, 3])
  }
  out
}

post_L_draws <- reconstruct_L_draws(draws, n_id)
post_L_mean <- apply(post_L_draws, c(2, 3), mean)
post_L_q10 <- apply(post_L_draws, c(2, 3), stats::quantile, probs = 0.1)
post_L_q90 <- apply(post_L_draws, c(2, 3), stats::quantile, probs = 0.9)

post_tf_vcov_draws <- post_L_draws
post_tf_vcov_draws[, , 1] <- softplus(post_tf_vcov_draws[, , 1]) - softplus(0)
post_tf_vcov_draws[, , 2] <- softplus(post_tf_vcov_draws[, , 2]) - softplus(0)
post_tf_vcov_draws[, , 3] <- softplus(post_tf_vcov_draws[, , 3]) - softplus(0)
post_tf_vcov_mean <- apply(post_tf_vcov_draws, c(2, 3), mean)

truth_L_mat <- cbind(
  L11 = truth_L[, 1, 1],
  L21 = truth_L[, 2, 1],
  L22 = truth_L[, 2, 2]
)

truth_slope_sd_scaled <- sqrt(truth_L_mat[, "L21"]^2 + truth_L_mat[, "L22"]^2)
truth_slope_sd_original <- truth_slope_sd_scaled / tmax
post_slope_sd_scaled_draws <- sqrt(post_L_draws[, , 2]^2 + post_L_draws[, , 3]^2)
post_slope_sd_original_draws <- post_slope_sd_scaled_draws / tmax
post_slope_sd_original_mean <- apply(post_slope_sd_original_draws, 2, mean)

truth_w_slope_original <- as.numeric(truth_w[, , 2]) / tmax
post_w_names <- as.vector(outer(rep(seq_len(n_id), each = n_marker), rep(seq_len(n_marker), times = n_id), function(i, d) paste0("w_idm[", i, ",", d, ",2]")))
post_w_available <- post_w_names[post_w_names %in% draw_names]
post_w_slope_original <- if (length(post_w_available) > 0) {
  as.numeric(draws[, post_w_available, drop = FALSE]) / tmax
} else {
  numeric(0)
}

raw_corr_tbl <- data.frame(
  component = component_labels,
  cor_truth_post_mean = vapply(seq_along(component_labels), function(j) stats::cor(truth_L_mat[, j], post_L_mean[, j]), numeric(1)),
  truth_mean = vapply(seq_along(component_labels), function(j) mean(truth_L_mat[, j]), numeric(1)),
  post_mean = vapply(seq_along(component_labels), function(j) mean(post_L_mean[, j]), numeric(1)),
  truth_sd = vapply(seq_along(component_labels), function(j) stats::sd(truth_L_mat[, j]), numeric(1)),
  post_sd = vapply(seq_along(component_labels), function(j) stats::sd(post_L_mean[, j]), numeric(1)),
  row.names = NULL,
  check.names = FALSE
)

tf_corr_tbl <- data.frame(
  component = paste0("vcov_tf[", seq_len(sd_check$M_vcov_tf), "]"),
  cor_truth_post_mean = vapply(seq_len(sd_check$M_vcov_tf), function(j) stats::cor(truth_tf_vcov[, j], post_tf_vcov_mean[, j]), numeric(1)),
  truth_mean = vapply(seq_len(sd_check$M_vcov_tf), function(j) mean(truth_tf_vcov[, j]), numeric(1)),
  post_mean = vapply(seq_len(sd_check$M_vcov_tf), function(j) mean(post_tf_vcov_mean[, j]), numeric(1)),
  truth_sd = vapply(seq_len(sd_check$M_vcov_tf), function(j) stats::sd(truth_tf_vcov[, j]), numeric(1)),
  post_sd = vapply(seq_len(sd_check$M_vcov_tf), function(j) stats::sd(post_tf_vcov_mean[, j]), numeric(1)),
  row.names = NULL,
  check.names = FALSE
)

subject_tbl <- data.frame(
  id = seq_len(n_id),
  truth_L11 = truth_L_mat[, "L11"],
  post_L11 = post_L_mean[, "L11"],
  truth_L21 = truth_L_mat[, "L21"],
  post_L21 = post_L_mean[, "L21"],
  truth_L22 = truth_L_mat[, "L22"],
  post_L22 = post_L_mean[, "L22"],
  truth_slope_sd_original = truth_slope_sd_original,
  post_slope_sd_original = post_slope_sd_original_mean,
  row.names = NULL,
  check.names = FALSE
)
subject_tbl$abs_err_L22 <- abs(subject_tbl$post_L22 - subject_tbl$truth_L22)
subject_tbl <- subject_tbl[order(subject_tbl$abs_err_L22, decreasing = TRUE), ]

alpha_draws <- draws[, paste0("alpha_L[", 1:3, "]"), drop = FALSE]
lambda_draws <- draws[, paste0("lambda_L[", 1:3, "]"), drop = FALSE]

cat("== Time basis ==\n")
cat("tmax:", format(tmax, digits = 6), "\n")
cat("observed time range:", paste(format(range(sim$dataLong$time), digits = 6), collapse = " to "), "\n")
cat("longitudinal designs use the observed time range directly\n\n")

cat("== Alpha/Lambda posterior ==\n")
print(data.frame(
  parameter = c(paste0("alpha_L[", 1:3, "]"), paste0("lambda_L[", 1:3, "]")),
  truth = c(sim$truth$id_marker_cov_effective$alpha, sim$truth$id_marker_cov_effective$lambda),
  rbind(
    t(apply(alpha_draws, 2, summarise_vec)),
    t(apply(lambda_draws, 2, summarise_vec))
  ),
  row.names = NULL,
  check.names = FALSE
), row.names = FALSE, digits = 4)

cat("\n== Raw L truth vs posterior mean across subjects ==\n")
print(raw_corr_tbl, row.names = FALSE, digits = 4)

cat("\n== Transformed vcov truth vs posterior mean across subjects ==\n")
print(tf_corr_tbl, row.names = FALSE, digits = 4)

cat("\n== Original-time slope SD implied by L ==\n")
print(data.frame(
  quantity = c("truth slope SD", "posterior mean slope SD"),
  rbind(summarise_vec(truth_slope_sd_original), summarise_vec(post_slope_sd_original_mean)),
  row.names = NULL,
  check.names = FALSE
), row.names = FALSE, digits = 4)

cat("\n== Realized marker-by-id slope coefficient on original time scale ==\n")
truth_w_summary <- summarise_vec(truth_w_slope_original)
if (length(post_w_slope_original) > 0) {
  post_w_summary <- summarise_vec(post_w_slope_original)
  print(data.frame(
    quantity = c("truth realized slope coef", "posterior realized slope coef draws"),
    rbind(truth_w_summary, post_w_summary),
    row.names = NULL,
    check.names = FALSE
  ), row.names = FALSE, digits = 4)
} else {
  print(data.frame(
    quantity = "truth realized slope coef",
    t(truth_w_summary),
    row.names = NULL,
    check.names = FALSE
  ), row.names = FALSE, digits = 4)
  cat("posterior w_idm[,,2] draws were not available in the fit output.\n")
}

cat("\n== Subjects with largest posterior L22 discrepancy ==\n")
print(utils::head(subject_tbl[, c(
  "id", "truth_L11", "post_L11", "truth_L21", "post_L21",
  "truth_L22", "post_L22", "truth_slope_sd_original", "post_slope_sd_original", "abs_err_L22"
)], 12), row.names = FALSE, digits = 4)

if (all(c("L_i[1,1,1]", "L_i[1,2,1]", "L_i[1,2,2]") %in% draw_names)) {
  direct_check <- data.frame(
    component = component_labels,
    max_abs_diff = c(
      max(abs(draws[, "L_i[1,1,1]"] - post_L_draws[, 1, 1])),
      max(abs(draws[, "L_i[1,2,1]"] - post_L_draws[, 1, 2])),
      max(abs(draws[, "L_i[1,2,2]"] - post_L_draws[, 1, 3]))
    ),
    row.names = NULL,
    check.names = FALSE
  )
  cat("\n== Direct Stan L_i draw check (subject 1) ==\n")
  print(direct_check, row.names = FALSE, digits = 6)
}

fit_summary <- suppressWarnings(summary(fit, draws = 100, seed = seed))
cat("\n== Association summary ==\n")
print(subset(fit_summary$tables$assoc, grepl("^vcov\\[", term)), row.names = FALSE)
