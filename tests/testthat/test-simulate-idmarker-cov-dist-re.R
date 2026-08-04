test_that("simulate_joinme supports id-specific id:marker covariance via formulaVCov", {
  set.seed(4401)
  sim <- simulate_joinme(
    n_id = 12,
    families = c("gaussian", "gaussian", "gaussian"),
    n_obs_per_marker_per_id = 4,
    times_obs = seq(0, 2, length.out = 5),
    seed = 4401,
    formulaVCov = ~ x1 + x2,
    re_params = list(
      id = list(sd = c(0.5, 0.2)),
      marker = list(sd = 0.3),
      id_marker_cov = list(
        latent = list(sd = c(0.5, 0.25)),
        alpha = c(-0.3, 0.15, -0.1),
        beta = matrix(0.1, nrow = 3, ncol = 2),
        lambda = 0.6,
        diag_link = "softplus"
      )
    )
  )

  raw1 <- sim$helpers$assoc_components_raw(1, 1.0)
  raw2 <- sim$helpers$assoc_components_raw(2, 1.0)

  expect_identical(names(sim$truth$formulaVCov), c("sd", "corr"))
  expect_equal(as.character(sim$truth$formulaVCov$sd), as.character(~ x1 + x2))
  expect_equal(as.character(sim$truth$formulaVCov$corr), as.character(~ x1 + x2))
  expect_true(length(raw1$corr_vals) > 0)
  expect_false(isTRUE(all.equal(raw1$corr_vals, raw2$corr_vals)))
})

test_that("simulate_joinme covariance regression uses component-specific subject latents across L entries", {
  set.seed(4403)
  sim <- simulate_joinme(
    n_id = 40,
    families = c("gaussian", "gaussian", "gaussian"),
    n_obs_per_marker_per_id = 4,
    times_obs = seq(0, 2, length.out = 5),
    seed = 4403,
    formulaVCov = ~ 1,
    re_params = list(
      id_marker_cov = list(
        latent = list(sd = c(0.5, 0.25)),
        alpha = c(-0.3, 0.15, -0.1),
        lambda = c(0.8, -0.4, 0.6),
        diag_link = "softplus"
      )
    )
  )

  Li <- sim$truth$L_i
  lp11 <- log(expm1(Li[, 1, 1]))
  sd2 <- sqrt(
    Li[, 2, 1]^2 + Li[, 2, 2]^2
  ) # reconstructed second-row standard deviation before the Cholesky-correlation factor
  partial_correlation21 <-
    Li[, 2, 1] / sd2 # row-two partial correlation implied by the Cholesky entries
  lp21 <- atanh(pmax(
    -1 + 1e-12,
    pmin(1 - 1e-12, partial_correlation21)
  ))
  lp22 <- log(expm1(sd2))
  corr_vals <- stats::cor(cbind(lp11, lp21, lp22))

  expect_lt(max(abs(corr_vals[upper.tri(corr_vals)])), 0.8)
})

test_that("simulate_joinme canonicalizes covariance-regression loadings to lambda >= 0", {
  input_lambda <- c(0.8, -0.4, 0.6)
  sim <- simulate_joinme(
    n_id = 30,
    families = c("gaussian", "gaussian", "gaussian"),
    n_obs_per_marker_per_id = 4,
    times_obs = seq(0, 2, length.out = 5),
    seed = 4410,
    formulaVCov = ~ 1,
    re_params = list(
      id_marker_cov = list(
        alpha = c(-0.3, 0.15, -0.1),
        lambda = input_lambda,
        diag_link = "softplus"
      )
    )
  )

  eff <- sim$truth$id_marker_cov_effective
  expect_equal(eff$lambda, abs(input_lambda))
  expect_equal(eff$lambda_sign, c(1, -1, 1))
  expect_true(all(eff$lambda >= 0))

  lp11 <- eff$alpha[1] + eff$lambda[1] * eff$z[, 1]
  lp21 <- eff$alpha[2] + eff$lambda[2] * eff$z[, 2]
  lp22 <- eff$alpha[3] + eff$lambda[3] * eff$z[, 3]
  sd2 <- log1p(
    exp(lp22)
  ) # second-coordinate standard deviation after the configured softplus link
  partial_correlation21 <- tanh(
    lp21
  ) # bounded row-two partial correlation after the off-diagonal tanh link

  expect_equal(sim$truth$L_i[, 1, 1], log1p(exp(lp11)), tolerance = 1e-10)
  expect_equal(
    sim$truth$L_i[, 2, 1],
    sd2 * partial_correlation21,
    tolerance = 1e-10
  )
  expect_equal(
    sim$truth$L_i[, 2, 2],
    sd2 * sqrt(1 - partial_correlation21^2),
    tolerance = 1e-10
  )
})

test_that("simulate_joinme default covariance regression uses component-specific loadings", {
  set.seed(4404)
  sim <- simulate_joinme(
    n_id = 12,
    families = c("gaussian", "gaussian", "gaussian"),
    n_obs_per_marker_per_id = 4,
    times_obs = seq(0, 2, length.out = 5),
    seed = 4404,
    assoc = "vcov"
  )

  eff <- sim$truth$id_marker_cov_effective

  expect_length(eff$alpha, 3)
  expect_length(eff$lambda, 3)
  expect_false(length(unique(round(eff$alpha, 8))) == 1L)
  expect_false(length(unique(round(eff$lambda, 8))) == 1L)
})

test_that("covariance-regression lambda is scalar or coordinate vector, not a matrix", {
  expect_error(
    simulate_joinme(
      formulaLong = y ~ 1 + time +
        (1 + time | id) +
        (0 + (1 + time | id) | marker),
      formulaEvent = survival::Surv(time, event) ~ 1,
      n_id = 4,
      families = c("gaussian", "gaussian"),
      times_obs = c(0, 0.2),
      time_cens = 0.3,
      re_params = list(
        id_marker_cov = list(
          lambda = matrix(c(0.5, 0, 0, 0.5), nrow = 2)
        )
      ),
      use_mirai = FALSE,
      seed = 4411
    ),
    "scalar or numeric vector"
  )
})

test_that("simulate_joinme default covariance regression initialises without length recycling", {
  set.seed(4405)
  sim <- simulate_joinme(
    formulaLong = y ~ 1 + time + (1 + time | id) + (0 + (1 + time | id) | marker),
    formulaEvent = survival::Surv(time, event) ~ 1,
    families = c("gaussian", "gaussian", "gaussian"),
    n_id = 8,
    n_obs_per_marker_per_id = 3,
    times_obs = seq(0, 2, length.out = 4),
    assoc = "vcov",
    seed = 4405
  )

  expect_equal(dim(sim$truth$L_i)[2:3], c(2, 2))
  expect_true(all(is.finite(sim$truth$L_i)))
})

test_that("simulate_joinme randomises omitted coefficients reproducibly and exposes effective values", {
  sim1 <- simulate_joinme(
    n_id = 10,
    families = c("gaussian", "gaussian", "gaussian"),
    n_obs_per_marker_per_id = 3,
    times_obs = seq(0, 2, length.out = 4),
    assoc = c("cv_mean", "vcov"),
    seed = 4406
  )

  sim2 <- simulate_joinme(
    n_id = 10,
    families = c("gaussian", "gaussian", "gaussian"),
    n_obs_per_marker_per_id = 3,
    times_obs = seq(0, 2, length.out = 4),
    assoc = c("cv_mean", "vcov"),
    seed = 4406
  )

  eff <- sim1$truth$re_structure

  expect_equal(sim1$truth$beta_long, sim2$truth$beta_long)
  expect_equal(sim1$truth$beta_event, sim2$truth$beta_event)
  expect_equal(sim1$truth$assoc_coefs, sim2$truth$assoc_coefs)
  expect_equal(eff$id$sd, sim2$truth$re_structure$id$sd)
  expect_equal(eff$id$corr, sim2$truth$re_structure$id$corr)
  expect_equal(eff$id$Lcorr, sim2$truth$re_structure$id$Lcorr)
  expect_equal(eff$id$cov, sim2$truth$re_structure$id$cov)
  expect_equal(eff$id_marker_cov$latent$sd, sim2$truth$re_structure$id_marker_cov$latent$sd)
  expect_equal(eff$id_marker_cov$alpha, sim2$truth$re_structure$id_marker_cov$alpha)
  expect_true(all(eff$id$sd > 0))
  expect_true(all(eff$marker$sd > 0))
  expect_identical(eff$id_marker_cov$latent$mode, "iid_standard_normal")
  expect_true(all(eff$id_marker_cov$latent$sd > 0))
  expect_equal(eff$id_marker_cov$latent$corr, diag(length(eff$id_marker_cov$latent$sd)))
  expect_equal(diag(eff$id$corr), rep(1, ncol(eff$id$corr)))
  expect_true(isTRUE(all.equal(eff$id$corr, t(eff$id$corr), tolerance = 1e-8)))
  expect_true(all(eff$id_marker_cov$lambda >= 0))
  expect_true(all(is.finite(sim1$truth$assoc_coefs)))
  expect_true(any(grepl("^vcov\\[", names(sim1$truth$assoc_coefs))))
})

test_that("simulate_joinme translates legacy latent covariance input into baseline alpha", {
  sim <- simulate_joinme(
    n_id = 8,
    families = c("gaussian", "gaussian", "gaussian"),
    n_obs_per_marker_per_id = 3,
    times_obs = seq(0, 2, length.out = 4),
    seed = 4408,
    assoc = "vcov",
    re_params = list(
      id_marker_cov = list(
        latent = list(sd = c(0.8, 0.4), corr = matrix(c(1, 0.2, 0.2, 1), 2, 2)),
        alpha = c(-0.2, 0.05, -0.1),
        lambda = 0
      )
    )
  )

  eff <- sim$truth$id_marker_cov_effective

  expect_identical(eff$latent$mode, "iid_standard_normal")
  expect_false(is.null(eff$legacy_latent_translation$input))
  expect_identical(eff$legacy_latent_translation$applied_to, "alpha")
  expect_true(all(dim(eff$legacy_latent_translation$factor) == c(2, 2)))
  expect_false(isTRUE(all.equal(eff$alpha, c(-0.2, 0.05, -0.1))))
})

test_that("simulate_joinme rejects invalid supplied random-effect correlation matrices", {
  expect_error(
    simulate_joinme(
      n_id = 6,
      families = c("gaussian", "gaussian"),
      n_obs_per_marker_per_id = 3,
      times_obs = seq(0, 2, length.out = 4),
      seed = 4407,
      re_params = list(
        id = list(
          sd = c(0.4, 0.2),
          corr = matrix(c(1, 1.2, 0, 1), nrow = 2, byrow = TRUE)
        )
      )
    ),
    "must be symmetric|must lie in \\[-1, 1\\]|must be positive definite"
  )
})

test_that("simulate_joinme supports distributional random effects", {
  set.seed(4402)
  sim <- simulate_joinme(
    n_id = 10,
    families = c("student_t", "student_t"),
    n_obs_per_marker_per_id = 5,
    times_obs = seq(0, 2, length.out = 6),
    seed = 4402,
    formulaDist = list(sigma ~ 1 + (1 | id)),
    re_params = list(
      id = list(sd = c(0.4, 0.2)),
      marker = list(sd = 0.25),
      id_marker_cov = list(latent = list(sd = c(0.4, 0.2))),
      dist = list(sigma = list(sd = 0.8))
    )
  )

  expect_true(is.list(sim$truth$dist_re_params))
  expect_true("sigma" %in% names(sim$truth$dist_re_params))
  expect_equal(sim$truth$dist_re_params$sigma$sd, 0.8)
  expect_true(all(is.finite(sim$dataLong$y)))
})

test_that("simulate_joinme saves realized distributional parameter truth and avoids duplicate beta aliases", {
  sim <- simulate_joinme(
    n_id = 8,
    families = c("student_t", "skew_normal", "beta"),
    n_obs_per_marker_per_id = 4,
    times_obs = seq(0, 2, length.out = 5),
    seed = 4411,
    formulaDist = list(
      sigma ~ 1 + time,
      nu[family = student_t] ~ 1,
      alpha[family = skew_normal] ~ 1,
      kappa[family = beta] ~ 1
    )
  )

  dp <- sim$truth$distributional_params

  expect_true(is.list(dp))
  expect_equal(nrow(dp$rowwise), nrow(sim$dataLong))
  expect_equal(length(dp$sigma), nrow(sim$dataLong))
  expect_equal(length(dp$nu), nrow(sim$dataLong))
  expect_equal(length(dp$alpha_skew), nrow(sim$dataLong))
  expect_equal(length(dp$kappa), nrow(sim$dataLong))
  expect_equal(as.character(dp$rowwise$marker), as.character(sim$dataLong$marker))
  expect_equal(as.numeric(dp$rowwise$time), as.numeric(sim$dataLong$time))
  expect_true(all(is.finite(dp$sigma)))
  expect_true(all(is.finite(dp$nu)))
  expect_true(all(is.finite(dp$alpha_skew)))
  expect_true(all(is.finite(dp$kappa)))
  expect_false("beta" %in% names(sim$truth))
  expect_false("gamma_w" %in% names(sim$truth))
})

test_that("simulate_joinme exposes fit-aligned truth for defaults and scaled parameters", {
  sim <- simulate_joinme(
    n_id = 8,
    families = c("gaussian", "student_t", "beta", "cumulative_logit"),
    n_obs_per_marker_per_id = 4,
    times_obs = seq(0, 2, length.out = 5),
    seed = 4412
  )

  fit_truth <- sim$truth$stan_fit
  fam_truth <- sim$truth$family

  expect_true(is.list(fit_truth))
  expect_equal(unname(fit_truth$tau_u), sim$truth$re_structure$id$sd)
  expect_equal(unname(fit_truth$Lcorr_u), sim$truth$re_structure$id$Lcorr)
  expect_equal(unname(fit_truth$Sigma_u), sim$truth$re_structure$id$cov)
  expect_equal(length(fit_truth$beta_scaled), length(sim$truth$beta_long))
  expect_true(is.finite(fit_truth$time_scale_generation))
  expect_true(is.finite(fit_truth$time_scale_observed_max))
  expect_gt(fit_truth$time_scale_generation, 0)
  expect_gt(fit_truth$time_scale_observed_max, 0)
  expect_equal(length(fit_truth$marker_id_row_scale_eff), length(sim$truth$re_structure$id_marker_cov$latent$sd))
  expect_true("student_t" %in% names(fit_truth$sigma_family))
  expect_true("student_t" %in% names(fit_truth$nu_family))
  expect_true("beta" %in% names(fit_truth$kappa_family))
  expect_equal(unname(fit_truth$cutpoints_ord), c(-1, 1))
  expect_equal(nrow(fam_truth$by_marker), 4)
  expect_equal(fam_truth$by_marker$marker_to_sigma_family, c(1L, 2L, 0L, 0L))
  expect_equal(fam_truth$by_marker$marker_to_nu_family, c(0L, 1L, 0L, 0L))
  expect_equal(fam_truth$by_marker$marker_to_kappa_family, c(0L, 0L, 1L, 0L))
})

test_that("simulate_joinme exposes public random-effect draws on the original-time scale", {
  sim <- simulate_joinme(
    formulaLong = y ~ 1 + time + (1 + time | id) + (1 + time + (1 + time | id) | marker),
    n_id = 10,
    families = rep("gaussian", 2),
    n_obs_per_marker_per_id = 4,
    times_obs = seq(0, 4, length.out = 5),
    seed = 4421
  )

  scale_factor <- unname(sim$truth$stan_fit$tau_v_eff[2] / sim$truth$re_structure$marker$sd[2])
  expect_true(is.finite(scale_factor))
  expect_gt(scale_factor, 0)

  expect_equal(sim$truth$re_draws$marker[, 1], sim$truth$re_draws_likelihood$marker[, 1])
  expect_equal(
    sim$truth$re_draws$marker[, 2] * scale_factor,
    sim$truth$re_draws_likelihood$marker[, 2],
    tolerance = 1e-8
  )
  expect_equal(
    sim$truth$re_draws$id[, 2] * scale_factor,
    sim$truth$re_draws_likelihood$id[, 2],
    tolerance = 1e-8
  )
})

test_that("joinme_standata time scale is determined from observed event times", {
  sim <- simulate_joinme(
    n_id = 10,
    families = rep("gaussian", 2),
    n_obs_per_marker_per_id = 3,
    times_obs = seq(0, 8, length.out = 10),
    time_cens = 8,
    seed = 4424
  )

  sd <- joinme_standata(
    formulaLong = sim$truth$formulaLong,
    dataLong = sim$dataLong,
    formulaEvent = sim$truth$formulaEvent,
    dataEvent = sim$dataEvent,
    assoc = sim$truth$assoc,
    families = sim$marker_info$families
  )

  expect_equal(as.numeric(sd$tmax), max(sim$dataEvent$time_stop, na.rm = TRUE))

  idx_time_beta <- as.integer(sim$truth$stan_fit$idx_time_beta %||% integer(0))
  if (length(idx_time_beta) > 0L) {
    expect_equal(
      sim$truth$stan_fit$beta_scaled[idx_time_beta],
      sim$truth$beta_long[idx_time_beta] * as.numeric(sd$tmax),
      tolerance = 1e-8
    )
  }
})

test_that("simulate_joinme draws id random effects from zero-mean MVN", {
  sim <- simulate_joinme(
    formulaLong = y ~ 1 + time + (1 + time | id) + (0 + (1 + time | id) | marker),
    n_id = 600,
    families = rep("gaussian", 2),
    n_obs_per_marker_per_id = 2,
    times_obs = seq(0, 2, length.out = 3),
    seed = 4423,
    re_params = list(
      id = list(sd = c(0.6, 0.3), corr = matrix(c(1, 0.25, 0.25, 1), 2, 2))
    )
  )

  u_draw <- sim$truth$re_draws_likelihood$id
  sigma_u <- sim$truth$stan_fit$Sigma_u

  expect_true(is.matrix(u_draw))
  expect_equal(colMeans(u_draw), c(0, 0), tolerance = 0.07)
  expect_equal(stats::cov(u_draw), sigma_u, tolerance = 0.08)
})
