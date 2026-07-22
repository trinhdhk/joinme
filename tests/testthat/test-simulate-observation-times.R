test_that("simulate_joinme clips observation-time jitter to the study window", {
  sim <- simulate_joinme(
    n_id = 3,
    families = rep("gaussian", 2),
    times_obs = seq(0, 4, length.out = 5),
    obs_time_noise_sd = 0.2,
    time_cens = 4,
    censor_longitudinal_after_event = FALSE,
    seed = 4201,
    use_mirai = FALSE,
    assoc = c("cv_total")
  )

  time_template <- seq(0, 4, length.out = 5)
  by_block <- split(sim$dataLong$time, interaction(sim$dataLong$id, sim$dataLong$marker, drop = TRUE))
  deviations <- unlist(lapply(by_block, function(tt) sort(tt) - time_template), use.names = FALSE)

  expect_true(all(sim$dataLong$time >= -1e-8))
  expect_true(all(sim$dataLong$time <= 4 + 1e-8))
  expect_true(any(abs(deviations) > 1e-8))
})

test_that("simulate_joinme evaluates measurements at the jittered observation times", {
  sim <- simulate_joinme(
    formulaLong = y ~ 1 + time + (1 + time | id),
    formulaEvent = survival::Surv(time, event) ~ 1,
    n_id = 2,
    families = "gaussian",
    times_obs = seq(0, 3, length.out = 4),
    obs_time_noise_sd = 0.25,
    censor_longitudinal_after_event = FALSE,
    beta_long = c(0, 1),
    re_params = list(id = list(sd = c(1e-12, 1e-12), corr = diag(2))),
    family_params = list(gaussian = list(sigma = 0)),
    seed = 4204,
    use_mirai = FALSE,
    assoc = c("cv_total")
  )

  expect_equal(sim$dataLong$y, sim$dataLong$time, tolerance = 1e-8)
})

test_that("simulate_joinme honors raw-time fixed-effect slopes", {
  sim <- simulate_joinme(
    formulaLong = y ~ 1 + time + (1 + time | id),
    formulaEvent = survival::Surv(time, event) ~ 1,
    n_id = 200,
    families = "gaussian",
    times_obs = seq(0, 8, length.out = 9),
    obs_time_noise_sd = 0,
    censor_longitudinal_after_event = FALSE,
    beta_long = c(0, 0.5),
    re_params = list(id = list(sd = c(1e-8, 1e-8), corr = diag(2))),
    family_params = list(gaussian = list(sigma = 0)),
    seed = 4205,
    use_mirai = FALSE,
    assoc = c("cv_total")
  )

  coef_hat <- stats::coef(stats::lm(y ~ time, data = sim$dataLong))
  expect_equal(unname(coef_hat[["time"]]), 0.5, tolerance = 1e-8)
})

test_that("simulate_joinme keeps vcov association features on the user scale", {
  sim <- simulate_joinme(
    formulaLong = y ~ 1 + time + (1 + time | id) + (0 + (1 + time | id) | marker),
    formulaEvent = survival::Surv(time, event) ~ x1,
    formulaVCov = ~ 1,
    families = rep("gaussian", 8),
    n_id = 80,
    times_obs = seq(0, 10, length.out = 10),
    obs_time_noise_sd = 0.25,
    time_cens = 10,
    assoc = "vcov",
    assoc_coefs = list(vcov = c(1, -1, -1)),
    transforms = joinme_tf(vcov = ~ softplus(x)),
    baseline_hazard = list(type = "weibull", shape = 1, scale = 6),
    beta_long = c(-1, 0.5),
    beta_event = c(x1 = 0.5),
    re_params = list(
      id = list(sd = c(1, 0.25), corr = matrix(c(1, -0.5, -0.5, 1), ncol = 2)),
      id_marker_cov = list(alpha = c(0.5, -0.5, -0.5), lambda = c(1, 0.15, 0.2))
    ),
    seed = 4206,
    use_mirai = FALSE
  )

  vcov_vals <- t(vapply(
    seq_len(nrow(sim$dataEvent)),
    function(i) sim$helpers$assoc_components_raw(i, 1)$vcov_vals,
    numeric(3)
  ))

  expect_true(mean(vcov_vals[, 3]) < 1)
  expect_true(mean(sim$dataEvent$event) > 0.1)
})

test_that("simulate_joinme ignores legacy count args silently", {
  expect_no_warning({
    sim <- simulate_joinme(
      n_id = 2,
      families = rep("gaussian", 2),
      times_obs = seq(0, 3, length.out = 4),
      n_obs_per_marker_per_id = 99,
      n_t = 123,
      seed = 4203,
      use_mirai = FALSE,
      assoc = c("cv_total")
    )
  })

  expect_equal(nrow(sim$dataLong), 2 * 2 * 4)
})

test_that("simulate_joinme can retain longitudinal measurements after the event time", {
  sim_censored <- simulate_joinme(
    n_id = 8,
    families = rep("gaussian", 2),
    times_obs = seq(0, 5, length.out = 6),
    seed = 4202,
    use_mirai = FALSE,
    assoc = c("cv_total"),
    formulaEvent = survival::Surv(time, event) ~ 1,
    baseline_hazard = list(type = "piecewise", breaks = c(1, 2), rates = c(2, 2, 2)),
    censor_longitudinal_after_event = TRUE
  )

  sim_uncensored <- simulate_joinme(
    n_id = 8,
    families = rep("gaussian", 2),
    times_obs = seq(0, 5, length.out = 6),
    seed = 4202,
    use_mirai = FALSE,
    assoc = c("cv_total"),
    formulaEvent = survival::Surv(time, event) ~ 1,
    baseline_hazard = list(type = "piecewise", breaks = c(1, 2), rates = c(2, 2, 2)),
    censor_longitudinal_after_event = FALSE
  )

  expect_equal(sim_censored$dataEvent$time, sim_uncensored$dataEvent$time, tolerance = 1e-10)
  expect_equal(sim_censored$dataEvent$event, sim_uncensored$dataEvent$event)

  event_ids <- sim_censored$dataEvent$id[sim_censored$dataEvent$event == 1L & sim_censored$dataEvent$time < 5]
  expect_true(length(event_ids) > 0)

  max_time_by_id_censored <- tapply(sim_censored$dataLong$time, sim_censored$dataLong$id, max)
  max_time_by_id_uncensored <- tapply(sim_uncensored$dataLong$time, sim_uncensored$dataLong$id, max)
  n_obs_by_id_censored <- tapply(sim_censored$dataLong$time, sim_censored$dataLong$id, length)
  n_obs_by_id_uncensored <- tapply(sim_uncensored$dataLong$time, sim_uncensored$dataLong$id, length)
  event_time_by_id <- setNames(sim_censored$dataEvent$time, sim_censored$dataEvent$id)

  expect_true(all(max_time_by_id_censored[as.character(event_ids)] <= event_time_by_id[as.character(event_ids)] + 1e-8))
  expect_true(any(max_time_by_id_uncensored[as.character(event_ids)] > event_time_by_id[as.character(event_ids)] + 1e-8))
  expect_true(any(n_obs_by_id_censored[as.character(event_ids)] < length(seq(0, 5, length.out = 6)) * 2))
  expect_true(all(n_obs_by_id_uncensored == length(seq(0, 5, length.out = 6)) * 2))
  expect_equal(sort(unique(sim_censored$dataLong$id)), sort(sim_censored$dataEvent$id))
})
