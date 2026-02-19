test_that("normalize_formula_dist parses family-scoped syntax", {
  spec <- joinme:::.normalize_formula_dist(list(
    sigma[family = student] ~ 1 + time,
    sigma[family = normal] ~ 1 + x1,
    alpha_skew[family = skew_normal] ~ 1,
    phi_beta ~ 1
  ))

  expect_true(joinme:::.is_dist_scope(spec$sigma))
  expect_setequal(names(spec$sigma$by_family), c("student_t", "gaussian"))
  expect_true(inherits(spec$sigma$by_family$student_t, "formula"))
  expect_true(inherits(spec$alpha$by_family$skew_normal, "formula"))
  expect_true(inherits(spec$phi_beta$default, "formula"))
})


test_that("joinme_standata builds family-gated distributional matrices", {
  set.seed(1901)
  sim <- simulate_joinme(
    n_id = 6,
    families = c("gaussian", "student_t"),
    n_obs_per_marker_per_id = 3,
    times_obs = seq(0, 2, length.out = 5),
    seed = 1901,
    assoc = c("cv_total"),
    assoc_coefs = c(cv_total = 0.2)
  )

  formulaLong <- y ~ 1 + time + x1 +
    (1 + time | id) +
    (0 + x1 + (1 + time | id) | marker)
  formulaEvent <- survival::Surv(time, event) ~ 1 + x1 + x2

  sd <- joinme_standata(
    formulaLong = formulaLong,
    dataLong = sim$dataLong,
    formulaEvent = formulaEvent,
    dataEvent = sim$dataEvent,
    families = c("gaussian", "student_t"),
    formulaDist = list(
      sigma[family = gaussian] ~ 1 + x1,
      sigma[family = student_t] ~ 1 + time,
      nu[family = student_t] ~ 1
    ),
    assoc = c("cv_total")
  )

  expect_equal(sd$P_sigma, 4L)
  expect_true(any(grepl("^family=gaussian::", sd$dist_cols$sigma)))
  expect_true(any(grepl("^family=student_t::", sd$dist_cols$sigma)))
  expect_true(any(grepl("^family=student_t::", sd$dist_cols$nu)))
  expect_equal(sd$P_nu, 1L)
})


test_that("joinme_standata rejects unsupported parameter-family scopes", {
  set.seed(1902)
  sim <- simulate_joinme(
    n_id = 5,
    families = c("gaussian", "student_t"),
    n_obs_per_marker_per_id = 3,
    times_obs = seq(0, 2, length.out = 5),
    seed = 1902,
    assoc = c("cv_total"),
    assoc_coefs = c(cv_total = 0.2)
  )

  formulaLong <- y ~ 1 + time + x1 +
    (1 + time | id) +
    (0 + x1 + (1 + time | id) | marker)
  formulaEvent <- survival::Surv(time, event) ~ 1 + x1 + x2

  expect_error(
    joinme_standata(
      formulaLong = formulaLong,
      dataLong = sim$dataLong,
      formulaEvent = formulaEvent,
      dataEvent = sim$dataEvent,
      families = c("gaussian", "student_t"),
      formulaDist = list(nu[family = gaussian] ~ 1),
      assoc = c("cv_total")
    ),
    "not used by family"
  )
})


test_that("simulate_joinme supports family-scoped formulaDist", {
  set.seed(1903)
  sim <- simulate_joinme(
    n_id = 8,
    families = c("gaussian", "student_t", "skew_normal"),
    n_obs_per_marker_per_id = 3,
    times_obs = seq(0, 2, length.out = 5),
    seed = 1903,
    formulaDist = list(
      sigma[family = gaussian] ~ 1 + x1,
      sigma[family = student_t] ~ 1 + time,
      alpha_skew[family = skew_normal] ~ 1
    )
  )

  expect_true(is.data.frame(sim$dataLong))
  expect_true(is.data.frame(sim$dataEvent))
  expect_gt(nrow(sim$dataLong), 0)
})


test_that("joinme fits and predicts with family-scoped formulaDist", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  has_cmdstan <- FALSE
  tryCatch({
    has_cmdstan <- !is.null(cmdstanr::cmdstan_version())
  }, error = function(e) {
    has_cmdstan <- FALSE
  })
  if (!has_cmdstan) skip("CmdStan is not installed.")

  set.seed(1904)
  sim <- simulate_joinme(
    n_id = 5,
    families = c("gaussian", "student_t"),
    n_obs_per_marker_per_id = 3,
    times_obs = seq(0, 2, length.out = 5),
    seed = 1904,
    assoc = c("cv_total"),
    assoc_coefs = c(cv_total = 0.2)
  )

  formulaLong <- y ~ 1 + time + x1 +
    (1 + time | id) +
    (0 + x1 + (1 + time | id) | marker)
  formulaEvent <- survival::Surv(time, event) ~ 1 + x1 + x2

  fit <- joinme(
    formulaLong = formulaLong,
    dataLong = sim$dataLong,
    formulaEvent = formulaEvent,
    dataEvent = sim$dataEvent,
    families = c("gaussian", "student_t"),
    formulaDist = list(
      sigma[family = gaussian] ~ 1 + x1,
      sigma[family = student_t] ~ 1 + time,
      nu[family = student_t] ~ 1
    ),
    assoc = c("cv_total"),
    transforms = list(cv_total = list(type = "identity")),
    control = list(
      engine = "cmdstanr",
      chains = 1,
      parallel_chains = 1,
      iter_warmup = 60,
      iter_sampling = 60,
      refresh = 0,
      seed = 1904
    )
  )

  expect_s3_class(fit, "JoinMeFit")
  expect_true(fit$stan_data$P_sigma >= 2)
  expect_true(any(grepl("^family=student_t::", fit$stan_data$dist_cols$sigma)))

  ndL <- sim$dataLong[sim$dataLong$id == 1, ]
  ndE <- sim$dataEvent[sim$dataEvent$id == 1, ]

  pred <- posterior_epred(
    fit,
    newdataLong = ndL,
    newdataEvent = ndE,
    time_start = max(ndL$time),
    n_samples = 10,
    n_times = 50,
    control = list(
      engine = "cmdstanr",
      chains = 1,
      iter_warmup = 20,
      iter_sampling = 1,
      refresh = 0
    ),
    seed = 1905
  )

  expect_s3_class(pred, "JoinMeDynPred")
})
