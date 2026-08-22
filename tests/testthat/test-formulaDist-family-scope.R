test_that("normalize_formula_dist parses family-scoped syntax", {
  spec <- .normalize_formula_dist(list(
    sigma[family = "student"] ~ 1 + time,
    sigma[family = normal] ~ 1 + x1,
    alpha_skew[family = skew_normal] ~ 1,
    kappa ~ 1
  ))

  expect_true(.is_dist_scope(spec$sigma))
  expect_setequal(names(spec$sigma$by_family), c("student_t", "gaussian"))
  expect_true(inherits(spec$sigma$by_family$student_t, "formula"))
  expect_true(inherits(spec$alpha$by_family$skew_normal, "formula"))
  expect_true(inherits(spec$kappa$default, "formula"))
})


test_that("joinme_standata builds family-gated distributional matrices", {
  set.seed(1901)
  sim <- simulate_joinme(
    n_id = 6,
    families = c("gaussian", "student_t"),
    times_obs = seq(0, 2, length.out = 5),
    seed = 1901,
    assoc = c("cv_total"),
    truth = jm_truth(assoc_coef = list(slope = c(cv_total = 0.2)))
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

test_that("joinme_standata assigns family- and marker-scoped distributional priors", {
  set.seed(19011)
  sim <- simulate_joinme(
    n_id = 6,
    families = c("gaussian", "student_t"),
    times_obs = seq(0, 2, length.out = 5),
    seed = 19011,
    assoc = "cv_total"
  )

  common_arguments <- list(
    formulaLong = y ~ 1 + time + x1 +
      (1 + time | id) +
      (0 + x1 + (1 + time | id) | marker),
    dataLong = sim$dataLong,
    formulaEvent = survival::Surv(time, event) ~ 1 + x1 + x2,
    dataEvent = sim$dataEvent,
    families = c("gaussian", "student_t"),
    formulaDist = list(
      sigma[family = gaussian] ~ 1 + x1,
      sigma[family = student_t] ~ 1 + time,
      nu[family = student_t] ~ 1
    ),
    assoc = "cv_total"
  )

  family_data <- do.call(joinme_standata, c(common_arguments, list(
    prior_specification = jm_prior(
      `sigma[family='student']` = list(
        intercept = prior_laplace(scale = 0.7),
        slope = prior_normal(scale = 0.2)
      )
    )
  )))
  sigma_positions <- family_data$prior_start_sigma + seq_len(family_data$P_sigma) - 1L
  student_positions <- grepl("^family=student_t::", family_data$dist_cols$sigma)
  expect_equal(
    as.integer(family_data$prior_regression_family[sigma_positions][student_positions]),
    c(3L, 2L)
  )
  expect_equal(
    family_data$prior_regression_scale[sigma_positions][student_positions],
    c(0.7, 0.2)
  )

  marker_data <- do.call(joinme_standata, c(common_arguments, list(
    prior_specification = jm_prior(
      `nu[marker='m2']` = list(intercept = prior_laplace(scale = 0.3))
    )
  )))
  nu_position <- marker_data$prior_start_nu + seq_len(marker_data$P_nu) - 1L
  expect_equal(as.integer(marker_data$prior_regression_family[nu_position]), 3L)
  expect_equal(marker_data$prior_regression_scale[nu_position], 0.3)
})

test_that("marker-scoped priors reject coefficients shared by several markers", {
  set.seed(19012)
  sim <- simulate_joinme(
    n_id = 5,
    families = c("student_t", "student_t"),
    times_obs = c(0, 1, 2),
    seed = 19012,
    assoc = "cv_total"
  )

  expect_error(
    joinme_standata(
      formulaLong = y ~ 1 + time + x1 +
        (1 + time | id) +
        (0 + x1 + (1 + time | id) | marker),
      dataLong = sim$dataLong,
      formulaEvent = survival::Surv(time, event) ~ 1 + x1 + x2,
      dataEvent = sim$dataEvent,
      families = c("student_t", "student_t"),
      formulaDist = list(sigma[family = student_t] ~ 1),
      prior_specification = jm_prior(
        `sigma[marker='m1']` = list(intercept = prior_normal())
      ),
      assoc = "cv_total"
    ),
    "share.*student_t.*distributional coefficients"
  )
})

test_that("family-scoped horseshoes retain separate global regularisation groups", {
  pri <- jm_prior(
    `sigma[family='gaussian']` = list(
      slope = prior_horseshoe(global_scale = 0.1)
    ),
    `sigma[family='student']` = list(
      slope = prior_horseshoe(global_scale = 0.2)
    )
  )
  block <- .distributional_prior_block(
    parameter = "sigma",
    columns = c(
      "family=gaussian::(Intercept)", "family=gaussian::x1",
      "family=student_t::(Intercept)", "family=student_t::time"
    ),
    roles = c("intercept", "slope", "intercept", "slope"),
    specification = pri$distributional$sigma,
    marker_levels = c("m1", "m2"),
    family_codes = c(1L, 2L)
  )
  packed <- .pack_regression_priors(list(sigma = block))

  expect_equal(
    as.integer(packed$prior_regression_horseshoe_group_index),
    c(0L, 1L, 0L, 2L)
  )
  expect_equal(packed$prior_regression_horseshoe_global_scale, c(0.1, 0.2))
})


test_that("joinme_standata rejects unsupported parameter-family scopes", {
  set.seed(1902)
  sim <- simulate_joinme(
    n_id = 5,
    families = c("gaussian", "student_t"),
    times_obs = seq(0, 2, length.out = 5),
    seed = 1902,
    assoc = c("cv_total"),
    truth = jm_truth(assoc_coef = list(slope = c(cv_total = 0.2)))
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
    times_obs = seq(0, 2, length.out = 5),
    seed = 1903,
    formulaDist = list(
      sigma[family = gaussian] ~ 1 + x1,
      sigma[family = student_t] ~ 1 + time,
      alpha_skew[family = skew_normal] ~ 1
    ),
    truth = jm_truth(
      `sigma[family='student']` = list(slope = prior_normal(scale = 0.25))
    )
  )

  expect_true(is.data.frame(sim$dataLong))
  expect_true(is.data.frame(sim$dataEvent))
  expect_gt(nrow(sim$dataLong), 0)
  expect_equal(
    sim$truth$priors$distributional$sigma$by_family$student_t$slope$scale,
    0.25
  )
})

test_that("simulate_joinme supports family-specific dist_coefs for scoped formulas", {
  set.seed(19031)
  sim <- simulate_joinme(
    n_id = 120,
    families = c("gaussian", "skew_normal"),
    times_obs = seq(0, 2, length.out = 8),
    seed = 19031,
    formulaDist = list(
      sigma[family = gaussian] ~ 1,
      sigma[family = skew_normal] ~ 1
    ),
    truth = jm_truth(
      longitudinal = c("(Intercept)" = 0, "time" = 0, "x1" = 0),
      `sigma[family='gaussian']` = list(intercept = c("(Intercept)" = log(0.15))),
      `sigma[family='skew_normal']` = list(intercept = c("(Intercept)" = log(1.10))),
      vcov = list(
        sd = list(intercept = rep(log(expm1(1e-8)), 2), latent = 0),
        corr = list(intercept = 0, latent = 0)
      ),
      re_params = list(
        id = list(sd = c(1e-8, 1e-8)),
        marker = list(sd = 1e-8)
      )
    )
  )

  y_sd_by_marker <- tapply(sim$dataLong$y, sim$dataLong$marker, stats::sd)
  expect_true(all(is.finite(y_sd_by_marker)))
  expect_true(unname(y_sd_by_marker["m2"]) > unname(y_sd_by_marker["m1"]) * 3)
})


test_that("JoiNMe fits and predicts with family-scoped formulaDist", {
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
    times_obs = seq(0, 2, length.out = 5),
    seed = 1904,
    assoc = c("cv_total"),
    truth = jm_truth(assoc_coef = list(slope = c(cv_total = 0.2)))
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

  expect_s3_class(fit, "JoiNMeFit")
  expect_true(fit$stan_data$P_sigma >= 2)
  expect_true(any(grepl("^family=student_t::", fit$stan_data$dist_cols$sigma)))

  ndL <- sim$dataLong[sim$dataLong$id == 1, ]
  ndE <- sim$dataEvent[sim$dataEvent$id == 1, ]

  pred <- posterior_epred(
    fit,
    newdataLong = ndL,
    newdataEvent = ndE,
    time_start = max(ndL$time),
    control = list(
      n_samples = 10,
      n_times = 50,
      engine = "cmdstanr",
      chains = 1,
      iter_warmup = 20,
      iter_sampling = 1,
      refresh = 0
    ),
    seed = 1905
  )

  expect_s3_class(pred, "JoiNMeDynPred")
})
