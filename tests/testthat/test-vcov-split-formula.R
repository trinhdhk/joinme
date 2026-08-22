test_that("formulaVCov builds independent SD and correlation designs", {
  subject_data <- data.frame(
    id = 1:4,
    treatment = c(0, 1, 0, 1),
    age = c(50, 60, 55, 65)
  )

  formulae <- .get_vcov_formula(
    list(sd = ~ treatment, corr = ~ age),
    context = "test"
  )
  design <- .build_vcov_design(
    formulae,
    dataEvent = subject_data,
    time_var = "time",
    context = "test"
  )

  expect_identical(names(formulae), c("sd", "corr"))
  expect_identical(design$K_cov_sd, 1L)
  expect_identical(design$K_cov_corr, 1L)
  expect_equal(as.numeric(design$Xcov_sd[, 1]), subject_data$treatment)
  expect_equal(as.numeric(design$Xcov_corr[, 1]), subject_data$age)
  expect_identical(colnames(design$Xcov_sd), "treatment")
  expect_identical(colnames(design$Xcov_corr), "age")
})

test_that("one formulaVCov remains an exact shared-design shorthand", {
  subject_data <- data.frame(id = 1:3, x = c(-1, 0, 1))
  formulae <- .get_vcov_formula(~ x, context = "test")
  design <- .build_vcov_design(
    formulae,
    dataEvent = subject_data,
    time_var = "time",
    context = "test"
  )

  expect_equal(as.character(formulae$sd), as.character(~ x))
  expect_equal(as.character(formulae$corr), as.character(~ x))
  expect_equal(design$Xcov_sd, design$Xcov_corr)
})

test_that("stored covariance designs retain read-only compatibility", {
  shared_matrix <- cbind(treatment = c(0, 1), age = c(50, 60))
  earlier <- .stored_vcov_design(list(
    K_cov = 2L,
    Xcov = shared_matrix,
    n_id = 2L
  ))
  expect_true(earlier$shared_format)
  expect_identical(earlier$k_sd, 2L)
  expect_identical(earlier$k_corr, 2L)
  expect_equal(earlier$x_sd, shared_matrix)
  expect_equal(earlier$x_corr, shared_matrix)

  current <- .stored_vcov_design(list(
    K_cov_sd = 1L,
    Xcov_sd = shared_matrix[, 1L, drop = FALSE],
    K_cov_corr = 1L,
    Xcov_corr = shared_matrix[, 2L, drop = FALSE],
    n_id = 2L
  ))
  expect_false(current$shared_format)
  expect_identical(colnames(current$x_sd), "treatment")
  expect_identical(colnames(current$x_corr), "age")
})

test_that("formulaVCov rejects partial or unknown component lists", {
  expect_error(
    .get_vcov_formula(list(sd = ~ x), context = "test"),
    "missing corr"
  )
  expect_error(
    .get_vcov_formula(list(sd = ~ x, corr = ~ z, covariance = ~ 1), context = "test"),
    "unknown"
  )
})

test_that("covariance prior blocks are assembled in their documented order", {
  priors <- jm_priors(vcov = list(
    sd = list(intercept = prior_laplace(
      mu = c(1, 2, 3, 4), scale = c(0.5, 0.6, 0.7, 0.8)
    )),
    corr = list(intercept = prior_normal(mu = c(-1, -2), scale = c(1, 2)))
  ))
  stan_prior <- .pack_regression_priors(
    list(
      vcov_sd = list(roles = rep("intercept", 4L), priors = priors$vcov$sd),
      vcov_corr = list(roles = rep("intercept", 2L), priors = priors$vcov$corr)
    )
  )

  expect_equal(as.integer(stan_prior$prior_regression_family), c(rep(3L, 4L), rep(2L, 2L)))
  expect_equal(stan_prior$prior_regression_mu, c(1, 2, 3, 4, -1, -2))
  expect_equal(stan_prior$prior_regression_scale, c(0.5, 0.6, 0.7, 0.8, 1, 2))
})

test_that("simulation and fitting share the split covariance contract", {
  skip_on_cran()
  supplied_sd_beta <- matrix(c(0.3, -0.2), nrow = 2L, ncol = 1L)
  supplied_corr_beta <- matrix(-0.4, nrow = 1L, ncol = 1L)
  simulation <- simulate_joinme(
    formulaLong = y ~ 1 + time + x_sd +
      (1 + time | id) +
      (1 + (1 + time | id) | marker),
    formulaEvent = survival::Surv(time, event) ~ x_sd,
    formulaVCov = list(sd = ~ x_sd, corr = ~ x_corr),
    n_id = 6,
    families = c("gaussian", "gaussian"),
    marker_levels = c("m1", "m2"),
    times_obs = c(0, 1),
    covariate_formulas = list(
      x_sd ~ rnorm(n_id),
      x_corr ~ rnorm(n_id)
    ),
    truth = jm_truth(vcov = list(
      sd = list(
        intercept = c(-0.2, -0.1),
        slope = as.vector(supplied_sd_beta),
        latent = c(0.2, 0.25)
      ),
      corr = list(
        intercept = 0.15,
        slope = as.vector(supplied_corr_beta),
        latent = 0.3
      )
    )),
    vcov_diag_link = "exp",
    re_params = list(
      id = list(sd = NULL, corr = NULL),
      marker = list(sd = NULL, corr = NULL),
      dist = list()
    ),
    h0 = function(time) rep(0.02, length(time)),
    time_cens = 2,
    seed = 718,
    n_workers = 1,
    use_mirai = FALSE
  )

  expect_equal(simulation$truth$id_marker_cov_effective$sd$beta, supplied_sd_beta)
  expect_equal(simulation$truth$id_marker_cov_effective$corr$beta, supplied_corr_beta)
  expect_equal(simulation$truth$id_marker_cov_effective$sd$alpha, c(-0.2, -0.1))
  expect_equal(simulation$truth$id_marker_cov_effective$corr$alpha, 0.15)
  expect_equal(simulation$truth$stan_fit$alpha_L, c(-0.2, 0.15, -0.1))
  expect_equal(simulation$truth$stan_fit$beta_L_sd, supplied_sd_beta)
  expect_equal(simulation$truth$stan_fit$beta_L_corr, supplied_corr_beta)
  expect_equal(
    simulation$truth$stan_fit$lambda_L,
    simulation$truth$id_marker_cov_effective$lambda
  )
  expect_equal(
    simulation$truth$stan_fit$z_L,
    simulation$truth$id_marker_cov_effective$z
  )

  # Reconstruct every subject's covariance factor from the two independent
  # designs. This checks the numerical generator, not merely stored metadata.
  effective_covariance <- simulation$truth$id_marker_cov_effective
  standard_deviation_predictor <-
    matrix(c(-0.2, -0.1), nrow = nrow(simulation$dataEvent), ncol = 2L, byrow = TRUE) +
    simulation$dataEvent$x_sd %o% as.numeric(supplied_sd_beta) +
    sweep(
      effective_covariance$sd$z,
      2L,
      effective_covariance$sd$lambda,
      `*`
    )
  correlation_predictor <-
    0.15 +
    simulation$dataEvent$x_corr * as.numeric(supplied_corr_beta) +
    effective_covariance$corr$z[, 1L] * effective_covariance$corr$lambda
  expected_standard_deviation <- exp(standard_deviation_predictor)
  expected_partial_correlation <- tanh(correlation_predictor)
  expect_equal(
    simulation$truth$L_i[, 1L, 1L],
    expected_standard_deviation[, 1L],
    tolerance = 1e-10
  )
  expect_equal(
    simulation$truth$L_i[, 2L, 1L],
    expected_standard_deviation[, 2L] * expected_partial_correlation,
    tolerance = 1e-10
  )
  expect_equal(
    simulation$truth$L_i[, 2L, 2L],
    expected_standard_deviation[, 2L] *
      sqrt(1 - expected_partial_correlation^2),
    tolerance = 1e-10
  )
  expect_identical(simulation$truth$recovery$entry_point, "joinme")
  expect_equal(
    simulation$truth$recovery$arguments$formulaVCov,
    simulation$truth$formulaVCov
  )
  expect_s3_class(
    simulation$truth$recovery$arguments$priors,
    "joinme_priors"
  )

  prepared <- do.call(
    joinme,
    c(simulation$truth$recovery$arguments, list(fit = FALSE))
  )
  expect_identical(prepared$stan_data$K_cov_sd, 1L)
  expect_identical(prepared$stan_data$K_cov_corr, 1L)
  expect_identical(prepared$stan_data$P_vcov_sd, 4L)
  expect_identical(prepared$stan_data$P_vcov_corr, 2L)
})
