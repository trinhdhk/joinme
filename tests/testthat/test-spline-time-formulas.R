test_that("standata accepts spline time terms without raw-time rescaling", {
  sim <- simulate_joinme_joint_student_t_cvtotal(
    n_id = 4,
    D = 2,
    n_t = 3,
    seed = 3101,
    include_marker_only = TRUE
  )

  formulaLong <- y ~ 1 + splines::ns(time, df = 3) + x1 +
    (0 + splines::bs(time, df = 3) | id) +
    (0 + x1 + (0 + splines::ns(time, df = 3) | id) | marker)

  sd <- joinme_standata(
    formulaLong = formulaLong,
    dataLong = sim$dataLong,
    formulaEvent = survival::Surv(time, event) ~ 1 + x1 + x2,
    dataEvent = sim$dataEvent,
    assoc = c("cv_total")
  )

  expect_s3_class(sd$design_blueprints$fixed, "JoiNMe_mm_blueprint")
  expect_true(length(sd$idx_time_beta) == 0L)
  expect_true(length(sd$idx_time_uid) == 0L)
  expect_true(length(sd$idx_time_idm) == 0L)
  expect_true(sd$P > 0L)
  expect_true(sd$R_id > 0L)
  expect_true(sd$Q_idm > 0L)

  probe <- sim$dataLong[seq_len(min(6L, nrow(sim$dataLong))), , drop = FALSE]
  probe[["time"]] <- seq(0, 1, length.out = nrow(probe))
  probe_matrix <- JoiNMe:::.mm(sd$design_blueprints$fixed, probe)
  expect_equal(ncol(probe_matrix), sd$P)
})

test_that("simulate_joinme keeps spline-expanded random-effect column labels in truth draws", {
  sim <- simulate_joinme(
    formulaLong = y ~ 1 + time +
      (1 | id) +
      (1 + splines::bs(time, knots = 0.5, degree = 2) + (1 | id) | marker),
    formulaEvent = survival::Surv(time, event) ~ 1,
    n_id = 4,
    families = rep("gaussian", 2),
    time_cens = 1,
    times_obs = seq(0, 1, length.out = 3),
    n_obs_per_marker_per_id = 3,
    assoc = c("cv_total"),
    assoc_coefs = c(cv_total = 0),
    seed = 31015,
    use_mirai = FALSE
  )

  expect_equal(
    colnames(sim$truth$re_draws$marker),
    c(
      "(Intercept)",
      "splines::bs(time, knots = 0.5, degree = 2)1",
      "splines::bs(time, knots = 0.5, degree = 2)2",
      "splines::bs(time, knots = 0.5, degree = 2)3"
    )
  )
  expect_equal(rownames(sim$truth$re_draws$marker), sim$marker_info$names)
})

test_that("fit and prediction support spline time terms in formulaLong", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  has_cmdstan <- FALSE
  tryCatch({
    has_cmdstan <- !is.null(cmdstanr::cmdstan_version())
  }, error = function(e) {
    has_cmdstan <- FALSE
  })
  if (!has_cmdstan) skip("CmdStan is not installed.")

  sim <- simulate_joinme(
    n_id = 3,
    families = rep("gaussian", 2),
    n_obs_per_marker_per_id = 3,
    times_obs = seq(0, 2, length.out = 4),
    quadrature_nodes = 7,
    seed = 3102,
    assoc = c("cv_total"),
    assoc_coefs = c(cv_total = 0.15)
  )

  formulaLong <- y ~ 1 + splines::ns(time, df = 3) + x1 +
    (0 + splines::bs(time, df = 3) | id) +
    (0 + x1 + (0 + splines::ns(time, df = 3) | id) | marker)
  formulaEvent <- survival::Surv(time, event) ~ 1 + x1 + x2

  fit <- joinme(
    formulaLong = formulaLong,
    dataLong = sim$dataLong,
    formulaEvent = formulaEvent,
    dataEvent = sim$dataEvent,
    formulaVCov = ~ 1,
    assoc = c("cv_total"),
    families = rep("gaussian", 2),
    transforms = list(cv_total = list(type = "identity")),
    control = list(
      engine = "cmdstanr",
      force_recompile = TRUE,
      quadrature_nodes = 7,
      chains = 1,
      parallel_chains = 1,
      iter_warmup = 10,
      iter_sampling = 10,
      seed = 3102,
      refresh = 0
    )
  )

  expect_s3_class(fit, "JoiNMeFit")
  expect_true(length(fit$stan_data$idx_time_beta) == 0L)

  pred <- posterior_epred(
    fit,
    newdataLong = sim$dataLong,
    newdataEvent = sim$dataEvent,
    time_start = min(sim$dataLong$time),
    times = seq(0, max(sim$dataLong$time) + 0.5, length.out = 15),
    control = list(
      engine = "cmdstanr",
      force_recompile = TRUE,
      quadrature_nodes = 7,
      n_samples = 10,
      chains = 1,
      iter_warmup = 5,
      iter_sampling = 5,
      refresh = 0
    ),
    seed = 3102
  )

  expect_s3_class(pred, "JoiNMeDynPred")
})
