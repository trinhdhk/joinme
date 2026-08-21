test_that("simulate_joinme_mix mirrors the fitting mixture syntax", {
  mixture_arguments <- c(
    "formulaLong",
    "formulaEvent",
    "formulaVCov",
    "formulaDist",
    "transforms",
    "families",
    "n_classes",
    "formulaClass",
    "class_type",
    "class_dimensions",
    "class_ordering",
    "seed"
  )

  expect_true(all(
    mixture_arguments %in% names(formals(simulate_joinme_mix))
  ))
  expect_true(all(
    mixture_arguments %in% names(formals(joinme_mix))
  ))
  expect_true("truth" %in% names(formals(simulate_joinme_mix)))
  expect_true("priors" %in% names(formals(joinme_mix)))
  expect_false("priors" %in% names(formals(simulate_joinme_mix)))
  expect_true(all(
    names(formals(simulate_joinme)) %in%
      c(names(formals(simulate_joinme_mix)), ".mixture_specification")
  ))
})

test_that("subject mixture draws follow their recorded class parameters", {
  simulation <- simulate_joinme_mix(
    formulaLong = y ~ 1 + time + x1 +
      (1 + time | id) +
      (0 + x1 + (1 + time | id) | marker),
    formulaEvent = survival::Surv(time, event) ~ x1,
    n_id = 16,
    families = c("gaussian", "gaussian"),
    times_obs = c(0, 0.25),
    time_cens = 0.5,
    n_classes = 2,
    formulaClass = ~x1,
    class_type = "subject",
    class_dimensions = 1,
    truth = jm_truth(class = list(
      slope = c("class_1:x1" = -0.8)
    )),
    class_parameters = list(
      probability = c(0.45, 0.55),
      location = matrix(c(-2.5, 2.5), nrow = 2),
      scale = 0.001
    ),
    seed = 8101,
    use_mirai = FALSE
  )

  mixture <- simulation$truth$mixture
  allocation <- unname(mixture$allocation$subject)
  realised <- mixture$standardised_draws$subject[, 1]
  expected <- mixture$location[allocation, 1]

  expect_s3_class(simulation$dataLong, "data.frame")
  expect_s3_class(simulation$dataEvent, "data.frame")
  expect_identical(mixture$n_classes, 2L)
  expect_identical(mixture$class_type, "subject")
  expect_equal(unname(mixture$coefficient$subject), -0.8)
  expect_equal(realised, expected, tolerance = 0.01)
  expect_equal(
    unname(rowSums(mixture$probability_by_unit$subject)),
    rep(1, 16)
  )
  expect_false(any(c("class", "latent_class") %in% names(simulation$dataLong)))
  expect_false(any(c("class", "latent_class") %in% names(simulation$dataEvent)))
  expect_false("true_params" %in% names(simulation))
})

test_that("omitted mixture population parameters are drawn once and stored", {
  simulation_arguments <- list(
    n_id = 5,
    families = c("gaussian", "gaussian"),
    times_obs = c(0, 0.1),
    time_cens = 0.2,
    n_classes = 2,
    formulaClass = ~x1,
    class_type = "subject",
    class_dimensions = 1,
    use_mirai = FALSE
  )
  simulation_a <- do.call(
    simulate_joinme_mix,
    c(simulation_arguments, list(seed = 8121))
  )
  simulation_b <- do.call(
    simulate_joinme_mix,
    c(simulation_arguments, list(seed = 8122))
  )

  truth_a <- simulation_a$truth$mixture
  expect_equal(sum(truth_a$probability), 1)
  expect_true(all(is.finite(truth_a$coefficient$subject)))
  expect_true(all(is.finite(truth_a$location)))
  expect_true(all(truth_a$scale > 0))
  expect_false(identical(truth_a$location, simulation_b$truth$mixture$location))
  expect_false(identical(
    truth_a$coefficient$subject,
    simulation_b$truth$mixture$coefficient$subject
  ))
})

test_that("compatible subject blocks use one simulated allocation", {
  simulation <- simulate_joinme_mix(
    formulaLong = y ~ 1 + time +
      (1 + time | id) +
      (0 + (1 + time | id) | marker),
    formulaEvent = survival::Surv(time, event) ~ 1,
    formulaVCov = ~1,
    n_id = 14,
    families = c("gaussian", "gaussian"),
    times_obs = c(0, 0.2),
    time_cens = 0.4,
    n_classes = 2,
    class_type = c("subject", "vcov"),
    class_dimensions = list(
      subject = 1,
      vcov = 1
    ),
    class_parameters = list(
      probability = c(0.5, 0.5),
      location = matrix(
        c(-2, -1, 2, 1),
        nrow = 2,
        byrow = TRUE
      ),
      scale = 0.001
    ),
    seed = 8102,
    use_mirai = FALSE
  )

  mixture <- simulation$truth$mixture
  allocation <- unname(mixture$allocation$subject)
  subject_draw <- mixture$standardised_draws$subject[, 1]
  covariance_draw <- mixture$standardised_draws$vcov[, 1]

  expect_identical(
    mixture$allocation_domains$subject,
    c("subject", "vcov")
  )
  expect_length(allocation, 14L)
  expect_equal(
    subject_draw,
    mixture$location[allocation, mixture$starts[["subject"]]],
    tolerance = 0.01
  )
  expect_equal(
    covariance_draw,
    mixture$location[allocation, mixture$starts[["vcov"]]],
    tolerance = 0.01
  )
})

test_that("marker weights are not accepted as a class type", {
  expect_error(
    simulate_joinme_mix(
      n_id = 4,
      families = c("gaussian", "gaussian"),
      times_obs = c(0, 0.1),
      time_cens = 0.2,
      class_type = "marker_weight",
      seed = 8103,
      use_mirai = FALSE
    ),
    "subject, marker, corr, or vcov"
  )
})

test_that("mixture simulation validates ordering and context-sensitive defaults", {
  expect_error(
    simulate_joinme_mix(
      n_id = 4,
      families = c("gaussian", "gaussian"),
      times_obs = c(0, 0.1),
      time_cens = 0.2,
      n_classes = 2,
      class_ordering = "probability",
      class_parameters = list(probability = c(0.7, 0.3)),
      seed = 8104,
      use_mirai = FALSE
    ),
    "strictly increasing"
  )

  class_specific <- simulate_joinme_mix(
    n_id = 6,
    families = c("gaussian", "gaussian"),
    times_obs = c(0, 0.1),
    time_cens = 0.2,
    n_classes = 2,
    formulaClass = list(~x1, ~1),
    class_parameters = list(scale = 0.5),
    seed = 8105,
    use_mirai = FALSE
  )
  expect_identical(class_specific$truth$mixture$ordering, "none")

  intercept_ordered <- simulate_joinme_mix(
    n_id = 6,
    families = c("gaussian", "gaussian"),
    times_obs = c(0, 0.1),
    time_cens = 0.2,
    n_classes = 2,
    class_type = "subject",
    class_dimensions = c(1, 2),
    class_ordering = "intercept",
    class_parameters = list(
      location = matrix(
        c(-1, 1, 1, -1),
        nrow = 2,
        byrow = TRUE
      ),
      scale = 0.5
    ),
    seed = 8106,
    use_mirai = FALSE
  )
  expect_gt(
    diff(intercept_ordered$truth$mixture$location[, 1]),
    0
  )
  expect_lt(
    diff(intercept_ordered$truth$mixture$location[, 2]),
    0
  )

  expect_error(
    simulate_joinme_mix(
      n_id = 4,
      families = c("gaussian", "gaussian"),
      times_obs = c(0, 0.1),
      time_cens = 0.2,
      class_type = "marker_weight",
      seed = 8110,
      use_mirai = FALSE
    ),
    "subject, marker, corr, or vcov"
  )
})

test_that("simulate_joinme_mix supports a longitudinal-only mixture", {
  simulation <- simulate_joinme_mix(
    formulaLong = y ~ 1 + time + x1 +
      (1 + time | id) +
      (0 + x1 + (1 + time | id) | marker),
    formulaEvent = NULL,
    n_id = 8,
    families = c("gaussian", "gaussian"),
    times_obs = c(0, 0.2, 0.4),
    time_cens = 0.5,
    n_classes = 2,
    formulaClass = ~x1,
    class_type = "subject",
    seed = 8107,
    use_mirai = FALSE
  )

  expect_null(simulation$dataEvent)
  expect_null(simulation$truth$formulaEvent)
  expect_equal(nrow(simulation$dataLong), 8L * 2L * 3L)
  expect_length(simulation$truth$mixture$allocation$subject, 8L)
  expect_true(all(is.finite(simulation$dataLong$y)))
  expect_length(simulation$truth$assoc, 0L)
  expect_length(simulation$truth$recovery$arguments$assoc, 0L)
  recovered_specification <- do.call(
    joinme_mix,
    c(simulation$truth$recovery$arguments, list(fit = FALSE))
  )
  expect_s3_class(recovered_specification, "JoiNMeMixStanData")
  expect_identical(recovered_specification$stan_data$include_survival, 0L)

  expect_error(
    simulate_joinme_mix(
      formulaEvent = NULL,
      assoc = "cv_total",
      n_id = 4,
      families = c("gaussian", "gaussian"),
      times_obs = c(0, 0.1),
      time_cens = 0.2,
      seed = 8108,
      use_mirai = FALSE
    ),
    "requires a survival process"
  )
})

test_that("simulated data enter joinme_mix through the matching syntax", {
  simulation <- simulate_joinme_mix(
    formulaLong = y ~ 1 + time + x1 +
      (1 + time | id) +
      (0 + x1 + (1 + time | id) | marker),
    formulaEvent = survival::Surv(time, event) ~ x1 + x2,
    formulaVCov = list(sd = ~ x1, corr = ~ x2),
    n_id = 8,
    families = c("gaussian", "student_t"),
    times_obs = c(0, 0.2, 0.4),
    time_cens = 0.5,
    n_classes = 2,
    formulaClass = ~x1,
    class_type = c("subject", "vcov"),
    class_dimensions = list(
      subject = c(1, 2),
      vcov = c(1, 2)
    ),
    seed = 8109,
    use_mirai = FALSE
  )

  prepared <- do.call(
    joinme_mix,
    c(
      simulation$truth$recovery$arguments,
      list(fit = FALSE, seed = 8109)
    )
  )

  expect_s3_class(prepared, "JoiNMeMixStanData")
  expect_identical(prepared$stan_data$K_cov_sd, 1L)
  expect_identical(prepared$stan_data$K_cov_corr, 1L)
  expect_identical(simulation$truth$recovery$entry_point, "joinme_mix")
  expect_equal(
    simulation$truth$recovery$arguments$formulaVCov,
    simulation$truth$formulaVCov
  )
  expect_equal(
    simulation$truth$recovery$arguments$formulaClass,
    simulation$truth$mixture$formulaClass
  )
  expect_identical(
    simulation$truth$recovery$arguments$class_type,
    simulation$truth$mixture$class_type
  )
  expect_identical(
    prepared$stan_data$mixture$n_classes,
    simulation$truth$mixture$n_classes
  )
  expect_identical(
    prepared$stan_data$mixture$class_type,
    simulation$truth$mixture$class_type
  )
  expect_identical(
    prepared$stan_data$mixture$dimensions,
    simulation$truth$mixture$dimensions
  )
  simulation_contract <-
    simulation$truth$mixture$stan_data_contract
  expect_true(is.list(simulation_contract))
  expect_true(length(simulation_contract) > 0L)
  for (field_name in names(simulation_contract)) {
    expect_true(
      field_name %in% names(prepared$stan_data),
      info = paste("missing fitted Stan field", field_name)
    )
    expect_equal(
      prepared$stan_data[[field_name]],
      simulation_contract[[field_name]],
      ignore_attr = TRUE,
      info = paste("simulation/fit mismatch for", field_name)
    )
  }
})
