test_that("simulate_joinme_mix mirrors the fitting mixture syntax", {
  mixture_arguments <- c(
    "formulaLong",
    "formulaEvent",
    "formulaVCov",
    "formulaDist",
    "transforms",
    "families",
    "fixed_marker_weights",
    "shared_marker_weights",
    "n_clusters",
    "formulaCluster",
    "cluster_type",
    "cluster_dimensions",
    "cluster_ordering",
    "seed"
  )

  expect_true(all(
    mixture_arguments %in% names(formals(simulate_joinme_mix))
  ))
  expect_true(all(
    mixture_arguments %in% names(formals(joinme_mix))
  ))
  expect_true(all(
    names(formals(simulate_joinme)) %in%
      c(names(formals(simulate_joinme_mix)), ".mixture_specification")
  ))
  expect_error(
    simulate_joinme_mix(cluster = "subject"),
    "no longer an argument"
  )
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
    n_clusters = 2,
    formulaCluster = ~x1,
    cluster_type = "subject",
    cluster_dimensions = 1,
    class_parameters = list(
      probability = c(0.45, 0.55),
      coefficient = list(
        subject = c("class_1:x1" = -0.8)
      ),
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
  expect_identical(mixture$n_clusters, 2L)
  expect_identical(mixture$cluster_type, "subject")
  expect_equal(realised, expected, tolerance = 0.01)
  expect_equal(
    unname(rowSums(mixture$probability_by_unit$subject)),
    rep(1, 16)
  )
  expect_false(any(c("class", "latent_class") %in% names(simulation$dataLong)))
  expect_false(any(c("class", "latent_class") %in% names(simulation$dataEvent)))
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
    n_clusters = 2,
    cluster_type = c("subject", "vcov"),
    cluster_dimensions = list(
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

test_that("marker weights are not accepted as a clustering type", {
  expect_error(
    simulate_joinme_mix(
      n_id = 4,
      families = c("gaussian", "gaussian"),
      times_obs = c(0, 0.1),
      time_cens = 0.2,
      cluster_type = "marker_weight",
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
      n_clusters = 2,
      cluster_ordering = "probability",
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
    n_clusters = 2,
    formulaCluster = list(~x1, ~1),
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
    n_clusters = 2,
    cluster_type = "subject",
    cluster_dimensions = c(1, 2),
    cluster_ordering = "intercept",
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
      fixed_marker_weights = TRUE,
      cluster_type = "marker_weight",
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
    n_clusters = 2,
    formulaCluster = ~x1,
    cluster_type = "subject",
    seed = 8107,
    use_mirai = FALSE
  )

  expect_null(simulation$dataEvent)
  expect_null(simulation$truth$formulaEvent)
  expect_equal(nrow(simulation$dataLong), 8L * 2L * 3L)
  expect_length(simulation$truth$mixture$allocation$subject, 8L)
  expect_true(all(is.finite(simulation$dataLong$y)))

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
    formulaVCov = ~x1,
    n_id = 8,
    families = c("gaussian", "student_t"),
    times_obs = c(0, 0.2, 0.4),
    time_cens = 0.5,
    n_clusters = 2,
    formulaCluster = ~x1,
    cluster_type = c("subject", "vcov"),
    cluster_dimensions = list(
      subject = 1,
      vcov = 1
    ),
    seed = 8109,
    use_mirai = FALSE
  )

  prepared <- joinme_mix(
    formulaLong = simulation$truth$formulaLong,
    dataLong = simulation$dataLong,
    formulaEvent = simulation$truth$formulaEvent,
    dataEvent = simulation$dataEvent,
    formulaVCov = simulation$truth$formulaVCov,
    families = simulation$marker_info$families,
    n_clusters = simulation$truth$mixture$n_clusters,
    formulaCluster = simulation$truth$mixture$formulaCluster,
    cluster_type = simulation$truth$mixture$cluster_type,
    cluster_dimensions = simulation$truth$mixture$dimensions[
      simulation$truth$mixture$cluster_type
    ],
    assoc = simulation$truth$assoc,
    shrinkage = simulation$truth$shrinkage,
    fit = FALSE,
    seed = 8109
  )

  expect_s3_class(prepared, "JoiNMeMixStanData")
  expect_identical(
    prepared$stan_data$mixture$n_clusters,
    simulation$truth$mixture$n_clusters
  )
  expect_identical(
    prepared$stan_data$mixture$cluster_type,
    simulation$truth$mixture$cluster_type
  )
  expect_identical(
    prepared$stan_data$mixture$dimensions,
    simulation$truth$mixture$dimensions
  )
})
