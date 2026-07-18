test_that("standardized shrinkage draws match all three Stan families", {
  set.seed(7301)
  got_t <- .sim_draw_standard_shrinkage(8, 0L)
  set.seed(7301)
  expect_equal(got_t, stats::rt(8, df = 6))

  set.seed(7302)
  got_laplace <- .sim_draw_standard_shrinkage(8, 1L)
  set.seed(7302)
  expected_laplace <- sample(c(-1, 1), 8, replace = TRUE) * stats::rexp(8, rate = 1)
  expect_equal(got_laplace, expected_laplace)

  set.seed(7303)
  got_normal <- .sim_draw_standard_shrinkage(8, 2L)
  set.seed(7303)
  expect_equal(got_normal, stats::rnorm(8))
})

test_that("simulate_joinme stores base, latent, and effective marker-weight truth", {
  sim <- simulate_joinme(
    n_id = 4,
    families = rep("gaussian", 3),
    times_obs = c(0, 0.5, 1),
    time_cens = 1,
    marker_weights = c(m1 = 0.5, m2 = -0.25, m3 = 1),
    fixed_marker_weights = FALSE,
    shrinkage = 1L,
    assoc = "cv_total",
    assoc_coefs = c(cv_total = 0),
    seed = 7304,
    use_mirai = FALSE
  )

  truth <- sim$truth
  expect_equal(truth$shrinkage, 1L)
  expect_equal(truth$shrinkage_distribution, "double_exponential(0, 1)")
  expect_false(truth$fixed_marker_weights)
  expect_equal(
    unname(truth$marker_weights),
    unname(truth$marker_weights_base + truth$marker_weights_latent)
  )
  expect_equal(
    unname(truth$stan_fit$z_marker_weights),
    unname(truth$marker_weights_latent)
  )
  expect_true(any(abs(truth$marker_weights_latent) > 0))
})

test_that("simulate_joinme uses zero bases for estimated weights and supplied fixed weights exactly", {
  estimated <- simulate_joinme(
    n_id = 3,
    families = rep("gaussian", 2),
    times_obs = c(0, 0.5),
    time_cens = 0.5,
    fixed_marker_weights = FALSE,
    shrinkage = 2L,
    assoc = "cv_total",
    assoc_coefs = c(cv_total = 0),
    seed = 7305,
    use_mirai = FALSE
  )
  expect_equal(unname(estimated$truth$marker_weights_base), c(0, 0))
  expect_equal(estimated$truth$shrinkage_distribution, "normal(0, 1)")

  supplied <- c(m1 = 2, m2 = -1)
  fixed <- simulate_joinme(
    n_id = 3,
    families = rep("gaussian", 2),
    times_obs = c(0, 0.5),
    time_cens = 0.5,
    marker_weights = supplied,
    fixed_marker_weights = TRUE,
    shrinkage = 0L,
    assoc = "cv_total",
    assoc_coefs = c(cv_total = 0),
    seed = 7306,
    use_mirai = FALSE
  )
  expect_equal(fixed$truth$marker_weights, supplied)
  expect_equal(unname(fixed$truth$marker_weights_latent), c(0, 0))
  expect_equal(fixed$truth$stan_fit$z_marker_weights, numeric(2))
})

test_that("term-specific simulated weights use Stan's flattened set layout", {
  sim <- simulate_joinme(
    n_id = 3,
    families = rep("gaussian", 2),
    times_obs = c(0, 0.5),
    time_cens = 0.5,
    marker_weights = list(cv_total = c(0.2, 0.4), cv_marker = c(-0.3, 0.1)),
    shared_marker_weights = FALSE,
    fixed_marker_weights = FALSE,
    shrinkage = 0L,
    assoc = c("cv_total", "cv_marker"),
    assoc_coefs = c(cv_total = 0, cv_marker = 0),
    seed = 7307,
    use_mirai = FALSE
  )

  truth <- sim$truth
  expect_equal(names(truth$marker_weights), c("cv_total", "cv_marker"))
  expect_equal(
    truth$stan_fit$z_marker_weights,
    as.numeric(t(truth$z_marker_weight_sets))
  )
  for (term in c("cv_total", "cv_marker")) {
    expect_equal(
      unname(truth$marker_weights[[term]]),
      unname(truth$marker_weights_base[[term]] + truth$marker_weights_latent[[term]])
    )
  }
})

test_that("simulate_joinme validates marker-weight simulation controls", {
  expect_error(
    simulate_joinme(n_id = 2, fixed_marker_weights = NA),
    "fixed_marker_weights"
  )
  expect_error(
    simulate_joinme(n_id = 2, shrinkage = 3L),
    "shrinkage"
  )
})
