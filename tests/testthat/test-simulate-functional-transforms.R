test_that("simulate_joinme functional transforms support common nonlinear functions", {
  cfgs <- list(
    list(expr = ~ inv_logit(x), ref = function(x) stats::plogis(x)),
    list(expr = ~ exp(0.2 * x), ref = function(x) exp(0.2 * x)),
    list(expr = ~ power(x, 2), ref = function(x) x^2),
    list(expr = ~ sqrt(abs(x) + 1), ref = function(x) sqrt(abs(x) + 1)),
    list(expr = ~ cbrt(x), ref = function(x) sign(x) * abs(x)^(1 / 3))
  )

  for (cfg in cfgs) {
    sim <- simulate_joinme(
      n_id = 3,
      families = c("gaussian", "gaussian"),
      n_obs_per_marker_per_id = 3,
      times_obs = seq(0, 2, length.out = 4),
      seed = 1301,
      assoc = c("cv_mean"),
      assoc_coefs = c(cv_mean = 0),
      transforms = list(
        cv_mean = list(type = "functional", expr = cfg$expr)
      )
    )

    t0 <- 0.8
    comp <- sim$helpers$assoc_components(1, t0)
    raw <- sim$helpers$assoc_components_raw(1, t0)
    expect_equal(comp$cv_mean, cfg$ref(raw$cv_mean), tolerance = 1e-8)
  }
})


test_that("simulate_joinme functional transform aliases behave identically", {
  sim_expit <- simulate_joinme(
    n_id = 3,
    families = c("gaussian", "gaussian"),
    n_obs_per_marker_per_id = 3,
    times_obs = seq(0, 2, length.out = 4),
    seed = 1302,
    assoc = c("cv_mean"),
    assoc_coefs = c(cv_mean = 0),
    transforms = list(cv_mean = list(type = "functional", expr = ~ expit(x)))
  )

  sim_sigmoid <- simulate_joinme(
    n_id = 3,
    families = c("gaussian", "gaussian"),
    n_obs_per_marker_per_id = 3,
    times_obs = seq(0, 2, length.out = 4),
    seed = 1302,
    assoc = c("cv_mean"),
    assoc_coefs = c(cv_mean = 0),
    transforms = list(cv_mean = list(type = "functional", expr = ~ sigmoid(x)))
  )

  sim_pow <- simulate_joinme(
    n_id = 3,
    families = c("gaussian", "gaussian"),
    n_obs_per_marker_per_id = 3,
    times_obs = seq(0, 2, length.out = 4),
    seed = 1302,
    assoc = c("cv_mean"),
    assoc_coefs = c(cv_mean = 0),
    transforms = list(cv_mean = list(type = "functional", expr = ~ pow(x, 2)))
  )

  t0 <- 1.0
  v_expit <- sim_expit$helpers$assoc_components(2, t0)$cv_mean
  v_sigmoid <- sim_sigmoid$helpers$assoc_components(2, t0)$cv_mean
  raw <- sim_pow$helpers$assoc_components_raw(2, t0)$cv_mean
  v_pow <- sim_pow$helpers$assoc_components(2, t0)$cv_mean

  expect_equal(v_expit, v_sigmoid, tolerance = 1e-10)
  expect_equal(v_pow, raw^2, tolerance = 1e-8)
})

test_that("simulate_joinme functional transform supports unary minus", {
  sim <- simulate_joinme(
    n_id = 3,
    families = c("gaussian", "gaussian"),
    n_obs_per_marker_per_id = 3,
    times_obs = seq(0, 2, length.out = 4),
    seed = 1303,
    assoc = c("cv_mean"),
    assoc_coefs = c(cv_mean = 0),
    transforms = list(cv_mean = list(type = "functional", expr = ~ -x))
  )

  t0 <- 0.9
  raw <- sim$helpers$assoc_components_raw(1, t0)$cv_mean
  tf <- sim$helpers$assoc_components(1, t0)$cv_mean

  expect_equal(tf, -raw, tolerance = 1e-10)
})
