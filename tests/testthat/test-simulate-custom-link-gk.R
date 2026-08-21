test_that("simulate_joinme supports custom inverse link and GK nodes", {
  skip_on_cran()

  sim <- simulate_joinme(
    n_id = 3,
    families = list(jm_family("bernoulli", inv_link = ~ inv_logit(x / 2))),
    times_obs = seq(0, 2, length.out = 3),
    quadrature_nodes = 31,
    integration_control = list(method = "gk", subdivisions = 8L, rel.tol = 1e-6),
    seed = 101
  )

  expect_equal(sim$truth$quadrature_nodes, 31L)
  expect_equal(length(sim$truth$gk_nodes), 31L)
  expect_equal(length(sim$truth$gk_weights), 31L)
  expect_equal(sum(sim$truth$gk_weights), 1, tolerance = 1e-8)
  expect_true(sim$truth$link_names[1] %in% c("custom", "logit"))
})
