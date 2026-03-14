test_that("joinme_tf validates and normalises shorthand declarations", {
  tf <- joinme_tf(
    cv_total = "identity",
    corr = ~ -x,
    cv_marker = list(type = "pwlin", x = c(-1, 0, 1), y = c(0.2, 1, 0.2))
  )

  expect_s3_class(tf, "joinme_tf")
  expect_identical(tf$cv_total$type, "identity")
  expect_identical(tf$corr$type, "functional")
  expect_true(inherits(tf$corr$expr, "formula"))
  expect_identical(tf$cv_marker$type, "pwlin")

  sd_tf <- build_standata_transforms(tf)
  expect_identical(sd_tf$tf_mode_cv_tot, 0)
  expect_identical(sd_tf$tf_mode_corr, 1)
  expect_identical(sd_tf$tf_mode_cv_marker, 3)
})

test_that("joinme_tf supports expit-based spline transforms", {
  expit_grid <- stats::plogis(seq(-4, 4, length.out = 21))
  tf <- joinme_tf(
    corr = list(
      type = "ispline_expit_penalised",
      x = expit_grid,
      n_knots = 5,
      degree = 2,
      lambda = 0.5
    )
  )

  expect_s3_class(tf, "joinme_tf")
  expect_identical(tf$corr$type, "ispline_expit_penalised")

  sd_tf <- build_standata_transforms(tf)
  expect_identical(sd_tf$tf_mode_corr, 4L)
  expect_true(all(sd_tf$knots_corr > 0 & sd_tf$knots_corr < 1))
  expect_equal(length(sd_tf$coeff_corr), sd_tf$n_coeff_corr)
})

test_that("joinme_tf rejects raw-scale spline inputs for expit-based transforms", {
  expect_error(
    joinme_tf(
      corr = list(
        type = "ispline_expit_penalised",
        x = seq(-4, 4, length.out = 21),
        n_knots = 5,
        degree = 2,
        lambda = 0.5
      )
    ),
    "must lie on the expit scale"
  )
})

test_that("joinme_tf rejects unknown channels", {
  expect_error(
    joinme_tf(not_a_channel = "identity"),
    "Unknown transform channel"
  )
})

test_that("joinme_priors validates exposed prior components", {
  pri <- joinme_priors(beta = list(scale = 2.5), alpha = 1.0, lkj = 2)

  expect_s3_class(pri, "joinme_priors")
  expect_equal(pri$beta$scale, 2.5)
  expect_equal(pri$alpha, 1.0)
  expect_equal(pri$lkj, 2)

  stan_pri <- .build_priors(beta_prior = pri$beta, alpha_prior = pri$alpha, lkj_prior = pri$lkj)
  expect_equal(stan_pri$alpha_scale, 1.0)
  expect_equal(stan_pri$lkj_eta, 2)
})

test_that("joinme_priors rejects invalid LKJ values", {
  expect_error(
    joinme_priors(lkj = 0),
    "positive numeric scalar"
  )
})
