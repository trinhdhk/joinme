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

test_that("joinme_tf records fit-only affine shift flags for supported functional transforms", {
  tf <- joinme_tf(cv_total = ~ SoftMax(x, intercept = TRUE, slope = TRUE))

  expect_identical(tf$cv_total$type, "functional")
  expect_identical(as.character(tf$cv_total$expr[[2]][[1]]), "softmax")
  expect_true(isTRUE(tf$cv_total$estimate_iota_intercept))
  expect_true(isTRUE(tf$cv_total$estimate_iota_slope))

  sd_tf <- build_standata_transforms(tf)
  expect_identical(sd_tf$estimate_iota_intercept_cv, 1L)
  expect_identical(sd_tf$estimate_iota_slope_cv, 1L)
})

test_that("joinme_tf tracks multiple affine-shifted nonlinear nodes", {
  tf <- joinme_tf(
    cv_total = ~ softplus(expit(x, intercept = TRUE, slope = TRUE), intercept = TRUE, slope = TRUE)
  )

  expect_identical(tf$cv_total$n_iota_intercept, 2L)
  expect_identical(tf$cv_total$n_iota_slope, 2L)
  expect_length(tf$cv_total$iota_nodes, 2L)
  expect_identical(vapply(tf$cv_total$iota_nodes, `[[`, character(1), "path"), c("root/1", "root"))

  sd_tf <- build_standata_transforms(tf)
  expect_identical(sd_tf$estimate_iota_intercept_cv, 2L)
  expect_identical(sd_tf$estimate_iota_slope_cv, 2L)
  expect_equal(sd_tf$functional_iota_intercept_idx_cv, c(0L, 1L, 2L))
  expect_equal(sd_tf$functional_iota_slope_idx_cv, c(0L, 1L, 2L))
})

test_that("joinme_tf rejects affine shift flags outside supported functional transforms", {
  expect_error(
    joinme_tf(cv_total = list(type = "pwlin", x = c(-1, 0, 1), y = c(0, 0.5, 1), intercept = TRUE)),
    "Fit-only affine-shift options"
  )
})

test_that("joinme_priors validates exposed prior components", {
  pri <- joinme_priors(beta = list(scale = 2.5), alpha = 1.0, iota = 0.75, lkj = 2)

  expect_s3_class(pri, "joinme_priors")
  expect_equal(pri$beta$scale, 2.5)
  expect_equal(pri$alpha, 1.0)
  expect_equal(pri$iota, 0.75)
  expect_equal(pri$lkj, 2)

  stan_pri <- .build_priors(beta_prior = pri$beta, alpha_prior = pri$alpha, iota_prior = pri$iota, lkj_prior = pri$lkj)
  expect_equal(stan_pri$alpha_scale, 1.0)
  expect_equal(stan_pri$iota_scale, 0.75)
  expect_equal(stan_pri$lkj_eta, 2)
})

test_that("joinme_priors rejects invalid LKJ values", {
  expect_error(
    joinme_priors(lkj = 0),
    "positive numeric scalar"
  )
})
