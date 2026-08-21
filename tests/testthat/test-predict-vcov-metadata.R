test_that("prediction draw extraction accepts split covariance dimensions without K_cov", {
  fitted_draws <- posterior::as_draws_matrix(matrix(
    c(
      0.2, 0.5, 1, -0.4, 0.3, -1.2, 0.1, 0.8,
      0.3, 0.6, 1, -0.3, 0.4, -1.1, 0.2, 0.9
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c(
      "beta[1]",
      "tau_u[1]",
      "Lcorr_u[1,1]",
      "alpha_L[1]",
      "lambda_L[1]",
      "bs_gamma_c[1,1]",
      "gamma_w[1,1]",
      "sigma_family[1]"
    ))
  )) # compact fitted draws sufficient to exercise every mandatory extractor block
  fitted_object <- structure(list(
    fit = structure(list(), class = "mock_fit"),
    stan_data = list(
      P = 1L,
      R_id = 1L,
      R_mk = 0L,
      Q_idm = 1L,
      D = 1L,
      indep_idmarker_cov = 0L,
      assoc_corr = 0L,
      assoc_vcov = 0L,
      K_cov_sd = 0L,
      Xcov_sd = matrix(0, nrow = 1L, ncol = 0L),
      K_cov_corr = 0L,
      Xcov_corr = matrix(0, nrow = 1L, ncol = 0L),
      n_id = 1L,
      K_event = 1L,
      Kbs = 1L,
      p_w = 1L,
      n_family_sigma = 1L,
      n_family_nu = 0L,
      n_family_phi = 0L,
      n_family_alpha = 0L,
      n_family_kappa = 0L,
      n_family_tau = 0L
    ) # current covariance metadata deliberately omits the superseded K_cov field
  ), class = "JoiNMeFit")

  testthat::local_mocked_bindings(
    .get_draws_obj = function(fit, variables = NULL, draws = NULL, seed = 1, keep_chains = FALSE) {
      fitted_draws
    },
    .package = "joinme"
  )

  extracted <- .extract_draws_for_pred(
    fitted_object,
    n_samples = 2L,
    seed = 4409
  ) # lower-level prediction input recovered from current fit metadata

  expect_equal(dim(extracted$beta_vcov_sd_flat), c(2L, 0L))
  expect_equal(dim(extracted$beta_vcov_corr_flat), c(2L, 0L))
  expect_equal(dim(extracted$alpha_vcov_reg), c(2L, 1L))
})
