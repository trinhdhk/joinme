test_that("summary reports family-aware distributional parameters with marker labels", {
  skip_on_cran()
  skip_if_not_installed("rstan")

  set.seed(431)
  sim <- simulate_joinme(
    n_id = 4,
    n_obs_per_marker_per_id = 4,
    families = c("gaussian", "negbin2", "skew_normal"),
    times_obs = seq(0, 2, length.out = 4),
    seed = 431,
    family_params = list(
      gaussian = list(sigma = 0.7),
      negbin2 = list(phi = 2.5),
      skew_normal = list(sigma = 1.0, alpha = 1.0)
    )
  )

  formulaLong <- y ~ 1 + time + x1 +
    (1 + time | id) +
    (0 + x1 + (1 + time | id) | marker)
  formulaEvent <- survival::Surv(time, event) ~ 1 + x1 + x2

  fit <- joinme(
    formulaLong = formulaLong,
    dataLong = sim$dataLong,
    formulaEvent = formulaEvent,
    dataEvent = sim$dataEvent,
    assoc = c("cv_total"),
    families = c("gaussian", "negbin2", "skew_normal"),
    transforms = list(cv_total = list(type = "identity")),
    control = list(
      engine = "rstan",
      chains = 1,
      parallel_chains = 1,
      iter_warmup = 30,
      iter_sampling = 30,
      refresh = 0,
      seed = 431
    )
  )

  sum_obj <- summary(fit, draws = 20, seed = 431)
  dist_tbl <- sum_obj$tables$distributional
  expect_true(is.data.frame(dist_tbl))

  dist_terms <- dist_tbl$term
  marker_levels <- fit$stan_data$marker_levels
  family_long <- fit$stan_data$family_long

  expected_terms <- character(0)
  forbidden_terms <- character(0)

  for (d in seq_along(marker_levels)) {
    req <- .family_distrib_params(family_long[d])
    mk <- marker_levels[d]

    sigma_term <- paste0("sigma[", mk, "]")
    nu_term <- paste0("nu_marker[", mk, "]")
    phi_term <- paste0("phi_nb_marker[", mk, "]")
    alpha_term <- paste0("alpha_skew_marker[", mk, "]")
    kappa_term <- paste0("kappa_marker[", mk, "]")
    tau_term <- paste0("tau_marker[", mk, "]")

    if ("sigma" %in% req) expected_terms <- c(expected_terms, sigma_term) else forbidden_terms <- c(forbidden_terms, sigma_term)
    if ("nu" %in% req) expected_terms <- c(expected_terms, nu_term) else forbidden_terms <- c(forbidden_terms, nu_term)
    if ("phi" %in% req) expected_terms <- c(expected_terms, phi_term) else forbidden_terms <- c(forbidden_terms, phi_term)
    if ("alpha" %in% req) expected_terms <- c(expected_terms, alpha_term) else forbidden_terms <- c(forbidden_terms, alpha_term)
    if ("kappa" %in% req) expected_terms <- c(expected_terms, kappa_term) else forbidden_terms <- c(forbidden_terms, kappa_term)
    if ("tau" %in% req) expected_terms <- c(expected_terms, tau_term) else forbidden_terms <- c(forbidden_terms, tau_term)
  }

  if (length(expected_terms) > 0) {
    expect_true(all(expected_terms %in% dist_terms))
  }
  if (length(forbidden_terms) > 0) {
    expect_false(any(forbidden_terms %in% dist_terms))
  }

  expect_false(any(grepl("_marker\\\\[[0-9]+\\\\]$", dist_terms)))

})
