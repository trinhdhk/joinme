test_that("distributional parameter names distinguish Beta, Student-t, and skew Laplace roles", {
  # Step 1: verify that each family advertises only its statistically relevant
  # distributional parameters.
  expect_equal(.family_distrib_params(2L), c("sigma", "nu"))
  expect_equal(.family_distrib_params(9L), c("sigma", "tau"))
  expect_equal(.family_distrib_params(10L), "kappa")

  # Step 2: verify that the formula parser retains the three distinct names and
  # their family scopes without conflating kappa with Student-t nu.
  specification <- .normalize_formula_dist(list(
    nu[family = student_t] ~ 1,
    kappa[family = beta] ~ 1,
    tau[family = skew_double_exponential] ~ 1
  ))

  expect_setequal(names(specification), c("nu", "kappa", "tau"))
  expect_true(inherits(specification$nu$by_family$student_t, "formula"))
  expect_true(inherits(specification$kappa$by_family$beta, "formula"))
  expect_true(inherits(
    specification$tau$by_family$skew_double_exponential,
    "formula"
  ))

  # Step 3: ensure the superseded public names are not silently retained. This
  # guards the fit, prediction, summary, and simulation interfaces against
  # diverging naming schemes.
  expect_error(
    .normalize_formula_dist(list(phi_beta ~ 1)),
    "Unknown distributional parameter"
  )
  expect_error(
    .normalize_formula_dist(list(tau_sde ~ 1)),
    "Unknown distributional parameter"
  )
  expect_error(
    .normalize_formula_dist(list(precision ~ 1)),
    "Unknown distributional parameter"
  )
  expect_error(
    .normalize_formula_dist(list(skew_sde ~ 1)),
    "Unknown distributional parameter"
  )
})

test_that("Beta simulation uses the mean and kappa sample-size parameterisation", {
  # Step 1: choose a mean and sample size for which the theoretical moments are
  # well separated from common shape/precision implementation errors.
  mu <- 0.25
  kappa <- 20

  # Step 2: draw enough observations for stable moment comparisons while keeping
  # the test independent of Stan compilation and posterior sampling.
  set.seed(2507)
  y <- .sample_from_family(
    n = 50000,
    mu = mu,
    family = 10L,
    kappa = kappa
  )

  # Step 3: compare empirical moments with the Beta mean/sample-size identities:
  # E(Y) = mu and Var(Y) = mu * (1 - mu) / (kappa + 1).
  expect_equal(mean(y), mu, tolerance = 0.003)
  expect_equal(
    stats::var(y),
    mu * (1 - mu) / (kappa + 1),
    tolerance = 0.0003
  )
})

test_that("skew-double-exponential tau retains its bounded interpretation", {
  # Step 1: use the symmetric value so the sample should be centred at mu.
  set.seed(2508)
  y <- .sample_from_family(
    n = 50000,
    mu = 1.5,
    family = 9L,
    sigma = 0.7,
    tau = 0.5
  )

  # Step 2: verify the statistical symmetry implied by tau = 0.5.
  expect_equal(stats::median(y), 1.5, tolerance = 0.02)
  expect_true(all(is.finite(y)))
})

test_that("standata wires kappa and tau through separate family-scoped regressions", {
  # Step 1: construct a small deterministic mixed-family data set. Beta outcomes
  # lie strictly inside (0, 1), while skew-double-exponential outcomes may lie
  # anywhere on the real line.
  data_long <- expand.grid(
    id = 1:4,
    marker = c("beta_marker", "sde_marker"),
    time = c(0, 1)
  )
  data_long <- data_long[order(
    data_long$id,
    data_long$marker,
    data_long$time
  ), ]
  data_long$x1 <- rep(c(-0.5, 0.5, -0.25, 0.25), each = 4)
  data_long$y <- ifelse(
    data_long$marker == "beta_marker",
    rep(c(0.25, 0.35), 4),
    rep(c(-0.2, 0.2), 4)
  )
  data_event <- data.frame(
    id = 1:4,
    time = c(1.5, 2, 2.5, 3),
    event = c(1, 0, 1, 0),
    x1 = c(-0.5, 0.5, -0.25, 0.25),
    x2 = c(0, 1, 0, 1)
  )

  # Step 2: build family-scoped regressions using the public parameter names.
  standata <- joinme_standata(
    formulaLong = y ~ 1 + time + x1 +
      (1 + time | id) +
      (0 + x1 + (1 + time | id) | marker),
    dataLong = data_long,
    formulaEvent = survival::Surv(time, event) ~ 1 + x1 + x2,
    dataEvent = data_event,
    families = c("beta", "skew_double_exponential"),
    formulaDist = list(
      kappa[family = beta] ~ 1 + time,
      tau[family = skew_double_exponential] ~ 1 + x1
    ),
    assoc = "cv_total"
  )

  # Step 3: verify separate design matrices and marker-to-family mappings. This
  # is the schema consumed unchanged by fitting and dynamic prediction Stan code.
  expect_equal(standata$P_kappa, 2L)
  expect_equal(standata$P_tau, 2L)
  expect_equal(standata$marker_to_kappa_family, c(1L, 0L))
  expect_equal(standata$marker_to_tau_family, c(0L, 1L))
  expect_true(all(c("kappa", "tau") %in% names(standata$dist_cols)))
  expect_false(any(c("P_phi_beta", "P_tau_sde") %in% names(standata)))

  # Step 4: reject simultaneous fixed and regressed tau specifications because
  # only one can define the likelihood's quantile/asymmetry parameter.
  expect_error(
    joinme_standata(
      formulaLong = y ~ 1 + time + x1 +
        (1 + time | id) +
        (0 + x1 + (1 + time | id) | marker),
      dataLong = data_long,
      formulaEvent = survival::Surv(time, event) ~ 1 + x1 + x2,
      dataEvent = data_event,
      families = c("beta", "skew_double_exponential"),
      formulaDist = list(tau[family = skew_double_exponential] ~ 1),
      tau_fixed = 0.3,
      assoc = "cv_total"
    ),
    "either.*tau_fixed.*tau"
  )
})
