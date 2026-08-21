longitudinal_only_example_data <- function() {
  design <- expand.grid(
    id = seq_len(6L),
    marker = paste0("m", seq_len(3L)),
    time = c(0, 1),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ) # balanced subject-by-marker observation schedule used by interface tests
  design$y <- 1 + 0.4 * design$time +
    0.1 * design$id +
    match(design$marker, paste0("m", seq_len(3L)))
  design
}

test_that("joinme prepares a likelihood-neutral longitudinal-only model", {
  skip_if_not_installed("cmdstanr")
  has_cmdstan <- tryCatch(
    !is.null(cmdstanr::cmdstan_version()),
    error = function(error) FALSE
  ) # whether the local environment can instantiate the compiled Stan programme
  if (!has_cmdstan) {
    skip("CmdStan is not installed.")
  }

  longitudinal_data <- longitudinal_only_example_data()
  prepared <- joinme(
    formulaLong = y ~ 1 + time + (1 + time | id) +
      (0 + (1 | id) | marker),
    dataLong = longitudinal_data,
    families = rep("gaussian", 3L),
    fit = FALSE,
    control = list(engine = "cmdstanr")
  )

  expect_s3_class(prepared, "JoiNMeStanData")
  expect_identical(prepared$stan_data$include_survival, 0L)
  expect_true(all(prepared$stan_data$S_entry == 0))
  expect_true(all(prepared$stan_data$S_event == 0))
  expect_true(all(prepared$stan_data$d_event == 0L))
  expect_true(all(prepared$stan_data$event_censor_type == 0L))
  expect_true(all(unname(prepared$make_cfg()$assoc) == 0L))
  expect_false(prepared$make_cfg()$include_survival)
  expect_equal(nrow(prepared$dataEvent), length(unique(longitudinal_data$id)))
  expect_true(all(prepared$dataEvent$.joinme_event == 0L))
  expect_match(
    paste(deparse(prepared$formulaEvent), collapse = " "),
    ".joinme_follow_up",
    fixed = TRUE
  )
})

test_that("joinme requires a complete event-process pair", {
  longitudinal_data <- longitudinal_only_example_data()
  event_data <- data.frame(
    id = seq_len(6L),
    time = 1,
    event = 0L
  ) # otherwise valid event records used to isolate pair validation

  expect_error(
    joinme(
      y ~ time,
      longitudinal_data,
      formulaEvent = survival::Surv(time, event) ~ 1,
      fit = FALSE
    ),
    "must be supplied together"
  )
  expect_error(
    joinme(
      y ~ time,
      longitudinal_data,
      dataEvent = event_data,
      fit = FALSE
    ),
    "must be supplied together"
  )
  expect_error(
    joinme(
      y ~ time,
      longitudinal_data,
      assoc = "cv_mean",
      fit = FALSE
    ),
    "Association terms require a survival process"
  )
})

test_that("longitudinal-only fits reject survival-specific requests early", {
  longitudinal_data <- longitudinal_only_example_data()
  fit_stub <- JoiNMeFit$new(
    fit = NULL,
    stan_data = list(
      include_survival = 0L,
      family_names = "gaussian"
    ),
    formulaLong = y ~ time,
    formulaEvent = survival::Surv(.joinme_follow_up, .joinme_event) ~ 1,
    formulaVCov = ~1,
    config = list(include_survival = FALSE),
    call = quote(joinme(y ~ time, longitudinal_data)),
    tmax = 1,
    dataLong = longitudinal_data,
    dataEvent = .longitudinal_only_event_scaffold(
      longitudinal_data,
      id_variable = "id",
      time_variable = "time"
    )
  )

  expect_false(.fit_includes_survival(fit_stub))
  expect_output(
    print(fit_stub),
    "Nested longitudinal mixed effects model summary"
  )
  expect_error(
    predict(
      fit_stub,
      newdataLong = longitudinal_data,
      process = "event"
    ),
    "unavailable for a longitudinal-only fit"
  )
  expect_error(
    concordance(fit_stub),
    "unavailable for a longitudinal-only fit"
  )
  expect_error(
    plot(fit_stub, type = "survival"),
    "unavailable for a longitudinal-only fit"
  )
})
