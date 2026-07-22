test_that("standata supports Surv(start, stop, event) interval rows", {
  sim <- simulate_joinme(
    n_id = 8,
    seed = 123,
    formulaEvent = survival::Surv(time_start, time_stop, event) ~ x1 + x2,
    covariate_formulas = list(
      x1 ~ time_varyring(rnorm, c(2, 4), mean = 0, sd = 1),
      x2 ~ rnorm(n_id)
    ),
    left_truncation_max = 2
  )

  sd <- joinme_standata(
    formulaLong = y ~ 1 + time + x1 +
      (1 + time | id) +
      (0 + x1 + (1 + time | id) | marker),
    dataLong = sim$dataLong,
    formulaEvent = survival::Surv(time_start, time_stop, event) ~ x1 + x2,
    dataEvent = sim$dataEvent,
    assoc = c("cv_total")
  )

  expect_true(sd$N_event >= sd$n_id)
  expect_equal(length(sd$event_id), sd$N_event)
  expect_equal(length(sd$event_start_idx), sd$n_id)
  expect_equal(length(sd$event_end_idx), sd$n_id)
  expect_true(all(sd$S_event > sd$S_entry))
  expect_true(all(sd$event_start_idx <= sd$event_end_idx))
  expect_true(any(sd$S_entry > 0))
  expect_true(anyDuplicated(sim$dataEvent$id) > 0)
  expect_true("x1" %in% names(sim$dataEvent))
})

test_that("time_varyring supports exact steps breakpoints", {
  sim <- simulate_joinme(
    n_id = 8,
    seed = 909,
    formulaEvent = survival::Surv(time_start, time_stop, event) ~ x + x2,
    covariate_formulas = list(
      x ~ time_varyring(rnorm, steps = c(1L, 3L, 5L), mean = 0, sd = 1),
      x2 ~ rnorm(n_id)
    ),
    time_cens = 6,
    left_truncation_max = 0.2
  )

  expect_true(anyDuplicated(sim$dataEvent$id) > 0)
  expect_true("x" %in% names(sim$dataEvent))

  # Split boundaries should include designated steps when they are inside each
  # subject interval; steps beyond stop time are naturally ignored.
  has_step_boundary <- any(abs(sim$dataEvent$time_start - 1) < 1e-8) ||
    any(abs(sim$dataEvent$time_stop - 1) < 1e-8) ||
    any(abs(sim$dataEvent$time_start - 3) < 1e-8) ||
    any(abs(sim$dataEvent$time_stop - 3) < 1e-8) ||
    any(abs(sim$dataEvent$time_start - 5) < 1e-8) ||
    any(abs(sim$dataEvent$time_stop - 5) < 1e-8)
  expect_true(has_step_boundary)
  expect_true(all(sim$dataEvent$time_stop > sim$dataEvent$time_start))
})

test_that("time_varyring rejects negative steps", {
  expect_error(
    simulate_joinme(
      n_id = 6,
      seed = 910,
      formulaEvent = survival::Surv(time_start, time_stop, event) ~ x + x2,
      covariate_formulas = list(
        x ~ time_varyring(rnorm, steps = c(1, -2, 4), mean = 0, sd = 1),
        x2 ~ rnorm(n_id)
      )
    ),
    regexp = "steps.*positive integers",
    fixed = FALSE
  )
})

test_that("standata supports Surv left and interval2 censoring", {
  sim <- simulate_joinme(
    n_id = 10,
    seed = 444
  )

  dataEvent_subject <- sim$dataEvent

  dataEvent_left <- dataEvent_subject
  dataEvent_left$status_left <- as.integer(dataEvent_left$event)

  sd_left <- joinme_standata(
    formulaLong = y ~ 1 + time + x1 +
      (1 + time | id) +
      (0 + x1 + (1 + time | id) | marker),
    dataLong = sim$dataLong,
    formulaEvent = survival::Surv(time, status_left, type = "left") ~ x1 + x2,
    dataEvent = dataEvent_left,
    assoc = c("cv_total")
  )

  expect_true(all(sd_left$event_censor_type %in% c(0L, 2L)))

  dataEvent_int2 <- dataEvent_subject
  dataEvent_int2$time1 <- ifelse(dataEvent_int2$event == 1L, pmax(dataEvent_int2$time - 0.2, 0), dataEvent_int2$time)
  dataEvent_int2$time2 <- ifelse(dataEvent_int2$event == 1L, dataEvent_int2$time, NA_real_)
  dataEvent_int2$time1[1] <- NA_real_
  dataEvent_int2$time2[1] <- max(0.2, min(dataEvent_int2$time[1], 1.0))

  sd_int2 <- joinme_standata(
    formulaLong = y ~ 1 + time + x1 +
      (1 + time | id) +
      (0 + x1 + (1 + time | id) | marker),
    dataLong = sim$dataLong,
    formulaEvent = survival::Surv(time1, time2, type = "interval2") ~ x1 + x2,
    dataEvent = dataEvent_int2,
    assoc = c("cv_total")
  )

  expect_true(all(sd_int2$event_censor_type %in% c(0L, 1L, 2L, 3L)))
})

test_that("legacy Surv interval type is rejected", {
  sim <- simulate_joinme(n_id = 6, seed = 777)

  de <- sim$dataEvent
  de$time1 <- ifelse(de$event == 1L, pmax(de$time - 0.3, 0), de$time)
  de$time2 <- ifelse(de$event == 1L, de$time, de$time)
  de$status <- ifelse(de$event == 1L, 3L, 0L)

  expect_error(
    joinme_standata(
      formulaLong = y ~ 1 + time + x1 +
        (1 + time | id) +
        (0 + x1 + (1 + time | id) | marker),
      dataLong = sim$dataLong,
      formulaEvent = survival::Surv(time1, time2, status, type = "interval") ~ x1 + x2,
      dataEvent = de,
      assoc = c("cv_total")
    ),
    regexp = "type='interval'|type = \"interval\"|not supported",
    fixed = FALSE
  )
})

test_that("left truncation and time-split run end-to-end", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  sim <- simulate_joinme(
    n_id = 6,
    seed = 321,
    formulaEvent = survival::Surv(time_start, time_stop, event) ~ x1 + x2,
    covariate_formulas = list(
      x1 ~ time_varyring(rnorm, c(2, 5), mean = 0, sd = 1),
      x2 ~ rnorm(n_id)
    ),
    left_truncation_max = 3,
    truncate_longitudinal_before_entry = TRUE
  )

  fit <- joinme(
    formulaLong = y ~ 1 + time + x1 +
      (1 + time | id) +
      (0 + x1 + (1 + time | id) | marker),
    dataLong = sim$dataLong,
    formulaEvent = survival::Surv(time_start, time_stop, event) ~ x1 + x2,
    dataEvent = sim$dataEvent,
    assoc = c("cv_total"),
    control = list(
      chains = 1,
      parallel_chains = 1,
      threads_per_chain = 1,
      iter_warmup = 1000,
      iter_sampling = 1000,
      adapt_delta = 0.8,
      refresh = 0,
      seed = 321
    )
  )

  expect_s3_class(fit, "JoiNMeFit")
  expect_true(inherits(summary(fit), "SummaryJoinMeFit"))

  dr <- draws(fit)
  expect_true(inherits(dr, c("draws_array", "draws_matrix", "draws_df")))

  ex <- extract(fit, what = "assoc")
  expect_true(is.list(ex) || is.matrix(ex))

  p_fit <- plot(fit, type = c("longitudinal", "survival"), subject=2)
  expect_s3_class(p_fit, "ggplot")

  pred <- predict(
    fit,
    newdataLong = sim$dataLong,
    newdataEvent = sim$dataEvent,
    process = c("longitudinal", "event"),
    time_start = 0,
    control = list(
      chains = 1,
      parallel_chains = 1,
      threads_per_chain = 1,
      iter_warmup = 200,
      iter_sampling = 200,
      n_samples = 50,
      n_pred_draws = 200,
      refresh = 0,
      seed = 321
    )
  )

  expect_s3_class(pred, "JoiNMeDynPred")
  p_pred <- plot(pred, type = "survival")
  expect_s3_class(p_pred, "ggplot")
})

test_that("longitudinal rows can be truncated before delayed entry", {
  sim <- simulate_joinme(
    n_id = 10,
    seed = 1201,
    formulaEvent = survival::Surv(time_start, time_stop, event) ~ x1 + x2,
    left_truncation_max = 2,
    truncate_longitudinal_before_entry = TRUE
  )

  entry_by_id <- tapply(sim$dataEvent$time_start, sim$dataEvent$id, min)
  long_entry <- as.numeric(entry_by_id[as.character(sim$dataLong$id)])
  expect_true(all(sim$dataLong$time + 1e-10 >= long_entry))
})
