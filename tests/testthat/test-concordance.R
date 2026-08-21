test_that("survival-curve concordance has the Antolini orientation", {
  outcomes <- data.frame(
    residual_time = c(1, 2, 3),
    cause_event = c(1L, 1L, 0L)
  )
  ordered_survival <- matrix(
    c(
      0.2, 0.6, 0.9,
      NA, 0.3, 0.8
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("1", "2"), c("1", "2", "3"))
  )

  ordered <- .concordance_from_survival_curves(
    outcomes,
    ordered_survival
  )
  reversed <- .concordance_from_survival_curves(
    outcomes,
    1 - ordered_survival
  )

  expect_equal(ordered$concordance, 1)
  expect_equal(reversed$concordance, 0)
  expect_equal(ordered$n_pairs, 3)
  expect_equal(ordered$n_events, 2)
  expect_equal(ordered$n_subjects, 3)
})

test_that("survival-curve concordance gives half credit to prediction ties", {
  outcomes <- data.frame(
    residual_time = c(1, 2, 3),
    cause_event = c(1L, 1L, 0L)
  )
  survival <- matrix(
    c(
      0.4, 0.4, 0.8,
      NA, 0.2, 0.7
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("1", "2"), c("1", "2", "3"))
  )

  out <- .concordance_from_survival_curves(outcomes, survival)

  expect_equal(out$concordant, 2)
  expect_equal(out$tied, 1)
  expect_equal(out$n_pairs, 3)
  expect_equal(out$concordance, 2.5 / 3)
})

test_that("premature censoring is incomparable and same-time censoring is comparable", {
  outcomes <- data.frame(
    residual_time = c(0.5, 1, 1, 2),
    cause_event = c(0L, 0L, 1L, 0L)
  )
  survival <- matrix(
    c(0.9, 0.7, 0.2, 0.8),
    nrow = 1,
    dimnames = list("1", as.character(seq_len(4)))
  )

  out <- .concordance_from_survival_curves(outcomes, survival)

  # Subject 2 was lost before the event and is not comparable.  Subject 3 was
  # censored at the event time and is known to have survived through that time.
  expect_equal(out$n_pairs, 2)
  expect_equal(out$concordance, 1)
})

test_that("concordance event weights follow survival::concordance", {
  outcomes <- data.frame(
    residual_time = c(1, 2, 3, 4, 5),
    cause_event = c(1L, 1L, 0L, 1L, 0L)
  )

  harrell <- .concordance_event_weights(outcomes, "none")
  uno <- .concordance_event_weights(outcomes, "n/G2")

  expect_equal(harrell[outcomes$cause_event == 1L], rep(1, 3))
  expect_true(all(uno[outcomes$cause_event == 1L] >= harrell[outcomes$cause_event == 1L]))
  expect_error(
    .concordance_event_weights(outcomes, "unsupported"),
    "type_weights"
  )
})

test_that("counting-process event rows reduce to one terminal outcome", {
  data_event <- data.frame(
    id = c(2, 1, 2, 1),
    start = c(0, 1, 1, 0),
    stop = c(1, 2, 3, 1),
    event = c(0L, 1L, 0L, 0L)
  )

  outcome <- .subject_event_outcomes(
    survival::Surv(start, stop, event) ~ 1,
    data_event,
    id_var = "id",
    context = "test"
  )

  expect_equal(outcome$id, c(1, 2))
  expect_equal(outcome$event_time, c(2, 3))
  expect_equal(outcome$event_status, c(1L, 0L))
})

test_that("subject-specific time columns remain aligned across interval rows", {
  data_event <- data.frame(
    id = c(1, 1, 2, 2),
    landmark = c(0.5, 0.5, 1, 1)
  )

  mapped <- .subject_time_map(
    "landmark",
    ids = c(1, 2),
    data_event = data_event,
    id_var = "id",
    argument = "time_start"
  )

  expect_equal(unname(mapped), c(0.5, 1))
  data_event$landmark[2] <- 0.75
  expect_error(
    .subject_time_map(
      "landmark",
      ids = c(1, 2),
      data_event = data_event,
      id_var = "id",
      argument = "time_start"
    ),
    "constant value"
  )
})

test_that("last measurement landmarks strictly precede observed outcomes", {
  data_long <- data.frame(
    id = c(1, 1, 1, 2, 2),
    time = c(0, 1, 2, 0, 3)
  )

  landmark <- .last_preoutcome_measurement(
    ids = c(1, 2),
    event_time = c(2, 4),
    data_long = data_long,
    id_var = "id",
    time_var = "time"
  )

  expect_equal(unname(landmark), c(1, 3))
})

test_that("concordance prediction grids contain every observable event time", {
  object <- structure(list(
    call = quote(joinme(id_var = "id", time_var = "time")),
    formulaEvent = survival::Surv(time, event) ~ 1,
    config = list(),
    stan_data = list()
  ), class = "JoiNMeFit")
  data_long <- data.frame(
    id = c(1, 1, 1, 2, 2, 3, 3),
    time = c(0, 1, 2, 0, 1.5, 0, 1),
    marker = "m1",
    y = 0
  )
  data_event <- data.frame(
    id = 1:3,
    time = c(2, 3, 3),
    event = c(1L, 0L, 1L)
  )
  testthat::local_mocked_bindings(
    predict.JoiNMeFit = function(
        object, newdataLong, newdataEvent, process, times, time_start,
        control, seed, ...) {
      expect_equal(unname(time_start), c(1, 1.5, 1))
      expect_true(2.5 %in% times[["2"]])
      expect_true(all(c(2, 3) %in% times[["3"]]))
      expect_lte(max(newdataLong$time[newdataLong$id == 1]), 1)

      rows <- lapply(names(times), function(id) {
        residual <- times[[id]] - time_start[[id]]
        data.frame(
          id = as.integer(id),
          time = times[[id]],
          Survival = exp(-as.integer(id) * residual / 10)
        )
      })
      structure(list(
        predictions = list(survival = do.call(rbind, rows))
      ), class = "JoiNMeDynPred")
    },
    .package = "joinme"
  )

  curves <- .concordance_survival_curves(
    object = object,
    newdataLong = data_long,
    newdataEvent = data_event,
    time_start = NULL,
    cause = 1L,
    predict_control = list(n_samples = 5, n_pred_draws = 5),
    seed = 1
  )

  expect_equal(curves$outcomes$residual_time, c(1, 1.5, 2))
  expect_equal(as.numeric(rownames(curves$survival)), c(1, 2))
  expect_equal(dim(curves$survival), c(2, 3))
})

test_that("concordance has no horizon and directs horizon analysis to tvROC", {
  object <- structure(list(
    dataLong = data.frame(id = 1:3),
    dataEvent = data.frame(id = 1:3),
    config = list()
  ), class = "JoiNMeFit")
  testthat::local_mocked_bindings(
    .concordance_survival_curves = function(
        object, newdataLong, newdataEvent, time_start, cause,
        predict_control, seed, ...) {
      expect_null(time_start)
      expect_equal(predict_control$n_samples, 200)
      list(
        outcomes = data.frame(
          residual_time = c(1, 2, 3),
          cause_event = c(1L, 1L, 0L)
        ),
        survival = matrix(
          c(0.2, 0.7, 0.9, NA, 0.3, 0.8),
          nrow = 2,
          byrow = TRUE,
          dimnames = list(c("1", "2"), as.character(1:3))
        )
      )
    },
    .package = "joinme"
  )

  out <- suppressWarnings(concordance(object))

  expect_s3_class(out, "concordance_JoiNMeFit")
  expect_equal(out$concordance, 1)
  expect_error(
    suppressWarnings(concordance(object, time_horizon = 3)),
    "tvROC\\(object"
  )
})

test_that("tvROC and tvAUC retain the established horizon parameter recipe", {
  object <- structure(list(
    dataLong = data.frame(id = 1:2),
    dataEvent = data.frame(id = 1:2),
    tmax = 3,
    config = list(),
    stan_data = list(tmax = 3)
  ), class = "JoiNMeFit")
  testthat::local_mocked_bindings(
    .dynamic_discrimination_risk_set = function(
        object, newdataLong, newdataEvent, time_start, time_horizon,
        cause, n_samples, seed, ...) {
      expect_equal(time_start, 0)
      expect_equal(time_horizon, 3)
      data.frame(
        risk = c(0.8, 0.2),
        event_window = c(1L, 0L),
        event_time = c(2, 4),
        event_status = c(1L, 0L),
        event_type = c(1L, 0L),
        time_horizon = c(3, 3)
      )
    },
    .package = "joinme"
  )

  roc <- suppressWarnings(tvROC(object))
  out <- suppressWarnings(tvAUC(object))

  expect_s3_class(roc, "tvROC")
  expect_s3_class(out, "tvAUC")
  expect_equal(roc$Tstart, 0)
  expect_equal(roc$Thoriz, 3)
  expect_equal(out$Tstart, 0)
  expect_equal(out$Thoriz, 3)
  expect_equal(out$auc, 1)
})
