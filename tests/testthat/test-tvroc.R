test_that("tvROC and tvAUC are re-exported JMbayes2 generics", {
  expect_identical(tvROC, JMbayes2::tvROC)
  expect_identical(tvAUC, JMbayes2::tvAUC)
  expect_true("tvROC.JoiNMeFit" %in% methods("tvROC"))
  expect_true("tvAUC.JoiNMeFit" %in% methods("tvAUC"))
})

test_that("tvROC time resolution preserves the former AUC parameter recipe", {
  object <- structure(
    list(tmax = 7, stan_data = list(tmax = 7)),
    class = "JoiNMeFit"
  )
  data_event <- data.frame(time = c(2, 7), event = c(1L, 0L))

  by_width <- .get_tvroc_times(
    object,
    data_event,
    time_start = c(1, 2),
    time_horizon = NULL,
    Dt = 3
  )
  by_horizon <- .get_tvroc_times(
    object,
    data_event,
    time_start = c(1, 2),
    time_horizon = 6,
    Dt = NULL
  )
  by_default <- .get_tvroc_times(
    object,
    data_event,
    time_start = 1,
    time_horizon = NULL,
    Dt = NULL
  )

  expect_equal(by_width$time_horizon, c(4, 5))
  expect_equal(by_horizon$time_horizon, c(6, 6))
  expect_equal(by_default$time_horizon, 7)
  expect_error(
    .get_tvroc_times(
      object,
      data_event,
      time_start = 2,
      time_horizon = 2,
      Dt = NULL
    ),
    "greater than"
  )
})

test_that("model-based ROC weights represent uncertain horizon status", {
  risk_set <- data.frame(
    risk = c(0.9, 0.2, 0.3, 0.4, 0.6),
    event_window = c(1L, 0L, 0L, 0L, 0L),
    event_time = c(1, 4, 1.5, 4, 1.7),
    event_status = c(1L, 1L, 0L, 0L, 1L),
    event_type = c(1L, 1L, 0L, 0L, 2L),
    time_horizon = rep(3, 5)
  )

  weights <- .tvroc_status_weights(
    risk_set,
    type_weights = "model-based",
    cause = 1L
  )

  expect_equal(weights$case, c(1, 0, 0.3, 0, 0.6))
  expect_equal(weights$control, 1 - weights$case)
})

test_that("IPCW ROC weights use observable cases and controls", {
  risk_set <- data.frame(
    risk = c(0.9, 0.2, 0.3, 0.4, 0.6, 0.1),
    event_window = c(1L, 0L, 0L, 0L, 0L, 0L),
    event_time = c(1, 4, 1.5, 5, 1.7, 6),
    event_status = c(1L, 1L, 0L, 0L, 1L, 1L),
    event_type = c(1L, 1L, 0L, 0L, 2L, 2L),
    time_horizon = rep(3, 6)
  )

  weights <- .tvroc_status_weights(
    risk_set,
    type_weights = "IPCW",
    cause = 1L
  )

  expect_gt(weights$case[[1]], 0)
  expect_equal(weights$case[-1], rep(0, 5))
  expect_true(all(weights$control[c(2, 4, 6)] > 0))
  expect_equal(weights$control[c(1, 3, 5)], rep(0, 3))
})

test_that("one ROC curve follows the JMbayes2 tvROC contract", {
  risk_set <- data.frame(
    risk = c(0.9, 0.7, 0.3, 0.1),
    event_window = c(1L, 1L, 0L, 0L),
    event_time = c(1.2, 1.7, 4, 5),
    event_status = c(1L, 1L, 0L, 0L),
    event_type = c(1L, 1L, 0L, 0L),
    time_horizon = rep(3, 4)
  )
  attr(risk_set, "risk_draws") <- cbind(
    c(0.92, 0.72, 0.32, 0.12),
    c(0.88, 0.68, 0.28, 0.08)
  )

  roc <- .tvroc_from_risk_set(
    risk_set,
    time_start = 0,
    time_horizon = 3,
    type_weights = "model-based",
    cause = 1L,
    object_name = "fit"
  )
  area <- JMbayes2::tvAUC(roc)

  expect_s3_class(roc, "tvROC")
  expect_s3_class(area, "tvAUC")
  expect_named(
    roc,
    c(
      "TP", "FP", "nTP", "nFN", "nFP", "nTN", "tp", "fp", "thrs",
      "thr", "F1score", "Youden", "Tstart", "Thoriz", "nr",
      "classObject", "type_weights", "nameObject", "cause"
    )
  )
  expect_equal(dim(roc$tp), c(101, 2))
  expect_equal(dim(roc$fp), c(101, 2))
  expect_true(all(diff(roc$TP) >= -sqrt(.Machine$double.eps)))
  expect_true(all(diff(roc$FP) >= -sqrt(.Machine$double.eps)))
  expect_equal(area$auc, 1)

  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  expect_invisible(plot(roc))
  grDevices::dev.off()
})

test_that("posterior event-risk draws are aligned at each exact horizon", {
  prediction <- list(draws = list(survival = list(
    "1" = list(
      matrix = rbind(c(1, 0.8, 0.6), c(1, 0.7, 0.5)),
      time = c(0, 1, 2)
    ),
    "2" = list(
      matrix = rbind(c(1, 0.9, 0.4), c(1, 0.85, 0.3)),
      time = c(0, 1, 3)
    )
  )))

  risks <- .discrimination_horizon_risk_draws(
    prediction,
    ids = 1:2,
    horizons = c("1" = 2, "2" = 3)
  )

  expect_equal(unname(risks), rbind(c(0.4, 0.5), c(0.6, 0.7)))
})

test_that("the cumulative case window includes an event at the horizon", {
  object <- structure(
    list(
      call = quote(joinme(id_var = "id", time_var = "time")),
      formulaEvent = survival::Surv(event_time, event) ~ 1
    ),
    class = "JoiNMeFit"
  )
  data_long <- data.frame(id = 1:2, time = c(0, 0), y = c(0, 0))
  data_event <- data.frame(
    id = 1:2,
    event_time = c(3, 4),
    event = c(1L, 0L)
  )
  testthat::local_mocked_bindings(
    predict.JoiNMeFit = function(
        object, newdataLong, newdataEvent, process, times, time_start,
        control, seed, ...) {
      prediction_rows <- data.frame(
        id = rep(1:2, each = 50),
        time = c(times[["1"]], times[["2"]]),
        Survival = c(
          seq(1, 0.2, length.out = 50),
          seq(1, 0.8, length.out = 50)
        )
      )
      list(
        predictions = list(survival = prediction_rows),
        draws = list(survival = list(
          "1" = list(
            matrix = rbind(
              seq(1, 0.25, length.out = 50),
              seq(1, 0.15, length.out = 50)
            ),
            time = times[["1"]]
          ),
          "2" = list(
            matrix = rbind(
              seq(1, 0.85, length.out = 50),
              seq(1, 0.75, length.out = 50)
            ),
            time = times[["2"]]
          )
        ))
      )
    },
    .package = "joinme"
  )

  risk_set <- .dynamic_discrimination_risk_set(
    object,
    newdataLong = data_long,
    newdataEvent = data_event,
    time_start = 0,
    time_horizon = 3,
    cause = 1L,
    n_samples = 2L,
    seed = 1
  )

  expect_equal(risk_set$event_window, c(1L, 0L))
  expect_equal(risk_set$risk, c(0.8, 0.2))
  expect_equal(dim(attr(risk_set, "risk_draws")), c(2, 2))
})

test_that("multiple landmark curves and areas retain their alignment", {
  object <- structure(list(), class = "JoiNMeFit")
  testthat::local_mocked_bindings(
    .get_train_data = function(object, newdataLong, newdataEvent, purpose) {
      list(
        newdataLong = data.frame(id = 1:4),
        newdataEvent = data.frame(id = 1:4)
      )
    },
    .dynamic_discrimination_risk_set = function(
        object, newdataLong, newdataEvent, time_start, time_horizon,
        cause, n_samples, seed, ...) {
      data.frame(
        risk = c(0.9, 0.7, 0.3, 0.1),
        event_window = c(1L, 1L, 0L, 0L),
        event_time = c(time_start + 0.5, time_start + 1, 5, 6),
        event_status = c(1L, 1L, 0L, 0L),
        event_type = c(1L, 1L, 0L, 0L),
        time_horizon = rep(time_horizon, 4)
      )
    },
    .package = "joinme"
  )

  roc <- suppressWarnings(tvROC(object, time_start = c(0, 1), Dt = 2))
  area <- suppressWarnings(tvAUC(object, time_start = c(0, 1), Dt = 2))

  expect_s3_class(roc, "tvROC_JoiNMeFit_list")
  expect_length(roc$curves, 2)
  expect_equal(roc$Tstart, c(0, 1))
  expect_equal(roc$Thoriz, c(2, 3))
  expect_s3_class(area, "tvAUC_JoiNMeFit")
  expect_equal(area$time_start, c(0, 1))
  expect_equal(area$time_horizon, c(2, 3))
  expect_equal(area$auc, c(1, 1))
})

test_that("ROC estimation rejects a window without requested-cause events", {
  risk_set <- data.frame(
    risk = c(0.4, 0.2),
    event_window = c(0L, 0L),
    event_time = c(4, 5),
    event_status = c(0L, 0L),
    event_type = c(0L, 0L),
    time_horizon = c(3, 3)
  )

  expect_error(
    .tvroc_from_risk_set(
      risk_set,
      time_start = 0,
      time_horizon = 3,
      type_weights = "model-based",
      cause = 1L,
      object_name = "fit"
    ),
    "No requested-cause event"
  )
})
