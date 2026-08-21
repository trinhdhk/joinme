test_that("simulated event endpoint ties receive a negligible positive separation", {
  interval_start <- c(0, 2, 3) # lower endpoints containing two exact ties and one genuinely reversed interval
  interval_stop <- c(0, 2, 2.5) # upper endpoints before numerical tie correction

  separated_stop <- .sim_separate_tied_event_endpoints(
    interval_start,
    interval_stop
  ) # corrected endpoints used by the simulation event data

  expect_equal(separated_stop[1], 1e-9)
  expect_equal(separated_stop[2], 2 + 1e-9)
  expect_equal(separated_stop[3], 2.5)
})

test_that("a boundary event time remains valid in the returned fitting data", {
  longitudinal_formula <- y ~ 1 + time +
    (1 + time | id) +
    (0 + (1 | id) | marker) # nested longitudinal structure used by simulation and standata
  event_formula <- survival::Surv(time, event) ~ 1 # ordinary subject-level right-censoring representation
  simulation <- simulate_joinme(
    formulaLong = longitudinal_formula,
    formulaEvent = event_formula,
    n_id = 3,
    families = rep("gaussian", 2),
    times_obs = 0,
    time_cens = 0,
    assoc = "cv_mean",
    seed = 4408,
    use_mirai = FALSE
  ) # zero administrative horizon forces the numerical boundary case deterministically

  standata <- joinme_standata(
    formulaLong = longitudinal_formula,
    dataLong = simulation$dataLong,
    formulaEvent = event_formula,
    dataEvent = simulation$dataEvent,
    assoc = "cv_mean",
    families = rep("gaussian", 2)
  ) # fitting representation reconstructed directly from the returned simulation

  expect_true(all(simulation$dataEvent$time > 0))
  expect_equal(simulation$dataEvent$time, simulation$dataEvent$time_stop)
  expect_true(all(standata$S_event > standata$S_entry))
})
