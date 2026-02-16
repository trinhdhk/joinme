devtools::load_all()
library(dplyr)

sim <- simulate_joinme(
  n_id = 20,
  families = c(rep("student_t", 3), 'bernoulli'),
  n_obs_per_marker_per_id = 6,
  times_obs = seq(0, 6, length.out = 10),
  beta_long = c(0.5, -2, 1),
  seed = 42,
  assoc = c("cv_total"),
  assoc_coefs = c(cv_total = 0.6)
)

formulaLong <- y ~ 1 +
  time +
  x1 +
  (1 + time | id) +
  (0 + x1 + (1 + time | id) | marker)
formulaEvent <- survival::Surv(time, event) ~ 1 + x1 + x2

fit <- joinme(
  formulaLong = formulaLong,
  dataLong = sim$dataLong,
  formulaEvent = formulaEvent,
  dataEvent = sim$dataEvent,
  assoc = c("cv_total"),
  families = c(rep("student_t", 3), 'bernoulli'),
  transforms = list(cv_total = list(type = "identity")),
  estimate_marker_weights = TRUE,
  marker_weight_scale = 1,
  control = list(
    threads_per_chain = 5,
    parallel_chains = 2,
    iter_warmup = 200,
    iter_sampling = 1000,
    refresh = 200,
    adapt_delta = 0.75,
    max_treedepth = 12,
    seed = 421
  )
)


sum_obj <- summary(fit)

lastTime <- sim$dataLong %>%
  group_by(id) %>%
  summarize(last_time = max(time)) %>%
  ungroup()

ndE <- sim$dataEvent %>%
  left_join(lastTime, by = "id") %>%
  mutate(time_start = as.numeric(last_time))

pred <- posterior_predict(
  fit,
  newdataLong = sim$dataLong |> filter(id%in% c(1, 2)),
  newdataEvent = ndE |> filter(id%in% c(1, 2)),
  time_start = "time_start",
  # times = seq(0, max(sim$dataLong$time) + 1, length.out = 20),
  n_samples = 100,
  n_times = 50,
  control = list(
    chains = 2,
    iter_warmup = 100,
    iter_sampling = 200,
    threads_per_chain = 6
  )
)

p <- plot(pred, which = c("longitudinal", "survival"), combined = TRUE)
