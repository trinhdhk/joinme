devtools::load_all()
library(dplyr)

sim <- simulate_joinme(
  n_id = 120,
  # formulaDist = list(
  #   sigma ~ 1 + (1 | marker)
  # ),
  families = c(rep("gaussian", 5)),
  n_obs_per_marker_per_id = 6,
  times_obs = seq(0, 6, length.out = 10),
  beta_long = c(0.5, -1, 1),
  seed = 112,
  marker_weights = c(0.5, 1, -0.5, 0.75, 1),
  assoc = c("cv_mean"),
  # transforms = list(cv_mean = list(type = "functional", expr =  ~ expit(x))),
  assoc_coefs = c(cv_mean = 0.25)
)

formulaLong <- y ~ 1 +
  time +
  x1 +
  (1 + time | id) +
  (0 + x1 + (1 + time | id) | marker)
formulaEvent <- survival::Surv(time, event) ~ x1 + x2

fit <- joinme(
  formulaLong = formulaLong,
  dataLong = sim$dataLong,
  formulaEvent = formulaEvent,
  # formulaDist = list(
  #   sigma ~ 1 + (1 | marker)
  #  ),
  dataEvent = sim$dataEvent,
  assoc = c("cv_mean"),
  families = c(rep("gaussian", 5)),
  # transforms = list(cv_mean = list(type = "functional", expr =  ~ expit(x))),
  # estimate_marker_weights = FALSE,
  # marker_weights = sim$truth$marker_weights,
  estimate_marker_weights = TRUE,
  # marker_weight_scale = 1,
  control = list(
    threads_per_chain = 6,
    parallel_chains = 2,
    iter_warmup = 300,
    iter_sampling = 1200,
    refresh = 200,
    adapt_delta = 0.85,
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
  newdataLong = sim$dataLong |> filter(id%in% c(5, 6, 7, 8)),
  newdataEvent = ndE |> filter(id%in% c(5, 6, 7, 8)),
  time_start = "time_start",
  # times = seq(0, max(sim$dataLong$time) + 1, length.out = 20),
  time_horizon = 10,
  control = list(
    n_samples = 50,
    n_times = 50,
    chains = 2, parallel_chains =2,
    iter_warmup = 200,
    iter_sampling = 800,
    show_messages = TRUE,
    refresh=50,
    threads_per_chain = 6
  )
)

p <- plot(pred, which = c("longitudinal", "survival"), combined = TRUE)

if (!dir.exists("tests/local")) dir.create("tests/local", recursive = TRUE)
p3 <- if (is.list(p) && "3" %in% names(p)) p[["3"]] else p
ggplot2::ggsave(
  filename = "tests/local/subject3_after_fix.png",
  plot = p3,
  width = 12,
  height = 6,
  dpi = 150
)
