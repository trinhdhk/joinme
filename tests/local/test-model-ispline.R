devtools::load_all()
library(dplyr)

expit_grid <- stats::plogis(seq(-4, 4, length.out = 31))

# Simulation uses the plug-in penalised spline constructor, which now matches
# the Stan-estimated fit convention:
# - 6 knots including boundaries
# - degree 2
# - 7 spline coefficients (= length(knots) + degree - 1)
# - first/last coefficients anchored at 0 and 1
# sim_cv_tf <- penalised_ispline_transform(
#   x = c(0, 0.25, 0.5, 0.6, 0.8, 1),
#   y = c(0, 0.5, 0.6, 0.75, 0.85, 1),
#   knots = seq(0, 1, 0.2),
#   degree = 2,
#   lambda = 1
# )



# Fitting can use the Stan-estimated penalised spline path by omitting y.
fit_cv_tf <- list(
  type = "ispline_penalised",
  knots = seq(-4, 4, 0.5),
  degree = 2,
  lambda = 1
)

sim <- simulate_joinme(
  formulaLong = y ~ 1 + time + x1 + (1 + time | id) + (0 + x1 + (1 + time | id) | marker),
  formulaEvent = survival::Surv(time, event) ~ x2,
  formulaCorr = ~ x1,
  families = c("gaussian", "gaussian"),
  n_id = 200,
  n_obs_per_marker_per_id = 6,
  times_obs = seq(0, 6, length.out = 8),
  assoc = c("cv_total", "corr"),
  beta_long = c('(Intercept)' = 0.5, time = 0.3, x1 = -0.2),
  beta_event = c(x2 = 0.2),
  assoc_coefs = list(cv_total = 0.25, corr = c(0.12)),
  transforms = joinme_tf(
    # cv_total = ~ expit(x), #sim_cv_tf,
    corr = list(
      type = "ispline_expit_penalised",
      x = expit_grid,
      y = stats::plogis(-3 * qlogis(expit_grid)),
      n_knots = 6,
      degree = 3,
      lambda = 1
    )
  ),
  re_params = list(
    id = list(sd = c(0.5, 0.25)),
    marker = list(sd = 0.3),
    id_marker_cov = list(
      latent = list(sd = c(0.6, 0.25)),
      alpha = c(-0.1, 0.05, -0.05),
      lambda = 0.3,
      sd_u = 0.4
    )
  ),
  use_mirai = TRUE,
  n_workers = 8,
  seed = 2026
)

engine <- if (requireNamespace("cmdstanr", quietly = TRUE)) "cmdstanr" else "rstan"
control <- list(
  engine = engine,
  chains = 2,
  iter_warmup = 1000,
  iter_sampling = 1000,
  parallel_chains = 2,
  threads_per_chain = 6,
  adapt_delta = 0.78,
  max_treedepth = 12,
  seed = 421
)


fit <- joinme(
  formulaLong = y ~ 1 + time + x1 + (1 + time | id) + (0 + x1 + (1 + time | id) | marker),
  formulaEvent = survival::Surv(time, event) ~ x2,
  formulaCorr = ~ x1,
  dataLong = sim$dataLong,
  dataEvent = sim$dataEvent,
  assoc = c("corr"),
  transforms = joinme_tf(
    # cv_total = fit_cv_tf,
    # cv_total = ~ expit(x),
    corr = list(
      type = "ispline_expit_penalised",
      x = expit_grid,
      degree = 3,
      lambda = 1
    )
  ),
  control = control
)

plot(fit, type = c("association"), show_data = TRUE)
sum_obj <- summary(fit)

lastTime <- sim$dataLong %>%
  group_by(id) %>%
  summarize(last_time = max(time)) %>%
  ungroup()

ndE <- sim$dataEvent %>%
  left_join(lastTime, by = "id") %>%
  mutate(time_start = as.numeric(last_time))

pred <- predict(
  fit,
  newdataLong = sim$dataLong |> filter(id %in% c(5,6,7,8)),
  newdataEvent = ndE |> filter(id %in% c(5,6,7,8)),
  process = c("longitudinal", "event"),
  time_start = "time_start",
  time_horizon = 10,
  control = list(
    n_samples = 100,
    n_times = 50,
    chains = 2, parallel_chains = 2,
    iter_warmup = 500,
    iter_sampling = 500,
    show_messages = TRUE,
    refresh = 50,
    threads_per_chain = 6
  )
)

p <- plot(pred, type = c("longitudinal", "survival"), combined = TRUE)

p3 <- if (is.list(p) && "3" %in% names(p)) p[["3"]] else p
ggplot2::ggsave(
  filename = "tests/local/subject3-ispline.png",
  plot = p3,
  width = 12,
  height = 6,
  dpi = 150
)
