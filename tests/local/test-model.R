devtools::load_all()
library(dplyr)

sim <- simulate_joinme(
  formulaLong = y ~ 1 + time + x1 + (1 + time | id) + (0 + x1 + (1 + time | id) | marker),
  formulaEvent = survival::Surv(time, event) ~ x2,
  formulaVCov = ~ x1,
  families = c("gaussian", "student_t", "student_t"),
  n_id = 100,
  times_obs = seq(0, 6, length.out = 10),
  assoc = c("cv_total", "corr"),
  truth = jm_truth(
    longitudinal = c('(Intercept)' = 0.5, time = 0.3, x1 = -0.2),
    survival = list(slope = c(x2 = 0.4)),
    assoc_coef = list(slope = c(cv_total = 0.35, "corr[1]" = 0.1)),
    vcov = list(
      sd = list(intercept = c(-0.2, -0.1), latent = 0.175),
      corr = list(intercept = 0, latent = 0.175)
    ),
    marker_weights = list(offset = c(1, -0.5, 2), family = "constant"),
    re_params = list(
      id = list(sd = c(0.6, 0.3)),
      marker = list(sd = 0.25),
      dist = list(sigma = list(sd = 0.4))
    )
  ),
  transforms = joinme_tf(corr = ~ -x),
  seed = 2026,
  use_mirai = TRUE,
  n_workers = 5
)

fit <- joinme(
  formulaLong = y ~ 1 + time + x1 + (1 + time | id) + (0 + x1 + (1 + time | id) | marker),
  formulaEvent = survival::Surv(time, event) ~ x2,
  formulaVCov = ~ x1,
  families = c("gaussian", "student_t", "student_t"),
  dataLong = sim$dataLong,
  dataEvent = sim$dataEvent,
  assoc = c("cv_total", "corr"),
  transforms = joinme_tf(corr = ~ -x),
  # transforms = list(cv_mean = list(type = "functional", expr =  ~ expit(x))),
  # priors = jm_prior(marker_weights = list(
  #   offset = sim$truth$marker_weights,
  #   family = "constant"
  # )),
  control = list(
    threads_per_chain = 6,
    parallel_chains = 2,
    iter_warmup = 1000,
    iter_sampling = 500,
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
    iter_warmup = 400,
    iter_sampling = 800,
    show_messages = TRUE,
    refresh=50,
    threads_per_chain = 6
  )
)

p <- plot(pred, type = c("longitudinal", "survival"), combined = TRUE)

p3 <- if (is.list(p) && "3" %in% names(p)) p[["3"]] else p
ggplot2::ggsave(
  filename = "tests/local/subject3.png",
  plot = p3,
  width = 12,
  height = 6,
  dpi = 150
)
