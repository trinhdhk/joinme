devtools::load_all()
library(ggplot2)

set.seed(42)
sim <- simulate_joinme(
  n_id = 40,
  families = c(rep("student_t", 3), "gaussian"),
  n_obs_per_marker_per_id = 6,
  times_obs = seq(0, 6, length.out = 10),
  beta_long = c(0.5, -2, 1),
  seed = 42,
  assoc = c("cv_total"),
  assoc_coefs = c(cv_total = 0.6)
)

formulaLong <- y ~ 1 + time + x1 +
  (1 + time | id) +
  (0 + x1 + (1 + time | id) | marker)
formulaEvent <- survival::Surv(time, event) ~ 1 + x1 + x2

fit <- joinme(
  formulaLong = formulaLong,
  dataLong = sim$dataLong,
  formulaEvent = formulaEvent,
  dataEvent = sim$dataEvent,
  assoc = c("cv_total"),
  families = c(rep("student_t", 3), "gaussian"),
  transforms = list(cv_total = list(type = "identity")),
  fixed_marker_weights = FALSE,
  marker_weight_scale = 1,
  control = list(
    threads_per_chain = 2,
    parallel_chains = 2,
    iter_warmup = 80,
    iter_sampling = 160,
    refresh = 0,
    adapt_delta = 0.9,
    max_treedepth = 12,
    seed = 421
  )
)

last_time <- aggregate(time ~ id, data = sim$dataLong, FUN = max)
names(last_time)[2] <- "time_start"
ndE <- merge(sim$dataEvent, last_time, by = "id", all.x = TRUE, sort = FALSE)

ndL3 <- sim$dataLong[sim$dataLong$id %in% c(1, 2, 3), , drop = FALSE]
ndE3 <- ndE[ndE$id %in% c(1, 2, 3), , drop = FALSE]

pred <- posterior_predict(
  fit,
  newdataLong = ndL3,
  newdataEvent = ndE3,
  time_start = "time_start",
  control = list(
    n_samples = 50,
    n_times = 50,
    chains = 2,
    parallel_chains = 2,
    iter_warmup = 80,
    iter_sampling = 160,
    show_messages = FALSE,
    refresh = 0,
    threads_per_chain = 2
  )
)

p <- plot(pred, type = c("longitudinal", "survival"), combined = TRUE)
p3 <- if (is.list(p) && "3" %in% names(p)) p[["3"]] else p

out <- "tests/local/subject3_after_fix_quick.png"
ggplot2::ggsave(filename = out, plot = p3, width = 12, height = 6, dpi = 150)

s3 <- pred$predictions$survival[pred$predictions$survival$id == 3, , drop = FALSE]
cat("Saved:", out, "\n")
cat("Subject3 survival range:", min(s3$Survival, na.rm = TRUE), max(s3$Survival, na.rm = TRUE), "\n")
