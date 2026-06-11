# End-to-end smoke for family-specific link handling.
# Steps:
# 1) simulate mixed-family data,
# 2) fit with per-family custom links,
# 3) predict with all longitudinal scales,
# 4) save plots,
# 5) run simple image-quality anomaly checks.

devtools::load_all(quiet = TRUE)

set.seed(2602)
sim <- simulate_joinme(
  n_id = 30,
  families = c("gaussian", "bernoulli", "poisson"),
  n_obs_per_marker_per_id = 4,
  times_obs = seq(0, 5, length.out = 8),
  assoc = c("cv_total"),
  assoc_coefs = c(cv_total = 0.2),
  seed = 2602
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
  families = list(
    jm_family("gaussian", link = "identity"),
    jm_family("bernoulli", link = "probit"),
    jm_family("poisson", inv_link = ~ exp(x))
  ),
  control = list(
    engine = "cmdstanr",
    chains = 1,
    parallel_chains = 1,
    iter_warmup = 80,
    iter_sampling = 80,
    refresh = 0,
    seed = 2602
  )
)

# Build landmark-time input for prediction
last_time <- aggregate(time ~ id, data = sim$dataLong, FUN = max)
names(last_time)[2] <- "time_start"
ndE <- merge(sim$dataEvent, last_time, by = "id", all.x = TRUE, sort = FALSE)
ndL <- sim$dataLong[sim$dataLong$id %in% c(1, 2, 3), , drop = FALSE]
ndE3 <- ndE[ndE$id %in% c(1, 2, 3), , drop = FALSE]

pred <- predict(
  fit,
  newdataLong = ndL,
  newdataEvent = ndE3,
  time_start = "time_start",
  scale = c("epred", "linpred", "predict"),
  control = list(
    n_samples = 40,
    n_pred_draws = 40,
    n_times = 30,
    chains = 1,
    iter_warmup = 40,
    iter_sampling = 40,
    refresh = 0
  ),
  seed = 2602
)

stopifnot(inherits(pred, "JoinMeDynPred"))
stopifnot(all(c("epred", "linpred", "predict") %in% pred$metadata$scales))

if (!dir.exists("tests/local")) dir.create("tests/local", recursive = TRUE)

p_long <- plot(pred, type = "longitudinal", scale = "epred")
p_surv <- plot(pred, type = "survival")

long_path <- "tests/local/e2e-link-longitudinal.png"
surv_path <- "tests/local/e2e-link-survival.png"

ggplot2::ggsave(long_path, plot = if (inherits(p_long, "ggplot")) p_long else p_long[[1]], width = 10, height = 6, dpi = 140)
ggplot2::ggsave(surv_path, plot = if (inherits(p_surv, "ggplot")) p_surv else p_surv[[1]], width = 10, height = 6, dpi = 140)

# Lightweight computer-vision-like anomaly checks:
# - image exists and has non-trivial size,
# - finite pixel matrix,
# - non-zero pixel variance (not blank/flat render).
if (!requireNamespace("png", quietly = TRUE)) {
  stop("Package 'png' is required for image anomaly checks.")
}

check_image <- function(path) {
  if (!file.exists(path)) stop("Missing output image: ", path)
  sz <- file.info(path)$size
  if (!is.finite(sz) || sz < 15000) stop("Image too small; possible render failure: ", path)
  arr <- png::readPNG(path)
  if (any(!is.finite(arr))) stop("Non-finite pixels detected: ", path)
  if (stats::var(as.numeric(arr)) < 1e-7) stop("Near-constant image detected (potential blank plot): ", path)
  TRUE
}

check_image(long_path)
check_image(surv_path)

cat("Saved:", long_path, "\n")
cat("Saved:", surv_path, "\n")
cat("Prediction rows (longitudinal):", nrow(pred$predictions$longitudinal), "\n")
cat("Prediction rows (survival):", nrow(pred$predictions$survival), "\n")
