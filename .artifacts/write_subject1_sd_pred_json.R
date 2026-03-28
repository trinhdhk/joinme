devtools::load_all(quiet = TRUE)

env <- new.env(parent = emptyenv())
load("tmp_fit-sim3.Rdata", envir = env)
fit <- env[["fit"]]
sim <- env[["sim"]]

meta <- joinme:::.recover_metadata(fit, NULL)
forms <- joinme:::.parse_formulas(fit)
seed <- 2401
n_samples <- 50L
n_pred_draws <- 50L
id_var <- eval(fit$call$id_var) %||% "id"
time_var <- eval(fit$call$time_var) %||% "time"

draws_list_raw <- joinme:::.extract_draws_for_pred(fit, n_samples, seed)
n_samples_extracted <- joinme:::.n_draws_in_prediction_list(draws_list_raw)
pred_draw_index <- joinme:::.prediction_draw_index(n_samples_extracted, n_pred_draws, seed)
draws_list <- joinme:::.subset_draws_for_prediction(draws_list_raw, pred_draw_index)

id <- unique(sim$dataEvent[[id_var]])[1]
dE <- sim$dataEvent[sim$dataEvent[[id_var]] == id, , drop = FALSE]
dL <- sim$dataLong[sim$dataLong[[id_var]] == id, , drop = FALSE]
sd_pred <- joinme:::.prepare_subject_standata(
  dE, dL, fit, meta$tmax, meta$knots, meta$col_means,
  1, numeric(0), seq(1, 3, length.out = 50),
  forms, draws_list,
  grainsize = NULL,
  control = list(quadrature_nodes = fit$stan_data$quadrature_nodes %||% 15L),
  degree = meta$degree
)
sd_pred <- joinme:::.coerce_rstan_time_indices(sd_pred)
sd_pred <- joinme:::.coerce_rstan_dist_arrays(sd_pred)
sd_pred <- joinme:::.coerce_rstan_vectors(sd_pred, c(
  "vec_cov_vcov",
  "const_data_cv",
  "const_data_cs",
  "const_data_corr",
  "const_data_vcov",
  "const_data_cv_mean",
  "const_data_cv_marker",
  "const_data_cs_mean",
  "const_data_cs_marker"
))
cmdstanr::write_stan_json(sd_pred, ".artifacts/subject1_sd_pred.json")
cat("WROTE\n")
