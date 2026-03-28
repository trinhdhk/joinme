devtools::load_all(quiet = TRUE)

env <- new.env(parent = emptyenv())
load("tmp_fit-sim3.Rdata", envir = env)
fit <- env[["fit"]]
sim <- env[["sim"]]
id_var <- eval(fit$call$id_var) %||% "id"
ids <- unique(sim$dataEvent[[id_var]])
ids <- ids[!is.na(ids)]
time_start_map <- stats::setNames(rep(1, length(ids)), ids)
time_horizon_map <- stats::setNames(rep(3, length(ids)), ids)
times_map <- stats::setNames(lapply(ids, function(id) seq(time_start_map[[as.character(id)]], time_horizon_map[[as.character(id)]], length.out = 50)), ids)

old_force <- getOption("joinme.force_recompile")
options(joinme.force_recompile = TRUE)
on.exit(options(joinme.force_recompile = old_force), add = TRUE)

out_file <- ".artifacts/check_dense_predict_recompile_result.txt"
result <- tryCatch({
  pred <- predict(
    fit,
    newdataLong = sim$dataLong,
    newdataEvent = sim$dataEvent,
    process = "event",
    times = times_map,
    time_start = time_start_map,
    control = list(n_samples = 50),
    seed = 2401
  )
  c("SUCCESS", paste("nrow_survival=", nrow(pred$predictions$survival), sep = ""))
}, error = function(e) {
  c("ERROR", conditionMessage(e))
})
writeLines(result, out_file)
