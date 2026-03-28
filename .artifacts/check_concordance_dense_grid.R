devtools::load_all(quiet = TRUE)

env <- new.env(parent = emptyenv())
load("tmp_fit-sim3.Rdata", envir = env)
fit <- env[["fit"]]
sim <- env[["sim"]]

id_var <- eval(fit$call$id_var) %||% "id"
time_var <- eval(fit$call$time_var) %||% "time"
ids <- unique(sim$dataEvent[[id_var]])
ids <- ids[!is.na(ids)]
time_start_map <- stats::setNames(rep(1, length(ids)), ids)
time_horizon_map <- stats::setNames(rep(3, length(ids)), ids)
times_map <- stats::setNames(lapply(ids, function(id) seq(time_start_map[[as.character(id)]], time_horizon_map[[as.character(id)]], length.out = 50)), ids)

res <- tryCatch(
  predict(
    fit,
    newdataLong = sim$dataLong,
    newdataEvent = sim$dataEvent,
    process = "event",
    times = times_map,
    time_start = time_start_map,
    control = list(n_samples = 50),
    seed = 2401
  ),
  error = function(e) {
    cat("ERROR:\n")
    cat(conditionMessage(e), "\n")
    NULL
  }
)

if (!is.null(res)) {
  cat("SUCCESS\n")
  print(head(res$predictions$survival))
}
