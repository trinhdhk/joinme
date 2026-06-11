devtools::load_all(quiet = TRUE)

env <- new.env(parent = emptyenv())
load("tmp_fit-sim3.Rdata", envir = env)
fit <- env[["fit"]]
sim <- env[["sim"]]

id_var <- eval(fit$call$id_var) %||% "id"
id1 <- unique(sim$dataEvent[[id_var]])[1]

dE <- sim$dataEvent[sim$dataEvent[[id_var]] == id1, , drop = FALSE]
dL <- sim$dataLong[sim$dataLong[[id_var]] == id1, , drop = FALSE]

run_case <- function(times_case, label) {
  cat("== ", label, " ==\n", sep = "")
  out <- tryCatch(
    predict(
      fit,
      newdataLong = dL,
      newdataEvent = dE,
      process = "event",
      times = stats::setNames(list(times_case), id1),
      time_start = stats::setNames(1, id1),
      control = list(n_samples = 50),
      seed = 2401
    ),
    error = function(e) e
  )
  if (inherits(out, "error")) {
    cat("ERROR\n")
    cat(conditionMessage(out), "\n")
  } else {
    cat("SUCCESS\n")
    print(utils::head(out$predictions$survival))
  }
}

run_case(3, "single-point horizon")
run_case(seq(1, 3, length.out = 50), "dense grid")