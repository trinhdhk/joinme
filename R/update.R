#' Update a JoiNMe fit
#'
#' @description
#' Refits a JoiNMe model using the stored call, with optional updates to formulas,
#' data, and control arguments.
#'
#' @param object A JoiNMe fit object.
#' @param formulaLong Optional updated longitudinal formula. Use an update formula
#'   (e.g., `~ . + x`) to modify the existing model.
#' @param dataLong Optional updated longitudinal dataset.
#' @param formulaEvent Optional updated survival formula (full or update form).
#' @param dataEvent Optional updated event dataset.
#' @param formulaVCov Optional updated covariance regression. Supply one formula
#'   to update both components, or `list(sd = ~ ..., corr = ~ ...)` to update
#'   the standard-deviation and off-diagonal correlation regressions
#'   independently. Update formulae containing `.` are supported within each
#'   component.
#' @param formulaDist Optional distributional regression formulas with parameter
#'   names on the LHS (e.g., `sigma ~ 1 + time`).
#' @param control Optional updated control list.
#' @param draws Optional draws override.
#' @param families Optional updated families specification.
#' @param transforms Optional updated transforms specification.
#' @param priors Optional updated priors list.
#' @param .env Optional environment for evaluating the updated call. 
#' If NULL, the parent frame is used. When failed, the environment of the original formulaLong is used.
#' @param ... Additional arguments passed to `joinme()`, or to `joinme_mix()`
#'   when `object` is a latent-progress mixture.
#'
#' @return A refitted JoiNMe object. Mixture fits retain their mixture entry
#'   point and class specification. Any longitudinal-only fit remains
#'   longitudinal only unless new event inputs are supplied explicitly.
#' @method update JoiNMeFit
#' @seealso [update()]
#' @export
update.JoiNMeFit <- function(
  object,
  formulaLong = NULL,
  dataLong = NULL,
  formulaEvent = NULL,
  dataEvent = NULL,
  formulaVCov = NULL,
  formulaDist = NULL,
  control = NULL,
  draws = NULL,
  families = NULL,
  transforms = NULL,
  priors = NULL,
  .env = NULL,
  ...
) {
  mixture_fit <- inherits(object, "JoiNMeMixFit")
  longitudinal_only_fit <- !.fit_includes_survival(
    object
  ) # whether the stored call intentionally omitted the event process
  call_obj <- object$call
  if (is.null(call_obj) || !is.call(call_obj)) {
    call_obj <- call(
      if (mixture_fit) "joinme_mix" else "joinme"
    )
  }
  call_obj[[1]] <- if (mixture_fit) {
    quote(joinme_mix)
  } else {
    quote(joinme)
  }

  update_formula <- function(current, updated, name) {
    # Support update formulas ("~ . + x") while preserving original
    if (is.null(updated)) return(current)
    if (!inherits(updated, "formula")) {
      cli::cli_abort("{.arg {name}} must be a formula.")
    }

    # Robust dot detection without regex; the dot symbol is not reported by all.vars().
    expr_has_dot <- function(expr) {
      if (is.null(expr)) return(FALSE)
      if (is.name(expr) && identical(as.character(expr), ".")) return(TRUE)
      if (!is.call(expr)) return(FALSE)
      any(vapply(as.list(expr)[-1], expr_has_dot, logical(1)))
    }

    replace_dot <- function(expr, replacement) {
      if (is.null(expr)) return(expr)
      if (is.name(expr) && identical(as.character(expr), ".")) return(replacement)
      if (!is.call(expr)) return(expr)
      as.call(c(expr[[1]], lapply(as.list(expr)[-1], replace_dot, replacement = replacement)))
    }

    has_lhs <- length(updated) == 3
    rhs <- if (has_lhs) updated[[3]] else updated[[2]]
    if (expr_has_dot(rhs) || !has_lhs) {
      current_exp <- reformulas::expandDoubleVerts(current)
      current_bars <- reformulas::findbars(current_exp)
      current_fix <- reformulas::nobars(current_exp)
      current_fix_rhs <- if (length(current_fix) == 3) current_fix[[3]] else current_fix[[2]]
      new_rhs <- if (expr_has_dot(rhs)) replace_dot(rhs, current_fix_rhs) else rhs

      tmp_form <- rlang::new_formula(NULL, new_rhs, env = rlang::`%||%`(environment(current), environment(updated)))
      has_bars <- length(reformulas::findbars(reformulas::expandDoubleVerts(tmp_form))) > 0
      if (!has_bars && length(current_bars) > 0) {
        new_rhs <- Reduce(function(acc, bar) call("+", acc, bar), current_bars, init = new_rhs)
      }

      lhs <- if (has_lhs) updated[[2]] else if (length(current) == 3) current[[2]] else NULL
      return(rlang::new_formula(lhs, new_rhs, env = rlang::`%||%`(environment(current), environment(updated))))
    }
    updated
  }

  call_obj$formulaLong <- update_formula(object$formulaLong, formulaLong, "formulaLong")
  call_obj$formulaEvent <- if (
    longitudinal_only_fit &&
      is.null(formulaEvent)
  ) {
    NULL
  } else {
    update_formula(
      object$formulaEvent,
      formulaEvent,
      "formulaEvent"
    )
  }
  update_vcov_formula <- function(current, updated) {
    current_components <- .get_vcov_formula(
      current,
      default = ~ 1,
      context = "update.JoiNMeFit()"
    ) # canonical SD/correlation formula pair stored by all current fit objects
    if (is.null(updated)) return(current_components)
    if (inherits(updated, "formula")) {
      return(list(
        sd = update_formula(current_components$sd, updated, "formulaVCov$sd"),
        corr = update_formula(current_components$corr, updated, "formulaVCov$corr")
      ))
    }
    updated_components <- .get_vcov_formula(
      updated,
      default = ~ 1,
      context = "update.JoiNMeFit()"
    )
    list(
      sd = update_formula(current_components$sd, updated_components$sd, "formulaVCov$sd"),
      corr = update_formula(current_components$corr, updated_components$corr, "formulaVCov$corr")
    )
  }
  call_obj$formulaVCov <- update_vcov_formula(object$formulaVCov, formulaVCov)
  if (!is.null(formulaDist)) call_obj$formulaDist <- formulaDist

  if (!is.null(call_obj$formulaLong) && inherits(call_obj$formulaLong, "formula")) {
    f_exp <- reformulas::expandDoubleVerts(call_obj$formulaLong)
    if (length(reformulas::findbars(f_exp)) == 0 && !is.null(object$formulaLong)) {
      current_exp <- reformulas::expandDoubleVerts(object$formulaLong)
      current_bars <- reformulas::findbars(current_exp)
      if (length(current_bars) > 0) {
        fe_form <- reformulas::nobars(f_exp)
        fe_rhs <- if (length(fe_form) == 3) fe_form[[3]] else fe_form[[2]]
        rhs <- Reduce(function(acc, bar) call("+", acc, bar), current_bars, init = fe_rhs)
        lhs <- if (length(call_obj$formulaLong) == 3) call_obj$formulaLong[[2]] else object$formulaLong[[2]]
        call_obj$formulaLong <- rlang::new_formula(lhs, rhs, env = rlang::`%||%`(environment(call_obj$formulaLong), environment(object$formulaLong)))
      }
    }
  }

  if (!is.null(dataLong)) call_obj$dataLong <- dataLong
  if (!is.null(dataEvent)) call_obj$dataEvent <- dataEvent
  if (!is.null(control)) call_obj$control <- control
  if (!is.null(draws)) call_obj$draws <- draws
  if (!is.null(families)) call_obj$families <- families
  if (!is.null(transforms)) call_obj$transforms <- transforms
  if (!is.null(priors)) call_obj$priors <- priors

  dots <- list(...)
  if (length(dots) > 0) {
    if (is.null(names(dots)) || any(names(dots) == "")) {
      cli::cli_abort("All additional arguments must be named.")
    }
    for (nm in names(dots)) {
      call_obj[[nm]] <- dots[[nm]]
    }
  }

  if (is.null(call_obj$dataLong)) {
    cli::cli_abort("{.arg dataLong} must be supplied when it is not available in the stored call.")
  }
  if (
    is.null(call_obj$dataEvent) &&
      !longitudinal_only_fit
  ) {
    cli::cli_abort("{.arg dataEvent} must be supplied when it is not available in the stored call.")
  }

  if (is.null(.env)) {
    .env <- parent.frame()
  }
  tryCatch(eval(call_obj, .env), 
    error = function(...){
        eval(call_obj, environment(object$call$formulaLong))
  })
}
