#' @name bytecode_interpreter_module
#' @title Portable stack-bytecode interpreter
#'
#' @description
#' A model-agnostic interpreter for scalar transformation programmes encoded in
#' reverse Polish notation.  This file deliberately depends only on base R and
#' the stats package.  It does not inspect JoiNMe formulae, fitted objects or
#' Stan data, which allows the file to be moved into a small standalone package
#' without changing its execution contract.
#'
#' The authoritative instruction protocol is documented in
#' `inst/stan/include/etc/bytecode/README.md` and mirrored by the Stan module under
#' `inst/stan/include/etc/bytecode`.
#'
#' @keywords internal
NULL

#' Bytecode instruction registry
#'
#' @description
#' Defines the language-neutral integer protocol shared by the R and Stan
#' interpreters.  A named registry makes callers independent of package model
#' objects and gives a future standalone package one place to extend the
#' instruction set.
#'
#' @return A named integer vector.
#' @keywords internal
#' @noRd
.bytecode_opcodes <- function() {
  c(
    PUSH_X = 0L, PUSH_CONST = 1L,
    ADD = 2L, SUB = 3L, MUL = 4L, DIV = 5L,
    LOG = 6L, EXP = 7L, SQRT = 8L, INV_LOGIT = 9L, LOGIT = 10L,
    RECIPROCAL = 11L, POW = 12L, SIN = 13L, COS = 14L, TAN = 15L,
    ABS = 16L, SQUARE = 17L, SINH = 18L, COSH = 19L, TANH = 20L,
    ASINH = 21L, ACOSH = 22L, ATANH = 23L, SOFTPLUS = 24L,
    CBRT = 25L, PHI = 26L, INV_PHI = 27L
  )
}

#' Standard-normal bytecode instruction identifiers
#'
#' @description
#' Provides the single R-side definition of the two standard-normal
#' transformation instructions. `PHI` maps a real-valued variate to a
#' probability through the standard normal cumulative distribution function.
#' `INV_PHI` maps a probability to its standard normal quantile. Keeping these
#' identifiers together prevents the forward probit link from being confused
#' with its inverse link.
#'
#' @return A named integer vector with elements `PHI = 26L` and
#'   `INV_PHI = 27L`.
#' @keywords internal
#' @noRd
.bytecode_normal_ops <- function() {
  .bytecode_opcodes()[c("PHI", "INV_PHI")]
}

#' Verify functional bytecode stack behaviour
#'
#' @description
#' Visits the instruction stream without evaluating numerical values. The
#' verifier checks that each unary or binary operation has enough operands,
#' that every constant instruction has a corresponding value, that every
#' instruction is recognised, and that exactly one result remains.
#'
#' @param bytecode Integer bytecode sequence.
#' @param const_data Numeric constants consumed successively by `PUSH_CONST`.
#'
#' @return `TRUE` invisibly when the program is structurally valid; otherwise
#'   an error is raised.
#' @keywords internal
#' @noRd
verify_bytecode <- function(bytecode, const_data) {
  stack_height <- 0
  const_idx <- 0
  unary_ops <- .bytecode_unary_ops()

  for (op in bytecode) {
    if (op == 0) {
      # PUSH_X
      stack_height <- stack_height + 1
    } else if (op == 1) {
      # PUSH_CONST
      const_idx <- const_idx + 1
      if (const_idx > length(const_data)) {
        stop("Functional bytecode references more constants than provided")
      }
      stack_height <- stack_height + 1
    } else if (op %in% unary_ops) {
      # Unary operations
      if (stack_height < 1) stop("Stack underflow: unary operation")
    } else if (op %in% c(2, 3, 4, 5, 12)) {
      # Binary operations
      if (stack_height < 2) stop("Stack underflow: binary operation")
      stack_height <- stack_height - 1
    } else {
      stop("Unknown functional bytecode instruction: ", op)
    }
  }

  if (stack_height != 1) {
    stop(paste("Final stack height is", stack_height, "; expected 1"))
  }

  invisible(TRUE)
}

#' Enumerate unary bytecode instructions
#'
#' @description
#' Returns the operation identifiers that consume one stack value and replace
#' it with one transformed value. The same list is used by validation and by
#' earlier one-operation program normalisation.
#'
#' @return Integer vector of unary operation identifiers.
#' @keywords internal
#' @noRd
.bytecode_unary_ops <- function() {
  c(
    6L, 7L, 8L, 9L, 10L, 11L, 13L, 14L, 15L, 16L, 17L, 18L, 19L,
    20L, 21L, 22L, 23L, 24L, 25L, unname(.bytecode_normal_ops())
  )
}

#' Evaluate bytecode for a scalar input
#'
#' @description
#' Executes a stack-based bytecode program for a single numeric input.
#' This evaluator mirrors the Stan-side instruction semantics exactly so
#' simulation, preprocessing checks, and Stan likelihood evaluation remain aligned.
#'
#' @param x Numeric scalar input.
#' @param bytecode Integer bytecode sequence.
#' @param const_data Numeric constants consumed by `PUSH_CONST` instructions.
#' @param iota_intercepts Numeric affine intercepts referenced by individual
#'   nonlinear instructions.
#' @param iota_slopes Numeric affine slopes referenced by individual nonlinear
#'   instructions.
#' @param op_iota_intercept_idx Integer instruction-aligned indices into
#'   `iota_intercepts`; zero means no fitted intercept.
#' @param op_iota_slope_idx Integer instruction-aligned indices into
#'   `iota_slopes`; zero means no fitted slope.
#'
#' @return Numeric scalar result.
#' @keywords internal
#' @noRd
eval_bytecode_scalar <- function(
  x,
  bytecode = NULL,
  const_data = numeric(),
  iota_intercepts = numeric(0),
  iota_slopes = numeric(0),
  op_iota_intercept_idx = NULL,
  op_iota_slope_idx = NULL
) {
  # Step 1: normalise inputs.
  program <- .normalize_bytecode_program(
    bytecode = bytecode,
    const_data = const_data,
    op_iota_intercept_idx = op_iota_intercept_idx,
    op_iota_slope_idx = op_iota_slope_idx
  )
  code <- program$bytecode
  constants <- program$const_data
  intercept_idx <- program$op_iota_intercept_idx
  slope_idx <- program$op_iota_slope_idx
  normal_ops <- .bytecode_normal_ops()

  # Identity behaviour for empty programs keeps backward compatibility.
  if (length(code) == 0L) {
    return(as.numeric(x))
  }

  # Step 2: validate stack behaviour before execution for reproducibility.
  verify_bytecode(code, constants)

  # Step 3: execute the bytecode stack machine instruction by instruction.
  stack <- numeric(0)
  const_idx <- 1L
  resolve_iota_shift <- function(value, intercept_ref, slope_ref) {
    intercept_val <- if (isTRUE(intercept_ref > 0L) && length(iota_intercepts) >= intercept_ref) iota_intercepts[[intercept_ref]] else 0
    slope_val <- if (isTRUE(slope_ref > 0L) && length(iota_slopes) >= slope_ref) iota_slopes[[slope_ref]] else 1
    as.numeric(intercept_val + slope_val * value)
  }
  for (op_pos in seq_along(code)) {
    op <- code[[op_pos]]
    op_intercept_idx <- intercept_idx[[op_pos]]
    op_slope_idx <- slope_idx[[op_pos]]
    if (op == 0L) {
      stack <- c(stack, as.numeric(x))
    } else if (op == 1L) {
      stack <- c(stack, constants[const_idx])
      const_idx <- const_idx + 1L
    } else if (op == 2L) {
      b <- stack[length(stack)]; a <- stack[length(stack) - 1L]
      stack <- c(stack[-c(length(stack) - 1L, length(stack))], a + b)
    } else if (op == 3L) {
      b <- stack[length(stack)]; a <- stack[length(stack) - 1L]
      stack <- c(stack[-c(length(stack) - 1L, length(stack))], a - b)
    } else if (op == 4L) {
      b <- stack[length(stack)]; a <- stack[length(stack) - 1L]
      stack <- c(stack[-c(length(stack) - 1L, length(stack))], a * b)
    } else if (op == 5L) {
      b <- stack[length(stack)]; a <- stack[length(stack) - 1L]
      stack <- c(stack[-c(length(stack) - 1L, length(stack))], a / b)
    } else if (op == 6L) {
      stack[length(stack)] <- log(resolve_iota_shift(stack[length(stack)], op_intercept_idx, op_slope_idx))
    } else if (op == 7L) {
      stack[length(stack)] <- exp(resolve_iota_shift(stack[length(stack)], op_intercept_idx, op_slope_idx))
    } else if (op == 8L) {
      stack[length(stack)] <- sqrt(resolve_iota_shift(stack[length(stack)], op_intercept_idx, op_slope_idx))
    } else if (op == 9L) {
      stack[length(stack)] <- stats::plogis(resolve_iota_shift(stack[length(stack)], op_intercept_idx, op_slope_idx))
    } else if (op == 10L) {
      stack[length(stack)] <- stats::qlogis(resolve_iota_shift(stack[length(stack)], op_intercept_idx, op_slope_idx))
    } else if (op == 11L) {
      stack[length(stack)] <- 1 / resolve_iota_shift(stack[length(stack)], op_intercept_idx, op_slope_idx)
    } else if (op == 12L) {
      b <- stack[length(stack)]; a <- stack[length(stack) - 1L]
      a <- resolve_iota_shift(a, op_intercept_idx, op_slope_idx)
      stack <- c(stack[-c(length(stack) - 1L, length(stack))], a^b)
    } else if (op == 13L) {
      stack[length(stack)] <- sin(resolve_iota_shift(stack[length(stack)], op_intercept_idx, op_slope_idx))
    } else if (op == 14L) {
      stack[length(stack)] <- cos(resolve_iota_shift(stack[length(stack)], op_intercept_idx, op_slope_idx))
    } else if (op == 15L) {
      stack[length(stack)] <- tan(resolve_iota_shift(stack[length(stack)], op_intercept_idx, op_slope_idx))
    } else if (op == 16L) {
      stack[length(stack)] <- abs(resolve_iota_shift(stack[length(stack)], op_intercept_idx, op_slope_idx))
    } else if (op == 17L) {
      stack[length(stack)] <- resolve_iota_shift(stack[length(stack)], op_intercept_idx, op_slope_idx)^2
    } else if (op == 18L) {
      stack[length(stack)] <- sinh(resolve_iota_shift(stack[length(stack)], op_intercept_idx, op_slope_idx))
    } else if (op == 19L) {
      stack[length(stack)] <- cosh(resolve_iota_shift(stack[length(stack)], op_intercept_idx, op_slope_idx))
    } else if (op == 20L) {
      stack[length(stack)] <- tanh(resolve_iota_shift(stack[length(stack)], op_intercept_idx, op_slope_idx))
    } else if (op == 21L) {
      stack[length(stack)] <- asinh(resolve_iota_shift(stack[length(stack)], op_intercept_idx, op_slope_idx))
    } else if (op == 22L) {
      stack[length(stack)] <- acosh(resolve_iota_shift(stack[length(stack)], op_intercept_idx, op_slope_idx))
    } else if (op == 23L) {
      stack[length(stack)] <- atanh(resolve_iota_shift(stack[length(stack)], op_intercept_idx, op_slope_idx))
    } else if (op == 24L) {
      shifted_value <- resolve_iota_shift(
        stack[length(stack)],
        op_intercept_idx,
        op_slope_idx
      )
      # Numerically stable log(1 + exp(x)); expressed with base R so this
      # interpreter has no dependency on a host package's utility functions.
      stack[length(stack)] <-
        log1p(exp(-abs(shifted_value))) + max(shifted_value, 0)
    } else if (op == 25L) {
      a <- resolve_iota_shift(stack[length(stack)], op_intercept_idx, op_slope_idx)
      stack[length(stack)] <- sign(a) * abs(a)^(1 / 3)
    } else if (op == normal_ops[["PHI"]]) {
      stack[length(stack)] <- stats::pnorm(resolve_iota_shift(stack[length(stack)], op_intercept_idx, op_slope_idx))
    } else if (op == normal_ops[["INV_PHI"]]) {
      stack[length(stack)] <- stats::qnorm(resolve_iota_shift(stack[length(stack)], op_intercept_idx, op_slope_idx))
    } else {
      stop("Unknown transform bytecode instruction: ", op, call. = FALSE)
    }
  }

  # Step 4: return final stack top as the transformation output.
  stack[length(stack)]
}

#' Evaluate a canonical bytecode program without scalar dispatch
#'
#' @description
#' Recognises the most frequent one-operation transformations and evaluates
#' them directly on a numeric vector. This avoids repeatedly normalising and
#' verifying the same program for every observation. The shortcut is used only
#' when the program has no constants or fitted affine shifts; all other
#' programs retain the general stack evaluator.
#'
#' @param x Numeric vector.
#' @param program Normalised bytecode program returned by
#'   `.normalize_bytecode_program()`.
#'
#' @return A numeric vector for a recognised canonical program, or `NULL` when
#'   the general evaluator is required.
#' @keywords internal
#' @noRd
.eval_canonical_bytecode_vector <- function(x, program) {
  code <- as.integer(program$bytecode)
  if (length(program$const_data) > 0L ||
      any(program$op_iota_intercept_idx > 0L) ||
      any(program$op_iota_slope_idx > 0L)) {
    return(NULL)
  }

  x <- as.numeric(x)
  if (identical(code, 0L)) return(x)
  if (identical(code, c(0L, 6L))) return(log(x))
  if (identical(code, c(0L, 7L))) return(exp(x))
  if (identical(code, c(0L, 9L))) return(stats::plogis(x))

  normal_ops <- .bytecode_normal_ops()
  if (identical(code, c(0L, unname(normal_ops[["PHI"]])))) {
    return(stats::pnorm(x))
  }
  if (identical(code, c(0L, unname(normal_ops[["INV_PHI"]])))) {
    return(stats::qnorm(x))
  }
  NULL
}

#' Evaluate bytecode for vector inputs
#'
#' @description
#' Applies `eval_bytecode_scalar()` element-wise to a numeric vector.
#'
#' @param x Numeric vector input.
#' @param bytecode Integer bytecode sequence.
#' @param const_data Numeric constants consumed by `PUSH_CONST` instructions.
#' @param iota_intercepts Numeric affine intercepts referenced by nonlinear
#'   instructions.
#' @param iota_slopes Numeric affine slopes referenced by nonlinear
#'   instructions.
#' @param op_iota_intercept_idx Integer instruction-aligned intercept indices.
#' @param op_iota_slope_idx Integer instruction-aligned slope indices.
#'
#' @return Numeric vector result.
#' @keywords internal
#' @noRd
eval_bytecode_vector <- function(
  x,
  bytecode = NULL,
  const_data = numeric(),
  iota_intercepts = numeric(0),
  iota_slopes = numeric(0),
  op_iota_intercept_idx = NULL,
  op_iota_slope_idx = NULL
) {
  program <- .normalize_bytecode_program(
    bytecode = bytecode,
    const_data = const_data,
    op_iota_intercept_idx = op_iota_intercept_idx,
    op_iota_slope_idx = op_iota_slope_idx
  )
  canonical <- .eval_canonical_bytecode_vector(x, program)
  if (!is.null(canonical)) {
    return(canonical)
  }
  vapply(
    as.numeric(x),
    eval_bytecode_scalar,
    numeric(1),
    bytecode = program$bytecode,
    const_data = program$const_data,
    iota_intercepts = iota_intercepts,
    iota_slopes = iota_slopes,
    op_iota_intercept_idx = program$op_iota_intercept_idx,
    op_iota_slope_idx = program$op_iota_slope_idx
  )
}

#' Normalise bytecode program inputs
#'
#' @param bytecode Integer bytecode sequence or `NULL`.
#' @param const_data Numeric constant vector.
#' @return List with normalised `bytecode` and `const_data`.
#' @keywords internal
#' @noRd
.normalize_bytecode_program <- function(
  bytecode = NULL,
  const_data = numeric(),
  op_iota_intercept_idx = NULL,
  op_iota_slope_idx = NULL
) {
  code <- bytecode
  if (is.null(code)) {
    code <- integer(0)
  }

  code <- as.integer(code)
  constants <- as.numeric(
    if (is.null(const_data)) numeric(0) else const_data
  )
  iota_intercepts <- as.integer(
    if (is.null(op_iota_intercept_idx)) {
      rep(0L, length(code))
    } else {
      op_iota_intercept_idx
    }
  )
  iota_slopes <- as.integer(
    if (is.null(op_iota_slope_idx)) {
      rep(0L, length(code))
    } else {
      op_iota_slope_idx
    }
  )

  # Backward compatibility: some earlier saved transforms encoded unary maps like
  # softplus(x) as [SOFTPLUS] rather than [PUSH_X, SOFTPLUS]. Prepend PUSH_X so
  # replay in R and Stan stays stable for old fitted objects used in prediction.
  if (length(code) > 0L && code[[1]] %in% .bytecode_unary_ops()) {
    code <- c(0L, code)
    iota_intercepts <- c(0L, iota_intercepts)
    iota_slopes <- c(0L, iota_slopes)
  }

  if (length(iota_intercepts) != length(code) || length(iota_slopes) != length(code)) {
    stop(
      "Bytecode affine-index arrays must align with the instruction count.",
      call. = FALSE
    )
  }

  list(
    bytecode = code,
    const_data = constants,
    op_iota_intercept_idx = iota_intercepts,
    op_iota_slope_idx = iota_slopes
  )
}
