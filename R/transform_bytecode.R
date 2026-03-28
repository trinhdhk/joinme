##' Functional Transform Builder for Joint Models
##'
##' Provides R utilities to parse user-friendly transformation expressions
##' and generate functional bytecode for Stan's arbitrary transformation evaluator.
##'
##' @description
##' The `parse_transform_expr` function takes a formula, quosure, quoted expression,
##' or character string and converts it to a functional bytecode representation compatible with the
##' Stan-side opcode evaluator.
##'
##' Supported operations:
##' - Arithmetic: +, -, *, /
##' - Power: ^
##' - Functions: log, exp, sqrt, sin, cos, tan, abs, sinh, cosh, tanh, asinh, acosh, atanh
##' - Sigmoid/Link: inv_logit, logit, sigmoid, expit
##' - Softplus: softplus, log1p_exp
##' - Other: cbrt, power
##' - Reciprocal: 1/x or rec(x)
##'
##' @param expr A formula (e.g. `~ x + 2`), quosure, quoted expression, or character string to parse.
##'
##' @note
##' Functional bytecode is a vector of operation codes (0-25) paired with a vector of
##' constant values. Constants are embedded using PUSH_CONST operations.
##'
##' @section Functional Bytecode Reference:
##' \describe{
##'   \item{0}{PUSH_X: Push input x onto stack}
##'   \item{1}{PUSH_CONST: Push next constant value}
##'   \item{2}{ADD: Pop b,a; push a+b}
##'   \item{3}{SUB: Pop b,a; push a-b}
##'   \item{4}{MUL: Pop b,a; push a*b}
##'   \item{5}{DIV: Pop b,a; push a/b}
##'   \item{6}{LOG: Pop a; push log(a)}
##'   \item{7}{EXP: Pop a; push exp(a)}
##'   \item{8}{SQRT: Pop a; push sqrt(a)}
##'   \item{9}{INV_LOGIT: Pop a; push inv_logit(a)}
##'   \item{10}{LOGIT: Pop a; push logit(a)}
##'   \item{11}{RECIPROCAL: Pop a; push 1/a}
##'   \item{12}{POW: Pop b,a; push a^b}
##'   \item{13}{SIN: Pop a; push sin(a)}
##'   \item{14}{COS: Pop a; push cos(a)}
##'   \item{15}{TAN: Pop a; push tan(a)}
##'   \item{16}{ABS: Pop a; push abs(a)}
##'   \item{17}{SQUARE: Pop a; push a^2}
##'   \item{18}{SINH: Pop a; push sinh(a)}
##'   \item{19}{COSH: Pop a; push cosh(a)}
##'   \item{20}{TANH: Pop a; push tanh(a)}
##'   \item{21}{ASINH: Pop a; push asinh(a)}
##'   \item{22}{ACOSH: Pop a; push acosh(a)}
##'   \item{23}{ATANH: Pop a; push atanh(a)}
##'   \item{24}{SOFTPLUS: Pop a; push log1p_exp(a)}
##'   \item{25}{CBRT: Pop a; push cbrt(a)}
##'   \item{26}{PROBIT: Pop a; push Phi(a)}
##' }
##'
##' @examples
##' # Simple: f(x) = x + 2
##' bc <- parse_transform_expr(~ x + 2)
##'
##' # Use a plain formula (recommended)
##' bc <- parse_transform_expr(~ log(x + 1))
##' # Returns: bytecode = c(0, 1, 2), const_data = c(2)
##'
##' # Complex: f(x) = (log(sqrt(x + 1/inv_logit(3*x - 3))))^2
##' bc <- parse_transform_expr(~ (log(sqrt(x + 1/inv_logit(3*x - 3))))^2)
##'
##' @export
# File overview:
# - Parse R expressions into portable bytecode programs.
# - Validate bytecode stack behaviour before handing programs to Stan.
# - Evaluate bytecode in R (shared by inverse-link and association transforms).
parse_transform_expr <- function(expr) {
  # Step 1: normalise any supported expression input into a language object.
  expr_call <- .coerce_transform_expr(expr)

  # Step 2: emit bytecode instructions and constant vector.
  bytecode <- integer()
  const_data <- numeric()
  result <- .emit_bytecode_expr(expr_call, bytecode, const_data)

  # Step 3: validate stack consistency so runtime evaluation is deterministic.
  verify_opcodes(result$bytecode, result$const_data)

  # Step 4: return both `bytecode` and legacy `opcodes` aliases for compatibility.
  list(
    opcodes = result$bytecode,
    bytecode = result$bytecode,
    const_data = result$const_data,
    n_ops = length(result$bytecode),
    n_bytecode = length(result$bytecode),
    n_const = length(result$const_data)
  )
}

##' @keywords internal
##' Coerce input to a language object using base R parsing.
.coerce_transform_expr <- function(expr) {
  # Normalise input to a single language object
  if (rlang::is_quosure(expr)) {
    rhs <- rlang::expr_text(rlang::get_expr(expr))
    expr <- stats::as.formula(paste0("~ ", rhs))
  }
  if (inherits(expr, "formula")) {
    return(expr[[2]])
  }
  if (is.expression(expr)) {
    if (length(expr) == 0) cli::cli_abort(c(x = "{.arg expr} is empty.", i = "Provide a valid expression or formula."))
    return(expr[[1]])
  }
  if (is.character(expr)) {
    parsed <- tryCatch(parse(text = expr), error = function(e) e)
    if (inherits(parsed, "error")) {
      cli::cli_abort(c(
        x = "Failed to parse {.arg expr} as an R expression.",
        i = "Provide a formula like ~ log(x) or a valid expression string."
      ))
    }
    if (length(parsed) == 0) cli::cli_abort(c(x = "{.arg expr} is empty.", i = "Provide a valid expression or formula."))
    return(parsed[[1]])
  }
  if (is.call(expr) || is.name(expr) || is.numeric(expr)) {
    return(expr)
  }
  cli::cli_abort(c(
    x = "Unsupported transform expression type.",
    i = "Use a formula (e.g., ~ log(x)), a quoted expression, or a character string."
  ))
}

##' @keywords internal
##' Emit bytecode from an R language object (calls, names, constants).
.emit_bytecode_expr <- function(node, bytecode, const_data) {
  # Recursive descent over AST nodes to emit bytecode
  if (is.numeric(node)) {
    bytecode <- c(bytecode, 1L)
    const_data <- c(const_data, as.numeric(node))
    return(list(bytecode = bytecode, const_data = const_data))
  }
  if (is.name(node)) {
    if (as.character(node) != "x") {
      cli::cli_abort(c(
        x = "Unknown symbol in transform expression: {as.character(node)}.",
        i = "Use 'x' as the transformation input variable."
      ))
    }
    bytecode <- c(bytecode, 0L)
    return(list(bytecode = bytecode, const_data = const_data))
  }
  if (!is.call(node)) {
    cli::cli_abort(c(
      x = "Unsupported node in transform expression.",
      i = "Use standard R arithmetic and function calls."
    ))
  }
  
  op <- as.character(node[[1]])
  args <- as.list(node[-1])

  if (op == "(") {
    if (length(args) != 1) {
      cli::cli_abort(c(
        x = "Parentheses in transform expression must wrap a single expression.",
        i = "Remove stray commas or extra arguments inside parentheses."
      ))
    }
    return(.emit_bytecode_expr(args[[1]], bytecode, const_data))
  }
  
  if (op %in% c("+", "-", "*", "/", "^")) {
    if (length(args) == 1 && op == "-") {
      res <- .emit_bytecode_expr(0, bytecode, const_data)
      bytecode <- res$bytecode
      const_data <- res$const_data
      res <- .emit_bytecode_expr(args[[1]], bytecode, const_data)
      bytecode <- res$bytecode
      const_data <- res$const_data
      bytecode <- c(bytecode, 3L)
      return(list(bytecode = bytecode, const_data = const_data))
    }
    if (length(args) != 2) {
      cli::cli_abort(c(
        x = "Binary operator {op} expects two arguments.",
        i = "Check the transform expression for missing operands."
      ))
    }
    res <- .emit_bytecode_expr(args[[1]], bytecode, const_data)
    bytecode <- res$bytecode
    const_data <- res$const_data
    res <- .emit_bytecode_expr(args[[2]], bytecode, const_data)
    bytecode <- res$bytecode
    const_data <- res$const_data
    op_id <- switch(op, "+" = 2L, "-" = 3L, "*" = 4L, "/" = 5L, "^" = 12L)
    bytecode <- c(bytecode, op_id)
    return(list(bytecode = bytecode, const_data = const_data))
  }

  if (op == "power") {
    if (length(args) != 2) {
      cli::cli_abort(c(
        x = "Function {op} expects two arguments.",
        i = "Check the transform expression for missing arguments."
      ))
    }
    res <- .emit_bytecode_expr(args[[1]], bytecode, const_data)
    bytecode <- res$bytecode
    const_data <- res$const_data
    res <- .emit_bytecode_expr(args[[2]], bytecode, const_data)
    bytecode <- res$bytecode
    const_data <- res$const_data
    bytecode <- c(bytecode, 12L)
    return(list(bytecode = bytecode, const_data = const_data))
  }
  
  func_op_id <- switch(op,
    log = 6L,
    exp = 7L,
    sqrt = 8L,
    inv_logit = 9L,
    sigmoid = 9L,
    expit = 9L,
    logit = 10L,
    rec = 11L,
    sin = 13L,
    cos = 14L,
    tan = 15L,
    abs = 16L,
    sinh = 18L,
    cosh = 19L,
    tanh = 20L,
    asinh = 21L,
    acosh = 22L,
    atanh = 23L,
    softplus = 24L,
    log1p_exp = 24L,
    cbrt = 25L,
    probit = 26L,
    NULL
  )
  if (!is.null(func_op_id)) {
    if (length(args) != 1) {
      cli::cli_abort(c(
        x = "Function {op} expects one argument.",
        i = "Check the transform expression for missing arguments."
      ))
    }
    res <- .emit_bytecode_expr(args[[1]], bytecode, const_data)
    bytecode <- res$bytecode
    const_data <- res$const_data
    bytecode <- c(bytecode, func_op_id)
    return(list(bytecode = bytecode, const_data = const_data))
  }
  
  if (op == "ISpline") {
    cli::cli_abort(c(
      x = "ISpline() is not supported in functional bytecode mode.",
      i = "Use the ispline transform mode and supply knots/coefficients instead."
    ))
  }
  
  cli::cli_abort(c(
    x = "Unsupported function in transform expression: {op}.",
    i = "Supported functions: log, exp, sqrt, inv_logit, logit, probit, sigmoid, expit, softplus, log1p_exp, cbrt, power, rec, sin, cos, tan, abs, sinh, cosh, tanh, asinh, acosh, atanh."
  ))
}

##' @keywords internal
##' Legacy tokenizer/parser removed in favor of R's formula parsing.

##' @keywords internal
##' Verify functional bytecode for basic sanity (stack under/overflow, etc.)
verify_opcodes <- function(opcodes, const_data) {
  # Accept both historical `opcodes` naming and modern `bytecode` naming.
  stack_height <- 0
  const_idx <- 0

  for (op in opcodes) {
    if (op == 0) {
      # PUSH_X
      stack_height <- stack_height + 1
    } else if (op == 1) {
      # PUSH_CONST
      const_idx <- const_idx + 1
      if (const_idx > length(const_data)) {
        stop("Functional opcodes reference more constants than provided")
      }
      stack_height <- stack_height + 1
    } else if (op %in% c(6, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26)) {
      # Unary operations
      if (stack_height < 1) stop("Stack underflow: unary operation")
    } else if (op %in% c(2, 3, 4, 5, 12)) {
      # Binary operations
      if (stack_height < 2) stop("Stack underflow: binary operation")
      stack_height <- stack_height - 1
    }
  }

  if (stack_height != 1) {
    stop(paste("Final stack height is", stack_height, "; expected 1"))
  }

  return(TRUE)
}

.bytecode_unary_ops <- function() {
  c(6L, 7L, 8L, 9L, 10L, 11L, 13L, 14L, 15L, 16L, 17L, 18L, 19L, 20L, 21L, 22L, 23L, 24L, 25L, 26L)
}

#' Evaluate bytecode for a scalar input
#'
#' @description
#' Executes a stack-based bytecode program for a single numeric input.
#' This evaluator mirrors the Stan-side instruction semantics exactly so
#' simulation, preprocessing checks, and Stan likelihood evaluation remain aligned.
#'
#' @param x Numeric scalar input.
#' @param bytecode Integer bytecode sequence. Legacy name `opcodes` is also accepted.
#' @param const_data Numeric constants consumed by `PUSH_CONST` instructions.
#'
#' @return Numeric scalar result.
#' @keywords internal
eval_bytecode_scalar <- function(x, bytecode = NULL, const_data = numeric(), opcodes = NULL) {
  # Step 1: normalise inputs and resolve legacy naming aliases.
  program <- .normalize_bytecode_program(bytecode = bytecode, opcodes = opcodes, const_data = const_data)
  code <- program$bytecode
  constants <- program$const_data

  # Identity behaviour for empty programs keeps backward compatibility.
  if (length(code) == 0L) {
    return(as.numeric(x))
  }

  # Step 2: validate stack behaviour before execution for reproducibility.
  verify_opcodes(code, constants)

  # Step 3: execute the bytecode stack machine instruction by instruction.
  stack <- numeric(0)
  const_idx <- 1L
  for (op in code) {
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
      stack[length(stack)] <- log(stack[length(stack)])
    } else if (op == 7L) {
      stack[length(stack)] <- exp(stack[length(stack)])
    } else if (op == 8L) {
      stack[length(stack)] <- sqrt(stack[length(stack)])
    } else if (op == 9L) {
      stack[length(stack)] <- stats::plogis(stack[length(stack)])
    } else if (op == 10L) {
      stack[length(stack)] <- stats::qlogis(stack[length(stack)])
    } else if (op == 11L) {
      stack[length(stack)] <- 1 / stack[length(stack)]
    } else if (op == 12L) {
      b <- stack[length(stack)]; a <- stack[length(stack) - 1L]
      stack <- c(stack[-c(length(stack) - 1L, length(stack))], a^b)
    } else if (op == 13L) {
      stack[length(stack)] <- sin(stack[length(stack)])
    } else if (op == 14L) {
      stack[length(stack)] <- cos(stack[length(stack)])
    } else if (op == 15L) {
      stack[length(stack)] <- tan(stack[length(stack)])
    } else if (op == 16L) {
      stack[length(stack)] <- abs(stack[length(stack)])
    } else if (op == 17L) {
      stack[length(stack)] <- stack[length(stack)]^2
    } else if (op == 18L) {
      stack[length(stack)] <- sinh(stack[length(stack)])
    } else if (op == 19L) {
      stack[length(stack)] <- cosh(stack[length(stack)])
    } else if (op == 20L) {
      stack[length(stack)] <- tanh(stack[length(stack)])
    } else if (op == 21L) {
      stack[length(stack)] <- asinh(stack[length(stack)])
    } else if (op == 22L) {
      stack[length(stack)] <- acosh(stack[length(stack)])
    } else if (op == 23L) {
      stack[length(stack)] <- atanh(stack[length(stack)])
    } else if (op == 24L) {
      stack[length(stack)] <- .softplus(stack[length(stack)])
    } else if (op == 25L) {
      a <- stack[length(stack)]
      stack[length(stack)] <- sign(a) * abs(a)^(1 / 3)
    } else if (op == 26L) {
      stack[length(stack)] <- stats::pnorm(stack[length(stack)])
    } else {
      cli::cli_abort("Unknown transform bytecode instruction: {op}.")
    }
  }

  # Step 4: return final stack top as the transformation output.
  stack[length(stack)]
}

#' Evaluate bytecode for vector inputs
#'
#' @description
#' Applies `eval_bytecode_scalar()` element-wise to a numeric vector.
#'
#' @param x Numeric vector input.
#' @param bytecode Integer bytecode sequence. Legacy name `opcodes` is also accepted.
#' @param const_data Numeric constants consumed by `PUSH_CONST` instructions.
#'
#' @return Numeric vector result.
#' @keywords internal
eval_bytecode_vector <- function(x, bytecode = NULL, const_data = numeric(), opcodes = NULL) {
  program <- .normalize_bytecode_program(bytecode = bytecode, opcodes = opcodes, const_data = const_data)
  vapply(as.numeric(x), eval_bytecode_scalar, numeric(1), bytecode = program$bytecode, const_data = program$const_data)
}

#' Normalise bytecode program inputs
#'
#' @param bytecode Integer bytecode sequence or `NULL`.
#' @param opcodes Legacy alias for bytecode.
#' @param const_data Numeric constant vector.
#' @return List with normalised `bytecode` and `const_data`.
#' @keywords internal
.normalize_bytecode_program <- function(bytecode = NULL, opcodes = NULL, const_data = numeric()) {
  code <- bytecode
  if (is.null(code) && !is.null(opcodes)) {
    code <- opcodes
  }
  if (is.null(code)) {
    code <- integer(0)
  }

  code <- as.integer(code)
  constants <- as.numeric(const_data %||% numeric(0))

  # Backward compatibility: some legacy saved transforms encoded unary maps like
  # softplus(x) as [SOFTPLUS] rather than [PUSH_X, SOFTPLUS]. Prepend PUSH_X so
  # replay in R and Stan stays stable for old fitted objects used in prediction.
  if (length(code) > 0L && code[[1]] %in% .bytecode_unary_ops()) {
    code <- c(0L, code)
  }

  list(
    bytecode = code,
    const_data = constants
  )
}

# Backward-compatible alias for older internal helper naming.
.emit_opcode_expr <- .emit_bytecode_expr
