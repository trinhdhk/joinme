##' Functional Transform Builder for Joint Models
##'
##' Provides R utilities to parse user-friendly transformation expressions
##' and generate functional opcodes for Stan's arbitrary transformation evaluator.
##'
##' @description
##' The `parse_transform_expr` function takes a formula, quosure, quoted expression,
##' or character string and converts it to a functional opcode representation compatible with the
##' Stan-side opcode evaluator.
##'
##' Supported operations:
##' - Arithmetic: +, -, *, /
##' - Power: ^
##' - Functions: log, exp, sqrt, sin, cos, tan, abs, sinh, cosh, tanh, asinh, acosh, atanh
##' - Sigmoid/Link: inv_logit, logit
##' - Reciprocal: 1/x or rec(x)
##'
##' @param expr A formula (e.g. `~ x + 2`), quosure, quoted expression, or character string to parse.
##'
##' @note
##' Functional opcodes are a vector of operation codes (0-17) paired with a vector of
##' constant values. Constants are embedded using PUSH_CONST operations.
##'
##' @section Functional Opcode Reference:
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
##' }
##'
##' @examples
##' # Simple: f(x) = x + 2
##' bc <- parse_transform_expr(~ x + 2)
##'
##' # Use a plain formula (recommended)
##' bc <- parse_transform_expr(~ log(x + 1))
##' # Returns: opcodes = c(0, 1, 2), const_data = c(2)
##'
##' # Complex: f(x) = (log(sqrt(x + 1/inv_logit(3*x - 3))))^2
##' bc <- parse_transform_expr(~ (log(sqrt(x + 1/inv_logit(3*x - 3))))^2)
##'
##' @export
parse_transform_expr <- function(expr) {
  expr_call <- .coerce_transform_expr(expr)
  
  opcodes <- integer()
  const_data <- numeric()
  result <- .emit_opcode_expr(expr_call, opcodes, const_data)
  
  verify_opcodes(result$opcodes, result$const_data)
  
  list(
    opcodes = result$opcodes,
    const_data = result$const_data,
    n_ops = length(result$opcodes),
    n_const = length(result$const_data)
  )
}

##' @keywords internal
##' Coerce input to a language object using base R parsing.
.coerce_transform_expr <- function(expr) {
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
##' Emit opcodes from an R language object (calls, names, constants).
.emit_opcode_expr <- function(node, opcodes, const_data) {
  if (is.numeric(node)) {
    opcodes <- c(opcodes, 1L)
    const_data <- c(const_data, as.numeric(node))
    return(list(opcodes = opcodes, const_data = const_data))
  }
  if (is.name(node)) {
    if (as.character(node) != "x") {
      cli::cli_abort(c(
        x = "Unknown symbol in transform expression: {as.character(node)}.",
        i = "Use 'x' as the transformation input variable."
      ))
    }
    opcodes <- c(opcodes, 0L)
    return(list(opcodes = opcodes, const_data = const_data))
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
    return(.emit_opcode_expr(args[[1]], opcodes, const_data))
  }
  
  if (op %in% c("+", "-", "*", "/", "^")) {
    if (length(args) == 1 && op == "-") {
      res <- .emit_opcode_expr(0, opcodes, const_data)
      opcodes <- res$opcodes
      const_data <- res$const_data
      res <- .emit_opcode_expr(args[[1]], opcodes, const_data)
      opcodes <- res$opcodes
      const_data <- res$const_data
      opcodes <- c(opcodes, 3L)
      return(list(opcodes = opcodes, const_data = const_data))
    }
    if (length(args) != 2) {
      cli::cli_abort(c(
        x = "Binary operator {op} expects two arguments.",
        i = "Check the transform expression for missing operands."
      ))
    }
    res <- .emit_opcode_expr(args[[1]], opcodes, const_data)
    opcodes <- res$opcodes
    const_data <- res$const_data
    res <- .emit_opcode_expr(args[[2]], opcodes, const_data)
    opcodes <- res$opcodes
    const_data <- res$const_data
    opcode <- switch(op, "+" = 2L, "-" = 3L, "*" = 4L, "/" = 5L, "^" = 12L)
    opcodes <- c(opcodes, opcode)
    return(list(opcodes = opcodes, const_data = const_data))
  }
  
  func_opcode <- switch(op,
    log = 6L,
    exp = 7L,
    sqrt = 8L,
    inv_logit = 9L,
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
    NULL
  )
  if (!is.null(func_opcode)) {
    if (length(args) != 1) {
      cli::cli_abort(c(
        x = "Function {op} expects one argument.",
        i = "Check the transform expression for missing arguments."
      ))
    }
    res <- .emit_opcode_expr(args[[1]], opcodes, const_data)
    opcodes <- res$opcodes
    const_data <- res$const_data
    opcodes <- c(opcodes, func_opcode)
    return(list(opcodes = opcodes, const_data = const_data))
  }
  
  if (op == "ISpline") {
    cli::cli_abort(c(
      x = "ISpline() is not supported in functional opcode mode.",
      i = "Use the ispline transform mode and supply knots/coefficients instead."
    ))
  }
  
  cli::cli_abort(c(
    x = "Unsupported function in transform expression: {op}.",
    i = "Supported functions: log, exp, sqrt, inv_logit, logit, rec, sin, cos, tan, abs, sinh, cosh, tanh, asinh, acosh, atanh."
  ))
}

##' @keywords internal
##' Legacy tokenizer/parser removed in favor of R's formula parsing.

##' @keywords internal
##' Verify functional opcodes for basic sanity (stack under/overflow, etc.)
verify_opcodes <- function(opcodes, const_data) {
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
    } else if (op %in% c(6, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23)) {
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
