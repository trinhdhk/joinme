#' Functional Transform Builder for Joint Models
#'
#' Provides R utilities to parse user-friendly transformation expressions
#' and generate functional bytecode for Stan's arbitrary transformation evaluator.
#'
#' @description
#' The `parse_transform_expr` function takes a formula, quosure, quoted expression,
#' or character string and converts it to a functional bytecode representation compatible with the
#' Stan-side bytecode evaluator.
#'
#' Supported operations:
#' - Arithmetic: +, -, *, /
#' - Power: ^, power
#' - Root: sqrt, cbrt
#' - Exponential transformations: log, exp, 
#' - Trigonometric functions: sin, cos, tan, abs, sinh, cosh, tanh, asinh, acosh, atanh
#' - Common links: inv_logit (softmax, SoftMax), logit, sigmoid, expit,
#'   softplus (log1p_exp), `Phi`/`pnorm`, and
#'   `inv_Phi`/`qnorm`/`probit`
#' - Reciprocal: 1/x or rec(x)
#'
#' @param expr A formula (e.g. `~ x + 2`), quosure, quoted expression, or character string to parse.
#' @param iota_nodes Optional internal description of transformation-specific
#'   affine shifts. Each declared node identifies whether an intercept, slope,
#'   or both are estimated for one nonlinear operation.
#'
#' @note
#' Functional bytecode is a vector of operation codes (0-27) paired with a vector of
#' constant values. Constants are embedded using PUSH_CONST operations.
#' `Phi` denotes the standard normal cumulative distribution function, whereas
#' `inv_Phi` and `probit` denote its quantile function. They are deliberately
#' assigned different instructions because their domains and ranges differ.
#'
#' @section Functional Bytecode Reference:
#' \describe{
#'   \item{0}{PUSH_X: Push input x onto stack}
#'   \item{1}{PUSH_CONST: Push next constant value}
#'   \item{2}{ADD: Pop b,a; push a+b}
#'   \item{3}{SUB: Pop b,a; push a-b}
#'   \item{4}{MUL: Pop b,a; push a*b}
#'   \item{5}{DIV: Pop b,a; push a/b}
#'   \item{6}{LOG: Pop a; push log(a)}
#'   \item{7}{EXP: Pop a; push exp(a)}
#'   \item{8}{SQRT: Pop a; push sqrt(a)}
#'   \item{9}{INV_LOGIT: Pop a; push inv_logit(a)}
#'   \item{10}{LOGIT: Pop a; push logit(a)}
#'   \item{11}{RECIPROCAL: Pop a; push 1/a}
#'   \item{12}{POW: Pop b,a; push a^b}
#'   \item{13}{SIN: Pop a; push sin(a)}
#'   \item{14}{COS: Pop a; push cos(a)}
#'   \item{15}{TAN: Pop a; push tan(a)}
#'   \item{16}{ABS: Pop a; push abs(a)}
#'   \item{17}{SQUARE: Pop a; push a^2}
#'   \item{18}{SINH: Pop a; push sinh(a)}
#'   \item{19}{COSH: Pop a; push cosh(a)}
#'   \item{20}{TANH: Pop a; push tanh(a)}
#'   \item{21}{ASINH: Pop a; push asinh(a)}
#'   \item{22}{ACOSH: Pop a; push acosh(a)}
#'   \item{23}{ATANH: Pop a; push atanh(a)}
#'   \item{24}{SOFTPLUS: Pop a; push log1p_exp(a)}
#'   \item{25}{CBRT: Pop a; push cbrt(a)}
#'   \item{26}{PHI: Pop a; push the standard normal CDF \eqn{\Phi(a)}}
#'   \item{27}{INV_PHI: Pop a; push the standard normal quantile
#'     \eqn{\Phi^{-1}(a)}}
#' }
#'
#' @examples
#' # Simple: f(x) = x + 2
#' bc <- parse_transform_expr(~ x + 2)
#'
#' # Use a plain formula (recommended)
#' bc <- parse_transform_expr(~ log(x + 1))
#' # Returns: bytecode = c(0, 1, 2), const_data = c(2)
#'
#' # Complex: f(x) = (log(sqrt(x + 1/inv_logit(3*x - 3))))^2
#' bc <- parse_transform_expr(~ (log(sqrt(x + 1/inv_logit(3*x - 3))))^2)
#'
#' # The normal CDF and its quantile use separate instructions
#' phi_bc <- parse_transform_expr(~ Phi(x))
#' probit_bc <- parse_transform_expr(~ inv_Phi(x))
#'
#' @export
# 
# - Parse R expressions into portable bytecode programs.
# - Validate bytecode stack behaviour before handing programs to Stan.
# - Evaluate bytecode in R (shared by inverse-link and association transforms).
parse_transform_expr <- function(expr, iota_nodes = NULL) {
  # Step 1: normalise any supported expression input into a language object.
  expr_call <- .coerce_transform_expr(expr)
  iota_node_map <- .iota_node_map(iota_nodes)

  # Step 2: emit bytecode instructions and constant vector.
  bytecode <- integer()
  const_data <- numeric()
  op_iota_intercept_idx <- integer()
  op_iota_slope_idx <- integer()
  result <- .emit_bytecode_expr(
    expr_call,
    bytecode,
    const_data,
    op_iota_intercept_idx,
    op_iota_slope_idx,
    iota_node_map = iota_node_map,
    path = "root"
  )

  # Step 3: validate stack consistency so runtime evaluation is deterministic.
  verify_bytecode(result$bytecode, result$const_data)

  # Step 4: return canonical bytecode payload.
  list(
    bytecode = result$bytecode,
    const_data = result$const_data,
    op_iota_intercept_idx = result$op_iota_intercept_idx,
    op_iota_slope_idx = result$op_iota_slope_idx,
    n_ops = length(result$bytecode),
    n_bytecode = length(result$bytecode),
    n_const = length(result$const_data),
    n_iota_intercept = max(c(0L, result$op_iota_intercept_idx)),
    n_iota_slope = max(c(0L, result$op_iota_slope_idx))
  )
}

#' Invert a one-to-one transformation expression
#'
#' @description
#' Derives the algebraic inverse of a transformation written as a one-sided
#' formula in `x`.  The result is another one-sided formula and can therefore
#' be passed directly to [parse_transform_expr()].  JoiNMe uses this function
#' when `link` is supplied as a formula to [joinme_family()].
#'
#' Inversion proceeds from the outermost operation towards `x`.  Composed
#' transformations are consequently reversed in the correct order.  For
#' example, `~ log(2 * x + 1)` becomes `~ (exp(x) - 1) / 2`.
#'
#' The supported one-to-one elementary pairs are `log`/`exp`,
#' `logit`/`inv_logit`, `Phi`/`inv_Phi` (with `pnorm`/`qnorm` and
#' `probit` aliases), `sinh`/`asinh`,
#' `tanh`/`atanh`, reciprocal, `sqrt`/square, and `cbrt`/cube.  Addition,
#' subtraction, multiplication, division, and powers by a finite non-zero
#' constant are also supported.  Arithmetic branches that do not contain `x`
#' must reduce to numeric constants.  Even integer powers are rejected because
#' they are not one-to-one on the real line; a deliberately restricted-domain
#' inverse can instead be declared explicitly through `inv_link`.
#'
#' Expressions containing `x` more than once, or globally non-injective
#' operations such as `abs`, `sin`, `cos`, and `cosh`, are rejected.  This is
#' deliberate: silently selecting one branch would not define a valid
#' inverse-link over the full response domain.
#'
#' @param expr A one-sided formula, quosure, quoted expression, or character
#'   string containing exactly one occurrence of `x`.
#'
#' @return A one-sided formula whose right-hand side is the algebraic inverse
#'   transformation, expressed in `x`.
#' @examples
#' invert_transform_expr(~ log(x))
#' invert_transform_expr(~ log(2 * x + 1))
#' invert_transform_expr(~ x^3)
#' @export
invert_transform_expr <- function(expr) {
  node <- .coerce_transform_expr(expr)
  n_x <- .transform_x_count(node)
  if (n_x != 1L) {
    cli::cli_abort(c(
      x = "A formula link must contain {.code x} exactly once.",
      i = "The supplied expression contains {n_x} occurrences."
    ))
  }

  inverse <- .invert_transform_node(node, quote(x))
  environment <- if (inherits(expr, "formula")) environment(expr) else parent.frame()
  rlang::new_formula(lhs = NULL, rhs = inverse, env = environment)
}

#' Count occurrences of the transformation variable
#'
#' @param node An R language object forming part of a transformation
#'   expression.
#'
#' @return A non-negative integer count of symbols named `x`.
#' @keywords internal
#' @noRd
.transform_x_count <- function(node) {
  if (is.name(node)) {
    return(as.integer(identical(as.character(node), "x")))
  }
  if (!is.call(node)) {
    return(0L)
  }
  sum(vapply(as.list(node)[-1L], .transform_x_count, integer(1)))
}

#' Reduce an `x`-free arithmetic branch to a numeric constant
#'
#' @description
#' Evaluates only literal numeric arithmetic needed while solving a formula
#' link.  Arbitrary R evaluation is intentionally excluded so inversion is
#' deterministic and does not depend on objects in a caller's environment.
#'
#' @param node An `x`-free arithmetic expression.
#'
#' @return One finite numeric scalar.
#' @keywords internal
#' @noRd
.transform_constant_value <- function(node) {
  if (is.numeric(node) && length(node) == 1L && is.finite(node)) {
    return(as.numeric(node))
  }
  if (!is.call(node)) {
    cli::cli_abort("Non-numeric constants are not supported in a formula link.")
  }

  operator <- as.character(node[[1L]])
  arguments <- as.list(node)[-1L]
  if (operator %in% c("+", "-") && length(arguments) == 1L) {
    value <- .transform_constant_value(arguments[[1L]])
    return(if (operator == "-") -value else value)
  }
  if (!(operator %in% c("+", "-", "*", "/", "^")) ||
      length(arguments) != 2L) {
    cli::cli_abort("Constant branches in a formula link may use only numeric arithmetic.")
  }

  left <- .transform_constant_value(arguments[[1L]])
  right <- .transform_constant_value(arguments[[2L]])
  value <- switch(
    operator,
    "+" = left + right,
    "-" = left - right,
    "*" = left * right,
    "/" = left / right,
    "^" = left^right
  )
  if (length(value) != 1L || !is.finite(value)) {
    cli::cli_abort("A constant branch in the formula link does not have a finite value.")
  }
  as.numeric(value)
}

#' Construct a binary transformation call
#'
#' @param operator Character arithmetic operator.
#' @param left,right Language objects or numeric constants.
#'
#' @return An unevaluated R call.
#' @keywords internal
#' @noRd
.transform_binary_call <- function(operator, left, right) {
  as.call(list(as.name(operator), left, right))
}

#' Solve one transformation expression for its input
#'
#' @description
#' Recursively solves `node = value` for the unique occurrence of `x`.
#' Each recursion removes one outer operation and passes the transformed
#' right-hand side inwards.  This produces the reverse composition required
#' for an inverse link without evaluating user-supplied code.
#'
#' @param node Current branch of the forward transformation.
#' @param value Language object representing the current right-hand side.
#'
#' @return A language object expressing the original `x` in terms of `value`.
#' @keywords internal
#' @noRd
.invert_transform_node <- function(node, value) {
  if (is.name(node) && identical(as.character(node), "x")) {
    return(value)
  }
  if (!is.call(node)) {
    cli::cli_abort("Could not isolate {.code x} in the formula link.")
  }

  operator <- tolower(as.character(node[[1L]]))
  arguments <- as.list(node)[-1L]

  if (operator == "(" && length(arguments) == 1L) {
    return(.invert_transform_node(arguments[[1L]], value))
  }
  if (operator %in% c("+", "-") && length(arguments) == 1L) {
    next_value <- if (operator == "-") {
      .transform_binary_call("-", 0, value)
    } else {
      value
    }
    return(.invert_transform_node(arguments[[1L]], next_value))
  }

  inverse_function <- switch(
    operator,
    log = "exp",
    exp = "log",
    logit = "inv_logit",
    inv_logit = "logit",
    sigmoid = "logit",
    expit = "logit",
    softmax = "logit",
    phi = "inv_Phi",
    pnorm = "qnorm",
    inv_phi = "Phi",
    qnorm = "pnorm",
    probit = "Phi",
    sinh = "asinh",
    asinh = "sinh",
    tanh = "atanh",
    atanh = "tanh",
    rec = "rec",
    NULL
  )
  if (!is.null(inverse_function)) {
    if (length(arguments) != 1L) {
      cli::cli_abort("Function {.fn {operator}} must have one argument in a formula link.")
    }
    return(.invert_transform_node(
      arguments[[1L]],
      as.call(list(as.name(inverse_function), value))
    ))
  }

  if (operator %in% c("sqrt", "cbrt", "softplus", "log1p_exp")) {
    if (length(arguments) != 1L) {
      cli::cli_abort("Function {.fn {operator}} must have one argument in a formula link.")
    }
    next_value <- switch(
      operator,
      sqrt = .transform_binary_call("^", value, 2),
      cbrt = .transform_binary_call("^", value, 3),
      softplus = as.call(list(
        as.name("log"),
        .transform_binary_call("-", as.call(list(as.name("exp"), value)), 1)
      )),
      log1p_exp = as.call(list(
        as.name("log"),
        .transform_binary_call("-", as.call(list(as.name("exp"), value)), 1)
      ))
    )
    return(.invert_transform_node(arguments[[1L]], next_value))
  }

  if (operator == "power") {
    operator <- "^"
  }
  if (operator %in% c("+", "-", "*", "/", "^") &&
      length(arguments) == 2L) {
    left_count <- .transform_x_count(arguments[[1L]])
    right_count <- .transform_x_count(arguments[[2L]])
    if (left_count + right_count != 1L) {
      cli::cli_abort("Each arithmetic operation in a formula link must have exactly one branch containing {.code x}.")
    }

    if (left_count == 1L) {
      constant <- .transform_constant_value(arguments[[2L]])
      next_value <- switch(
        operator,
        "+" = .transform_binary_call("-", value, constant),
        "-" = .transform_binary_call("+", value, constant),
        "*" = {
          if (constant == 0) cli::cli_abort("Multiplication by zero is not invertible.")
          .transform_binary_call("/", value, constant)
        },
        "/" = {
          if (constant == 0) cli::cli_abort("Division by zero is not a valid link.")
          .transform_binary_call("*", value, constant)
        },
        "^" = {
          if (constant == 0) cli::cli_abort("A zero power is not invertible.")
          if (constant > 0 && constant == round(constant) &&
              as.integer(constant) %% 2L == 0L) {
            cli::cli_abort(c(
              x = "An even integer power is not one-to-one on the real line.",
              i = "Supply {.arg inv_link} directly only when a scientifically justified restricted domain is intended."
            ))
          }
          if (constant == 3) {
            as.call(list(as.name("cbrt"), value))
          } else {
            .transform_binary_call("^", value, 1 / constant)
          }
        }
      )
      return(.invert_transform_node(arguments[[1L]], next_value))
    }

    constant <- .transform_constant_value(arguments[[1L]])
    next_value <- switch(
      operator,
      "+" = .transform_binary_call("-", value, constant),
      "-" = .transform_binary_call("-", constant, value),
      "*" = {
        if (constant == 0) cli::cli_abort("Multiplication by zero is not invertible.")
        .transform_binary_call("/", value, constant)
      },
      "/" = {
        if (constant == 0) cli::cli_abort("A zero numerator is not invertible.")
        .transform_binary_call("/", constant, value)
      },
      "^" = {
        if (constant <= 0 || constant == 1) {
          cli::cli_abort("A constant power base must be positive and different from one.")
        }
        .transform_binary_call(
          "/",
          as.call(list(as.name("log"), value)),
          log(constant)
        )
      }
    )
    return(.invert_transform_node(arguments[[2L]], next_value))
  }

  cli::cli_abort(c(
    x = "Function or operation {.fn {operator}} is not invertible.",
    i = "Use a one-to-one composition of supported arithmetic and link functions, or supply {.arg inv_link} directly."
  ))
}

#' @keywords internal
.iota_node_map <- function(iota_nodes = NULL) {
  if (is.null(iota_nodes) || !length(iota_nodes)) {
    return(list())
  }
  stats::setNames(iota_nodes, vapply(iota_nodes, function(node) node$path %||% "", character(1)))
}

#' @keywords internal
.bytecode_iota_indices <- function(iota_node_map, path, function_name) {
  node <- iota_node_map[[path]]
  if (is.null(node)) {
    return(list(intercept = 0L, slope = 0L))
  }
  node_fun <- tolower(as.character(node$function_name %||% ""))
  if (!identical(node_fun, tolower(function_name))) {
    cli::cli_abort(c(
      x = "Transform affine-shift metadata is inconsistent with the functional expression.",
      i = "Rebuild the transform specification before parsing bytecode."
    ))
  }
  list(
    intercept = as.integer(node$intercept_index %||% 0L),
    slope = as.integer(node$slope_index %||% 0L)
  )
}

#' Coerce input to a language object
#' @keywords internal
#' @noRd
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

#' Emit bytecode from an R language object (calls, names, constants).
#' @keywords internal
#' @noRd
.emit_bytecode_expr <- function(
  node,
  bytecode,
  const_data,
  op_iota_intercept_idx,
  op_iota_slope_idx,
  iota_node_map = list(),
  path = "root"
) {
  # Recursive descent over AST nodes to emit bytecode
  if (is.numeric(node)) {
    bytecode <- c(bytecode, 1L)
    const_data <- c(const_data, as.numeric(node))
    op_iota_intercept_idx <- c(op_iota_intercept_idx, 0L)
    op_iota_slope_idx <- c(op_iota_slope_idx, 0L)
    return(list(
      bytecode = bytecode,
      const_data = const_data,
      op_iota_intercept_idx = op_iota_intercept_idx,
      op_iota_slope_idx = op_iota_slope_idx
    ))
  }
  if (is.name(node)) {
    if (as.character(node) != "x") {
      cli::cli_abort(c(
        x = "Unknown symbol in transform expression: {as.character(node)}.",
        i = "Use 'x' as the transformation input variable."
      ))
    }
    bytecode <- c(bytecode, 0L)
    op_iota_intercept_idx <- c(op_iota_intercept_idx, 0L)
    op_iota_slope_idx <- c(op_iota_slope_idx, 0L)
    return(list(
      bytecode = bytecode,
      const_data = const_data,
      op_iota_intercept_idx = op_iota_intercept_idx,
      op_iota_slope_idx = op_iota_slope_idx
    ))
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
    return(.emit_bytecode_expr(
      args[[1]], bytecode, const_data,
      op_iota_intercept_idx, op_iota_slope_idx,
      iota_node_map = iota_node_map,
      path = paste0(path, "/1")
    ))
  }
  
  if (op %in% c("+", "-", "*", "/", "^")) {
    if (length(args) == 1 && op == "-") {
      res <- .emit_bytecode_expr(0, bytecode, const_data, op_iota_intercept_idx, op_iota_slope_idx, iota_node_map = iota_node_map, path = paste0(path, "/1"))
      bytecode <- res$bytecode
      const_data <- res$const_data
      op_iota_intercept_idx <- res$op_iota_intercept_idx
      op_iota_slope_idx <- res$op_iota_slope_idx
      res <- .emit_bytecode_expr(args[[1]], bytecode, const_data, op_iota_intercept_idx, op_iota_slope_idx, iota_node_map = iota_node_map, path = paste0(path, "/1"))
      bytecode <- res$bytecode
      const_data <- res$const_data
      op_iota_intercept_idx <- res$op_iota_intercept_idx
      op_iota_slope_idx <- res$op_iota_slope_idx
      bytecode <- c(bytecode, 3L)
      op_iota_intercept_idx <- c(op_iota_intercept_idx, 0L)
      op_iota_slope_idx <- c(op_iota_slope_idx, 0L)
      return(list(
        bytecode = bytecode,
        const_data = const_data,
        op_iota_intercept_idx = op_iota_intercept_idx,
        op_iota_slope_idx = op_iota_slope_idx
      ))
    }
    if (length(args) != 2) {
      cli::cli_abort(c(
        x = "Binary operator {op} expects two arguments.",
        i = "Check the transform expression for missing operands."
      ))
    }
    res <- .emit_bytecode_expr(args[[1]], bytecode, const_data, op_iota_intercept_idx, op_iota_slope_idx, iota_node_map = iota_node_map, path = paste0(path, "/1"))
    bytecode <- res$bytecode
    const_data <- res$const_data
    op_iota_intercept_idx <- res$op_iota_intercept_idx
    op_iota_slope_idx <- res$op_iota_slope_idx
    res <- .emit_bytecode_expr(args[[2]], bytecode, const_data, op_iota_intercept_idx, op_iota_slope_idx, iota_node_map = iota_node_map, path = paste0(path, "/2"))
    bytecode <- res$bytecode
    const_data <- res$const_data
    op_iota_intercept_idx <- res$op_iota_intercept_idx
    op_iota_slope_idx <- res$op_iota_slope_idx
    op_id <- switch(op, "+" = 2L, "-" = 3L, "*" = 4L, "/" = 5L, "^" = 12L)
    bytecode <- c(bytecode, op_id)
    op_iota_intercept_idx <- c(op_iota_intercept_idx, 0L)
    op_iota_slope_idx <- c(op_iota_slope_idx, 0L)
    return(list(
      bytecode = bytecode,
      const_data = const_data,
      op_iota_intercept_idx = op_iota_intercept_idx,
      op_iota_slope_idx = op_iota_slope_idx
    ))
  }

  if (op == "power") {
    if (length(args) != 2) {
      cli::cli_abort(c(
        x = "Function {op} expects two arguments.",
        i = "Check the transform expression for missing arguments."
      ))
    }
    res <- .emit_bytecode_expr(args[[1]], bytecode, const_data, op_iota_intercept_idx, op_iota_slope_idx, iota_node_map = iota_node_map, path = paste0(path, "/1"))
    bytecode <- res$bytecode
    const_data <- res$const_data
    op_iota_intercept_idx <- res$op_iota_intercept_idx
    op_iota_slope_idx <- res$op_iota_slope_idx
    res <- .emit_bytecode_expr(args[[2]], bytecode, const_data, op_iota_intercept_idx, op_iota_slope_idx, iota_node_map = iota_node_map, path = paste0(path, "/2"))
    bytecode <- res$bytecode
    const_data <- res$const_data
    op_iota_intercept_idx <- res$op_iota_intercept_idx
    op_iota_slope_idx <- res$op_iota_slope_idx
    bytecode <- c(bytecode, 12L)
    op_iota <- .bytecode_iota_indices(iota_node_map, path, "power")
    op_iota_intercept_idx <- c(op_iota_intercept_idx, op_iota$intercept)
    op_iota_slope_idx <- c(op_iota_slope_idx, op_iota$slope)
    return(list(
      bytecode = bytecode,
      const_data = const_data,
      op_iota_intercept_idx = op_iota_intercept_idx,
      op_iota_slope_idx = op_iota_slope_idx
    ))
  }
  
  op <- tolower(op)

  normal_ops <- .bytecode_normal_ops()
  func_op_id <- switch(op,
    log = 6L,
    exp = 7L,
    sqrt = 8L,
    inv_logit = 9L,
    sigmoid = 9L,
    expit = 9L,
    softmax = 9L,
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
    phi = unname(normal_ops[["PHI"]]),
    pnorm = unname(normal_ops[["PHI"]]),
    inv_phi = unname(normal_ops[["INV_PHI"]]),
    qnorm = unname(normal_ops[["INV_PHI"]]),
    probit = unname(normal_ops[["INV_PHI"]]),
    NULL
  )
  if (!is.null(func_op_id)) {
    if (length(args) != 1) {
      cli::cli_abort(c(
        x = "Function {op} expects one argument.",
        i = "Check the transform expression for missing arguments."
      ))
    }
    res <- .emit_bytecode_expr(args[[1]], bytecode, const_data, op_iota_intercept_idx, op_iota_slope_idx, iota_node_map = iota_node_map, path = paste0(path, "/1"))
    bytecode <- res$bytecode
    const_data <- res$const_data
    op_iota_intercept_idx <- res$op_iota_intercept_idx
    op_iota_slope_idx <- res$op_iota_slope_idx
    bytecode <- c(bytecode, func_op_id)
    op_iota <- .bytecode_iota_indices(iota_node_map, path, op)
    op_iota_intercept_idx <- c(op_iota_intercept_idx, op_iota$intercept)
    op_iota_slope_idx <- c(op_iota_slope_idx, op_iota$slope)
    return(list(
      bytecode = bytecode,
      const_data = const_data,
      op_iota_intercept_idx = op_iota_intercept_idx,
      op_iota_slope_idx = op_iota_slope_idx
    ))
  }
  
  if (op == "ispline") {
    cli::cli_abort(c(
      x = "ISpline() is not supported in functional bytecode mode.",
      i = "Use the ispline transform mode and supply knots/coefficients instead."
    ))
  }
  
  cli::cli_abort(c(
    x = "Unsupported function in transform expression: {op}.",
    i = "Supported functions: log, exp, sqrt, inv_logit, logit, Phi, pnorm, inv_Phi, qnorm, probit, sigmoid, expit, softmax, softplus, log1p_exp, cbrt, power, rec, sin, cos, tan, abs, sinh, cosh, tanh, asinh, acosh, atanh."
  ))
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
  c(PHI = 26L, INV_PHI = 27L)
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
#' legacy one-operation program normalisation.
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
      stack[length(stack)] <- softplus(resolve_iota_shift(stack[length(stack)], op_intercept_idx, op_slope_idx))
    } else if (op == 25L) {
      a <- resolve_iota_shift(stack[length(stack)], op_intercept_idx, op_slope_idx)
      stack[length(stack)] <- sign(a) * abs(a)^(1 / 3)
    } else if (op == normal_ops[["PHI"]]) {
      stack[length(stack)] <- stats::pnorm(resolve_iota_shift(stack[length(stack)], op_intercept_idx, op_slope_idx))
    } else if (op == normal_ops[["INV_PHI"]]) {
      stack[length(stack)] <- stats::qnorm(resolve_iota_shift(stack[length(stack)], op_intercept_idx, op_slope_idx))
    } else {
      cli::cli_abort("Unknown transform bytecode instruction: {op}.")
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
  constants <- as.numeric(const_data %||% numeric(0))
  iota_intercepts <- as.integer(op_iota_intercept_idx %||% rep(0L, length(code)))
  iota_slopes <- as.integer(op_iota_slope_idx %||% rep(0L, length(code)))

  # Backward compatibility: some legacy saved transforms encoded unary maps like
  # softplus(x) as [SOFTPLUS] rather than [PUSH_X, SOFTPLUS]. Prepend PUSH_X so
  # replay in R and Stan stays stable for old fitted objects used in prediction.
  if (length(code) > 0L && code[[1]] %in% .bytecode_unary_ops()) {
    code <- c(0L, code)
    iota_intercepts <- c(0L, iota_intercepts)
    iota_slopes <- c(0L, iota_slopes)
  }

  if (length(iota_intercepts) != length(code) || length(iota_slopes) != length(code)) {
    cli::cli_abort(c(
      x = "Functional transform bytecode metadata is malformed.",
      i = "The iota index arrays must align with the bytecode length."
    ))
  }

  list(
    bytecode = code,
    const_data = constants,
    op_iota_intercept_idx = iota_intercepts,
    op_iota_slope_idx = iota_slopes
  )
}
