# Functional Transform Builder for Joint Models

The `parse_transform_expr` function takes a formula, quosure, quoted
expression, or character string and converts it to a functional bytecode
representation compatible with the Stan-side bytecode evaluator.

Supported operations:

- Arithmetic: +, -, \*, /

- Power: ^, power

- Root: sqrt, cbrt

- Exponential transformations: log, exp,

- Trigonometric functions: sin, cos, tan, abs, sinh, cosh, tanh, asinh,
  acosh, atanh

- Common links: inv_logit (softmax, SoftMax), logit, sigmoid, expit,
  softplus (log1p_exp)

- Reciprocal: 1/x or rec(x)

## Usage

``` r
parse_transform_expr(expr, iota_nodes = NULL)
```

## Arguments

- expr:

  A formula (e.g. `~ x + 2`), quosure, quoted expression, or character
  string to parse.

## Details

Provides R utilities to parse user-friendly transformation expressions
and generate functional bytecode for Stan's arbitrary transformation
evaluator.

## Note

Functional bytecode is a vector of operation codes (0-25) paired with a
vector of constant values. Constants are embedded using PUSH_CONST
operations.

## Functional Bytecode Reference

- 0:

  PUSH_X: Push input x onto stack

- 1:

  PUSH_CONST: Push next constant value

- 2:

  ADD: Pop b,a; push a+b

- 3:

  SUB: Pop b,a; push a-b

- 4:

  MUL: Pop b,a; push a\*b

- 5:

  DIV: Pop b,a; push a/b

- 6:

  LOG: Pop a; push log(a)

- 7:

  EXP: Pop a; push exp(a)

- 8:

  SQRT: Pop a; push sqrt(a)

- 9:

  INV_LOGIT: Pop a; push inv_logit(a)

- 10:

  LOGIT: Pop a; push logit(a)

- 11:

  RECIPROCAL: Pop a; push 1/a

- 12:

  POW: Pop b,a; push a^b

- 13:

  SIN: Pop a; push sin(a)

- 14:

  COS: Pop a; push cos(a)

- 15:

  TAN: Pop a; push tan(a)

- 16:

  ABS: Pop a; push abs(a)

- 17:

  SQUARE: Pop a; push a^2

- 18:

  SINH: Pop a; push sinh(a)

- 19:

  COSH: Pop a; push cosh(a)

- 20:

  TANH: Pop a; push tanh(a)

- 21:

  ASINH: Pop a; push asinh(a)

- 22:

  ACOSH: Pop a; push acosh(a)

- 23:

  ATANH: Pop a; push atanh(a)

- 24:

  SOFTPLUS: Pop a; push log1p_exp(a)

- 25:

  CBRT: Pop a; push cbrt(a)

- 26:

  PROBIT: Pop a; push Phi(a)

## Examples

``` r
# Simple: f(x) = x + 2
bc <- parse_transform_expr(~ x + 2)

# Use a plain formula (recommended)
bc <- parse_transform_expr(~ log(x + 1))
# Returns: bytecode = c(0, 1, 2), const_data = c(2)

# Complex: f(x) = (log(sqrt(x + 1/inv_logit(3*x - 3))))^2
bc <- parse_transform_expr(~ (log(sqrt(x + 1/inv_logit(3*x - 3))))^2)
```
