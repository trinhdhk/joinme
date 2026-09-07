# Invert a one-to-one transformation expression

Derives the algebraic inverse of a transformation written as a one-sided
formula in `x`. The result is another one-sided formula and can
therefore be passed directly to
[`parse_transform_expr()`](https://trinhdhk.github.io/joinme/reference/parse_transform_expr.md).
JoiNMe uses this function when `link` is supplied as a formula to
[`joinme_family()`](https://trinhdhk.github.io/joinme/reference/joinme_family.md).

Inversion proceeds from the outermost operation towards `x`. Composed
transformations are consequently reversed in the correct order. For
example, `~ log(2 * x + 1)` becomes `~ (exp(x) - 1) / 2`.

The supported one-to-one elementary pairs are `log`/`exp`,
`logit`/`inv_logit`, `Phi`/`inv_Phi` (with `pnorm`/`qnorm` and `probit`
aliases), `sinh`/`asinh`, `tanh`/`atanh`, reciprocal, `sqrt`/square, and
`cbrt`/cube. Addition, subtraction, multiplication, division, and powers
by a finite non-zero constant are also supported. Arithmetic branches
that do not contain `x` must reduce to numeric constants. Even integer
powers are rejected because they are not one-to-one on the real line; a
deliberately restricted-domain inverse can instead be declared
explicitly through `inv_link`.

Expressions containing `x` more than once, or globally non-injective
operations such as `abs`, `sin`, `cos`, and `cosh`, are rejected. This
is deliberate: silently selecting one branch would not define a valid
inverse-link over the full response domain.

## Usage

``` r
invert_transform_expr(expr)
```

## Arguments

- expr:

  A one-sided formula, quosure, quoted expression, or character string
  containing exactly one occurrence of `x`.

## Value

A one-sided formula whose right-hand side is the algebraic inverse
transformation, expressed in `x`.

## Examples

``` r
invert_transform_expr(~ log(x))
#> ~exp(x)
#> <environment: 0x5624dc825148>
invert_transform_expr(~ log(2 * x + 1))
#> ~(exp(x) - 1)/2
#> <environment: 0x5624dc825148>
invert_transform_expr(~ x^3)
#> ~cbrt(x)
#> <environment: 0x5624dc825148>
```
