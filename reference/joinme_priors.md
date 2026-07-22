# Declare priors with validation

`joinme_priors()` is a user-facing wrapper for prior declarations passed
to [`joinme()`](https://trinhdhk.github.io/joinme/reference/joinme.md).
It validates the exposed prior components and returns a structured
object that can be supplied directly as the `priors` argument.

Exposed components are:

- `beta`: longitudinal/survival regression prior scale(s)

- `alpha`: association prior scale

- `iota`: fit-only affine-shift prior scale for functional association
  transforms

- `lkj`: LKJ concentration parameter

## Usage

``` r
joinme_priors(
  beta = NULL,
  alpha = NULL,
  iota = NULL,
  lkj = NULL,
  .validate = TRUE
)

jm_priors(beta = NULL, alpha = NULL, iota = NULL, lkj = NULL, .validate = TRUE)
```

## Arguments

- beta:

  Beta prior specification. Use either a numeric scale (or vector of
  scales) or a list with `scale`/`sd`.

- alpha:

  Alpha prior specification. Use either a numeric scale or a list with
  `scale`/`sd`.

- iota:

  Iota prior specification for fit-only affine-shift intercept and slope
  parameters in functional association transforms. Use either a numeric
  scale or a list with `scale`/`sd`.

- lkj:

  LKJ concentration parameter.

- .validate:

  Logical; if `TRUE` (default), validate the resulting prior
  declarations immediately.

## Value

An object of class `joinme_priors`.

## Examples

``` r
pri <- joinme_priors(
  beta = list(scale = 2.5),
  alpha = list(scale = 1.0),
  iota = list(scale = 1.0),
  lkj = 2
)
print(pri)
#> Prior specification for Joint Mixed Effects model
#>  component value
#>       beta   2.5
#>      alpha     1
#>       iota     1
#>        lkj     2
```
