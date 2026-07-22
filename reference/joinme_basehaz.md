# Build specs for baseline hazard

This function provides a formal mean to specify the baseline hazard
specification for the survival submodel. It evaluates nothing, but it
checks the validity of the input and returns a structured object that
can be passed to
[`joinme()`](https://trinhdhk.github.io/joinme/reference/joinme.md).

## Usage

``` r
joinme_basehaz(
  type = c("bs", "ns", "formula"),
  n_knots = 5L,
  knots = NULL,
  degree = 3L,
  formula = ~1 + time
)

jm_basehaz(
  type = c("bs", "ns", "formula"),
  n_knots = 5L,
  knots = NULL,
  degree = 3L,
  formula = ~1 + time
)
```

## Arguments

- type:

  Baseline hazard type: "bs", "ns", or "formula".

- n_knots:

  Number of internal knots for spline baseline hazards.

- knots:

  Optional numeric vector of internal knots for spline baseline hazards.

- degree:

  Degree of spline basis for baseline hazard.

- formula:

  Formula for baseline hazard when `type = "formula"`. Ignored
  otherwise.

## Value

An object of class `joinme_basehaz`.

## Examples

``` r
basehaz_spec <- joinme_basehaz(
 type = "bs",
 n_knots = 5,
 degree = 3
)
```
