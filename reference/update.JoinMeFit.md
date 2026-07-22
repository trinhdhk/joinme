# Update a JoiNMe fit

Refits a JoiNMe model using the stored call, with optional updates to
formulas, data, and control arguments. This mirrors the pattern used by
brms::update.

## Usage

``` r
# S3 method for class 'JoiNMeFit'
update(
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
  ...
)
```

## Arguments

- object:

  A JoiNMe fit object.

- formulaLong:

  Optional updated longitudinal formula. Use an update formula (e.g.,
  `~ . + x`) to modify the existing model.

- dataLong:

  Optional updated longitudinal dataset.

- formulaEvent:

  Optional updated survival formula (full or update form).

- dataEvent:

  Optional updated event dataset.

- formulaVCov:

  Optional updated covariance formula (full or update form).

- formulaDist:

  Optional distributional regression formulas with parameter names on
  the LHS (e.g., `sigma ~ 1 + time`).

- control:

  Optional updated control list.

- draws:

  Optional draws override.

- families:

  Optional updated families specification.

- transforms:

  Optional updated transforms specification.

- priors:

  Optional updated priors list.

- ...:

  Additional arguments passed to
  [`joinme()`](https://trinhdhk.github.io/joinme/reference/joinme.md).

## Value

A refitted JoiNMe object.
