# Update a JoiNMe fit

Refits a JoiNMe model using the stored call, with optional updates to
formulas, data, and control arguments.

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
  .env = NULL,
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

  Optional updated covariance regression. Supply one formula to update
  both components, or `list(sd = ~ ..., corr = ~ ...)` to update the
  standard-deviation and off-diagonal correlation regressions
  independently. Update formulae containing `.` are supported within
  each component.

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

- .env:

  Optional environment for evaluating the updated call. If NULL, the
  parent frame is used. When failed, the environment of the original
  formulaLong is used.

- ...:

  Additional arguments passed to
  [`joinme()`](https://trinhdhk.github.io/joinme/reference/joinme.md),
  or to
  [`joinme_mix()`](https://trinhdhk.github.io/joinme/reference/joinme_mix.md)
  when `object` is a latent-progress mixture.

## Value

A refitted JoiNMe object. Mixture fits retain their mixture entry point
and class specification. Any longitudinal-only fit remains longitudinal
only unless new event inputs are supplied explicitly.

## See also

[`update()`](https://rdrr.io/r/stats/update.html)
