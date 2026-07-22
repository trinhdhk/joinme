# Stan diagnostics for JoiNMe models

Stan diagnostics for JoiNMe models

## Usage

``` r
stan_rhat.JoiNMeFit(
  object,
  pars = NULL,
  regex_pars = NULL,
  draws = NULL,
  seed = 1,
  ...
)

stan_ess.JoiNMeFit(
  object,
  pars = NULL,
  regex_pars = NULL,
  draws = NULL,
  seed = 1,
  type = c("bulk", "tail"),
  ...
)

stan_mcse.JoiNMeFit(
  object,
  pars = NULL,
  regex_pars = NULL,
  draws = NULL,
  seed = 1,
  type = c("mean", "sd", "median"),
  ...
)
```

## Arguments

- object:

  A fitted object of class `JoiNMeFit`.

- pars:

  Optional character vector of parameter names to include.

- regex_pars:

  Optional regular expression for parameter selection.

- draws:

  Optional number of posterior draws to subset.

- seed:

  Random seed for draw subsetting.

- ...:

  Unused.

- type:

  Diagnostic type (for ESS and MCSE).

## Value

A data frame of diagnostics by parameter.
