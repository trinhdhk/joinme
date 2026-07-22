# Time-varying concordance for the survival component

Time-varying concordance for the survival component

## Usage

``` r
# S3 method for class 'JoiNMeFit'
concordance(
  object,
  newdataLong = NULL,
  newdataEvent = NULL,
  time_start,
  time_horizon = NULL,
  Dt = NULL,
  cause = 1,
  n_samples = 200,
  seed = 123,
  type_weights = "none",
  ...
)
```

## Arguments

- object:

  A fitted object of class `JoiNMeFit`.

- newdataLong:

  Longitudinal data for evaluation (defaults to training data).

- newdataEvent:

  Event data for evaluation (defaults to training data).

- time_start:

  Numeric landmark time(s) or a column name in `newdataEvent`.

- time_horizon:

  Numeric horizon time(s) or a column name in `newdataEvent`.

- Dt:

  Numeric horizon width; used when `time_horizon` is not supplied.

- cause:

  Integer; cause index for competing risks (default 1).

- n_samples:

  Number of posterior draws for prediction (default 200).

- seed:

  Random seed for draw subsetting.

- type_weights:

  Time-weighting for concordance; passed to
  [`survival::concordance`](https://rdrr.io/pkg/survival/man/concordance.html).

- ...:

  Additional arguments passed to
  [`predict.JoiNMeFit()`](https://trinhdhk.github.io/joinme/reference/predict.JoinMeFit.md).

## Value

A data frame with time-varying concordance at each landmark time.
