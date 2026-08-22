# Summarise a latent-progress mixture fit

Extends the ordinary JoiNMe summary with baseline class probabilities,
selected component locations and scales, class-membership regression,
and a compact posterior allocation-count distribution. The allocation
table contains one row per class and active allocation domain; use
[`posterior_class()`](https://trinhdhk.github.io/joinme/reference/posterior_class.md)
when subject-by-class or marker-by-class probabilities are required.

## Usage

``` r
# S3 method for class 'JoiNMeMixFit'
summary(
  object,
  draws = NULL,
  seed = .Random.seed[[1]],
  digits = 3,
  include_corr = TRUE,
  ...
)
```

## Arguments

- object:

  A JoiNMe fit object.

- draws:

  Number of draws to use for summaries.

- seed:

  Random seed for subsetting draws.

- digits:

  Number of digits to round summary values.

- include_corr:

  Logical; include covariance summaries.

- ...:

  Unused.

## Value

A `summary_JoiNMeMixFit` object inheriting from `summary_JoiNMeFit`.
