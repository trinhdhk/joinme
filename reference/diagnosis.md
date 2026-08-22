# Diagnostic summary for JoiNMe objects

Returns diagnostic summaries for fitted (`JoiNMeFit`) and dynamic
prediction (`JoiNMeDynPred`) objects.

For `JoiNMeFit`, diagnostics include both an overall summary table and a
parameter-level diagnostics table built from the same cached posterior
summaries used by `summary.JoiNMeFit()`. This avoids the slower
backend-wide diagnostic pass and keeps the reported metrics aligned with
the summary sections users already inspect.

For `JoiNMeDynPred`, diagnostics summarise posterior-draw quality for
predicted quantities and are aligned to the same metric schema used for
`JoiNMeFit`.

## Usage

``` r
diagnosis(object, ...)
```

## Arguments

- object:

  A JoiNMe object.

- ...:

  Additional arguments passed to class-specific methods.

## Value

For `JoiNMeFit`, a `JoiNMe_diagnosis` object with components `summary`
and `by_parameter`. For `JoiNMeDynPred`, a data frame with `metric` and
`value` columns.
