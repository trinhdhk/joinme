# Concordance of conditional survival predictions

Estimates a single, follow-up-wide concordance index from the
conditional survival curves of a fitted joint model. The estimator
follows the survival-curve ordering of Antolini and colleagues: when
subject \\i\\ experiences the event before subject \\j\\, their pair is
concordant when \\\widehat S_i(r_i) \< \widehat S_j(r_i)\\, where
\\r_i\\ is the observed residual event time and both curves are
evaluated at that same residual time.

## Usage

``` r
# S3 method for class 'JoiNMeFit'
concordance(
  object,
  newdataLong = NULL,
  newdataEvent = NULL,
  time_start = NULL,
  cause = 1,
  predict_control = list(n_samples = 200, n_pred_draws = 20),
  seed = .Random.seed[[1]],
  type_weights = "none",
  ...
)
```

## Arguments

- object:

  A fitted object of class `JoiNMeFit`.

- newdataLong:

  Longitudinal evaluation data. The fitted longitudinal data are used
  when this and `newdataEvent` are both omitted.

- newdataEvent:

  Event-process evaluation data. The fitted event data are used when
  this and `newdataLong` are both omitted.

- time_start:

  Optional conditioning time. A scalar applies to every subject; a named
  numeric vector is matched by subject identifier; and a character
  scalar names a subject-constant column in `newdataEvent`. When
  omitted, each subject's latest longitudinal measurement strictly
  before their observed event or censoring time is used.

- cause:

  Positive integer identifying the event cause of interest. Other
  observed causes are treated as censoring at their event times.

- predict_control:

  Named list for prediction configuration.

- seed:

  Integer seed used when posterior draws are subsampled.

- type_weights:

  Character event-time weighting rule. `"none"` and `"n"` give Harrell's
  pair weighting. `"S"`, `"S/G"`, `"n/G2"` (Uno's weighting), and `"I"`
  have the meanings used by
  [`survival::concordance()`](https://rdrr.io/pkg/survival/man/concordance.html).

- ...:

  Additional arguments passed to
  [`predict.JoiNMeFit()`](https://trinhdhk.github.io/joinme/reference/predict.JoiNMeFit.md).

## Value

A one-row data frame with the concordance estimate, weighted counts of
concordant, discordant, and predictor-tied pairs, total comparison
weight `n_pairs`, and the numbers of analysed subjects and
cause-specific events.

## Details

Residual follow-up is measured from the conditioning time:
\\R_i=T_i-T\_{0i}\\. Conditional survival is one at residual time zero.
For every observed event at \\R_i\\, the method compares subject \\i\\
with subjects known to have survived beyond \\R_i\\. A subject censored
before \\R_i\\ is not comparable; a subject censored exactly at \\R_i\\
is known to be event-free at that time and remains comparable. Thus
premature censoring is never interpreted as exceptionally long survival.

The censoring and event-time weights are obtained from
[`survival::concordance()`](https://rdrr.io/pkg/survival/man/concordance.html).
Under the usual independent-censoring assumption, `"n/G2"` reduces
sensitivity to the censoring distribution by applying Uno's
inverse-censoring weighting. No method can recover the ordering of a
pair whose event-time order was not observed without further
assumptions.

The comparison deliberately evaluates both predicted curves at the
earlier event time. Evaluating each subject's curve at their own
observed outcome time and then ranking those values would use the
response to construct the predictor and can produce optimistically
biased concordance.

Dynamic predictions contain the complete fitted event model, including
baseline hazard, event covariates, latent longitudinal trajectories, and
every identity, functional, monotone-spline, or ordered piecewise-linear
association. Consequently no separate piecewise-linear approximation is
made by this method.

## References

Antolini L, Boracchi P, Biganzoli E (2005). A time-dependent
discrimination index for survival data. *Statistics in Medicine*, 24,
3927–3944. [doi:10.1002/sim.2427](https://doi.org/10.1002/sim.2427) .

## See also

[`tvROC.JoiNMeFit()`](https://trinhdhk.github.io/joinme/reference/tvROC.md),
[`tvAUC.JoiNMeFit()`](https://trinhdhk.github.io/joinme/reference/tvROC.md),
[`predict.JoiNMeFit()`](https://trinhdhk.github.io/joinme/reference/predict.JoiNMeFit.md),
[`survival::concordance()`](https://rdrr.io/pkg/survival/man/concordance.html)
