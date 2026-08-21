# Introduction to Joint Nested Mixed Effects Models

## Installation

``` r

# install.packages("remotes")
# remotes::install_github("trinhdhk/JoiNMe")
```

For Stan backends:

``` r

# CmdStanR backend (recommended)
# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# cmdstanr::install_cmdstan()

# RStan backend (optional)
# install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
```

## Quick start

## Conceptual overview

The `joinme` model links two components:

1.  A multivariate longitudinal model for repeated biomarker
    measurements.
2.  A survival model whose hazard depends on summaries of the
    longitudinal process at time t.

For the longitudinal part, each marker can follow its own distribution
(Gaussian, Student-t, Poisson, negative binomial, Bernoulli, beta, skew
families, or ordinal). The linear predictor is decomposed into fixed
effects, subject-level random effects, marker-level random effects, and
optional marker-by-subject random effects. The survival part uses a
proportional hazards model with a flexible spline baseline hazard and
time-dependent association terms.

The association terms implemented in `joinme` package are:

- Current value (CV) of the latent linear predictor.
- Current slope (CS) of the latent linear predictor.
- Correlation association (CORR) derived from off-diagonal
  marker-by-subject correlations.
- Variance-covariance association (VCOV) derived from lower-triangular
  marker-by-subject Cholesky-factor entries taken directly from `L`.

Each association term can be left on its original scale or transformed
using a functional expression, a monotone I-spline, or an ordered
piecewise-linear mapping.

## Family-shared distributional semantics

For mixed longitudinal families, `joinme` treats distributional
parameters as **shared by family** when no distributional regression is
supplied. In brief,

- Markers with the same family share one latent baseline process for
  each applicable distributional parameter.
- Markers from different families do not share those baselines.
- Parameters are instantiated only when a family uses them (for example
  `nu` for Student-`t`, `phi` for negative-binomial).

Example interpretation:

- If markers `A` and `B` are both Student-`t` and marker `C` is
  Gaussian, then `A` and `B` share the same baseline `sigma` and `nu`
  processes, while `C` has its own `sigma` process and no `nu`.

### What changes when `formulaDist` is provided?

Supplying a distributional regression (for example `sigma ~ 1 + time`)
regresses the parameter from a family-level constant baseline to a
regression-driven process. Summaries report distributional regression
coefficients (`(Intercept)`, covariate effects, and random-effect SD
terms when present) rather than only baseline family constants.

Illustrative `formulaDist` specification:

``` r

formulaDist <- list(
  sigma ~ 1 + time,
  nu ~ 1
)
```

This means:

- `sigma` is modelled through regression (fixed/random effects if
  included),
- `nu` is regression-driven with an intercept-only structure,
- other distributional parameters without formulas continue using
  family-shared baselines when applicable.

## Data

You need two data frames:

- `dataLong`: long-format longitudinal data with columns for subject id,
  marker id, time, outcome, and covariates.
- `dataEvent`: one row per subject for an ordinary right-, left- or
  interval-censored outcome, or ordered risk intervals for a
  counting-process outcome, together with the event covariates.

The event formula accepts `Surv(time, status)`,
`Surv(start, stop, status)`, `Surv(time, status, type = "left")`, and
`Surv(lower, upper, type = "interval2")`. An interval-censored event in
((L,R\]) contributes (S(L)-S(R)) to the full event likelihood. Its lower
inspection limit is not interpreted as delayed entry. See [Interpreting
a fitted JoiNMe
model](https://trinhdhk.github.io/joinme/articles/joinme-interpretation.md)
for the censoring equations and worked interpretation of conditional
effects and contrasts.

Formulae and reported coefficients use the original study-time scale.

### Simulation example

We generate a modest dataset so the workflow runs quickly while still
containing both longitudinal and survival information.

``` r

sim <- simulate_joinme(
  n_id = 10,
  families = list(
    jm_family("student_t"),
    jm_family("student_t")
  ),
  times_obs = seq(0, 5, length.out = 8),
  quadrature_nodes = 31,
  seed = 1,
  assoc = c("cv_total"),
  truth = jm_truth(
    assoc_coef = list(slope = c(cv_total = 0.6)),
    basehaz = list(type = "weibull", shape = 2, scale = 8)
  )
)

head(sim$dataLong)
```

      id marker      time         x1       x2          y
    1  1     m1 0.0000000 -0.6264538 1.511781  4.6351130
    2  1     m1 0.7142857 -0.6264538 1.511781  1.0454417
    3  1     m1 1.4285714 -0.6264538 1.511781 -0.6578417
    4  1     m1 2.1428571 -0.6264538 1.511781  2.6007495
    5  1     m1 2.8571429 -0.6264538 1.511781  3.3550153
    6  1     m1 3.5714286 -0.6264538 1.511781  1.5404495

``` r

head(sim$dataEvent)
```

      id         x1          x2 time event time_start time_stop
    1  1 -0.6264538  1.51178117    8     0          0         8
    2  2  0.1836433  0.38984324    8     0          0         8
    3  3 -0.8356286 -0.62124058    8     0          0         8
    4  4  1.5952808 -2.21469989    8     0          0         8
    5  5  0.3295078  1.12493092    8     0          0         8
    6  6 -0.8204684 -0.04493361    8     0          0         8

The longitudinal formula defines fixed effects and random effects. The
survival formula defines baseline covariates in the hazard model.

Formula-scoped weighting can be declared directly in grouping terms via
`weighted(group, weights = <column>)`, including marker-level grouping
declarations.

``` r

formulaLong <- y ~ 1 + time + x1 +
  (1 + time || id) +
  (0 + x1 + (1 + time || id) || marker)

formulaEvent <- survival::Surv(time, event) ~ 1 + x1 + x2
```

``` r

formulaLong_weighted <- y ~ 1 + time + x1 +
  (1 + time | weighted(id, weights = id_w)) +
  (0 + x1 + (1 + time | weighted(id, weights = id_w)) | weighted(marker, weights = marker_w))
```

We link the survival hazard to the current value (CV) of the
longitudinal process using an identity transformation.

``` r

fit <- joinme(
  formulaLong = formulaLong,
  dataLong = sim$dataLong,
  formulaEvent = formulaEvent,
  dataEvent = sim$dataEvent,
  assoc = c("cv_total"),
  families = "student_t",
  transforms = joinme_tf(cv_total = "identity"),
  basehaz = joinme_basehaz('bs'), 
  control = list(
    engine = "rstan",
    chains = 1,
    parallel_chains = 1,
    iter_warmup = 500,
    iter_sampling = 500,
    seed = 1,
    refresh = 0,
    adapt_delta = 0.85
  )
)
```

To fit only the multivariate nested longitudinal process, omit both
event arguments. The formula continues to define population, subject,
marker, and subject-by-marker effects, but no survival likelihood or
association coefficient is introduced.

``` r

longitudinal_fit <- joinme(
  formulaLong = formulaLong,
  dataLong = sim$dataLong,
  families = "student_t",
  control = list(engine = "cmdstanr")
)
```

Longitudinal summaries, posterior predictions and trajectory plots
retain their ordinary interpretation. Event prediction, survival plots,
concordance, and time-dependent ROC/AUC are deliberately unavailable
because no event outcome contributed to the likelihood.

``` r

summary(fit) # or posterior_summary(fit)
```


    Joint mixed effects model summary
    • Call: joinme(formulaLong = formulaLong, dataLong = sim$dataLong, formulaEvent = formulaEvent, dataEvent = sim$dataEvent, control = list(engine = "rstan", chains = 1, parallel_chains = 1, iter_warmup = 500, iter_sampling = 500, seed = 1, refresh = 0, adapt_delta = 0.85), families = "student_t", transforms = joinme_tf(cv_total = "identity"), basehaz = joinme_basehaz("bs"), assoc = c("cv_total"))
    • Family: student_t, student_t
    • Event process: joint longitudinal-survival
    • Baseline hazard type: bs
    • tmax: 8
    • Marker-weight offsets [shared]: m1=0, m2=0

    Sampler diagnostics
    # A tibble: 11 × 2
       metric               value
       <chr>                <chr>
     1 draws                500
     2 divergences          11
     3 treedepth_hits       0
     4 ebfmi_min            0.865
     5 max_rhat             1.027
     6 min_ess_bulk         131.299
     7 min_ess_tail         118.267
     8 n_terms_total        28
     9 n_terms_bad_rhat     2
    10 n_terms_low_ess_bulk 0
    11 n_terms_low_ess_tail 0

    Transformations
    # A tibble: 1 × 2
      term     formula
      <chr>    <chr>
    1 cv_total identity

    Fixed effects (beta)
    # A tibble: 3 × 8
      term        Estimate Est.Error   Q2.5 Q97.5  Rhat ess_bulk ess_tail
      <chr>          <dbl>     <dbl>  <dbl> <dbl> <dbl>    <dbl>    <dbl>
    1 (Intercept)    1.81      0.668  0.408 2.84  0.998     356.     188.
    2 time           0.188     0.286 -0.369 0.739 1.00      286.     272.
    3 x1             0.81      0.529 -0.367 1.73  1.00      402.     320.

    Baseline hazard coefficients
    # A tibble: 6 × 11
      term    Estimate Hazard.Ratio Est.Error    Q2.5 Q97.5 HR.Q2.5 HR.Q97.5  Rhat
      <chr>      <dbl>        <dbl>     <dbl>   <dbl> <dbl>   <dbl>    <dbl> <dbl>
    1 basis_1   -5.79         0.003     3.48  -13.2   0.776   0         2.17 1.00
    2 basis_2   -0.588        0.555     0.915  -2.63  1.02    0.072     2.78 1.00
    3 basis_3   -0.303        0.739     0.631  -1.66  0.797   0.191     2.22 1.00
    4 basis_4   -0.039        0.962     0.537  -1.01  1.04    0.364     2.83 1.01
    5 basis_5    0.195        1.22      0.657  -0.983 1.43    0.374     4.19 0.999
    6 basis_6    0.416        1.52      0.97   -1.36  2.51    0.256    12.3  1.01
      ess_bulk ess_tail
         <dbl>    <dbl>
    1     530.     339.
    2     267.     154.
    3     256.     118.
    4     296.     143.
    5     385.     306.
    6     395.     247.

    Survival process (non-association covariates)
    # A tibble: 2 × 11
      term  Estimate Hazard.Ratio Est.Error  Q2.5 Q97.5 HR.Q2.5 HR.Q97.5  Rhat
      <chr>    <dbl>        <dbl>     <dbl> <dbl> <dbl>   <dbl>    <dbl> <dbl>
    1 x1      -0.052        0.949     0.829 -1.59  1.49   0.203     4.43  1.02
    2 x2      -0.205        0.815     0.94  -2.09  1.53   0.123     4.62  1.00
      ess_bulk ess_tail
         <dbl>    <dbl>
    1     855.     315.
    2     436.     348.

    Association parameters
    # A tibble: 1 × 8
      term     Estimate Est.Error  Q2.5 Q97.5  Rhat ess_bulk ess_tail
      <chr>       <dbl>     <dbl> <dbl> <dbl> <dbl>    <dbl>    <dbl>
    1 cv_total     1.02     0.574 0.191  2.51  1.01     227.     193.

    Marker-weight distribution
    # A tibble: 2 × 8
      term        Estimate Est.Error   Q2.5 Q97.5  Rhat ess_bulk ess_tail
      <chr>          <dbl>     <dbl>  <dbl> <dbl> <dbl>    <dbl>    <dbl>
    1 mean weight    -1.60     1.40  -4.90  0.365 1         268.     129.
    2 SD weight       1.03     0.536  0.263 2.26  0.998     258.     248.

    Distributional parameters
    # A tibble: 2 × 8
      term                  Estimate Est.Error  Q2.5 Q97.5  Rhat ess_bulk ess_tail
      <chr>                    <dbl>     <dbl> <dbl> <dbl> <dbl>    <dbl>    <dbl>
    1 sigma[family=student]     1.48     0.153  1.20  1.78  1.00     333.     472.
    2 nu[family=student]        2.61     0.503  2.02  3.77  1.01     256.     172.

    Covariance summaries (diagonal entries are variances)

    id
    # A tibble: 4 × 10
      block row         col         Estimate Est.Error     Q2.5 Q97.5   Rhat
      <chr> <chr>       <chr>          <dbl>     <dbl>    <dbl> <dbl>  <dbl>
    1 id    (Intercept) (Intercept)    0.190     0.312 0.000233  1.10  0.999
    2 id    (Intercept) time           0         0     0         0    NA
    3 id    time        (Intercept)    0         0     0         0    NA
    4 id    time        time           0.624     0.505 0.0479    1.97  1.00
      ess_bulk ess_tail
         <dbl>    <dbl>
    1     349.     333.
    2      NA       NA
    3      NA       NA
    4     243.     162.

    marker
    # A tibble: 4 × 10
      block  row         col         Estimate Est.Error     Q2.5 Q97.5   Rhat
      <chr>  <chr>       <chr>          <dbl>     <dbl>    <dbl> <dbl>  <dbl>
    1 marker x1          x1             0.608      1.30 0.000119  4.55  1.00
    2 marker x1          (Intercept)    0          0    0         0    NA
    3 marker (Intercept) x1             0          0    0         0    NA
    4 marker (Intercept) (Intercept)    0.731      2.14 0.000169  6.32  1.000
      ess_bulk ess_tail
         <dbl>    <dbl>
    1     459.     203.
    2      NA       NA
    3      NA       NA
    4     310.     214.

    id:marker covariance parameters

    covariance regression coefficients
    # A tibble: 4 × 11
      block         row         col         term               Estimate Est.Error
      <chr>         <chr>       <chr>       <chr>                 <dbl>     <dbl>
    1 SD[id:marker] (Intercept) (Intercept) (Intercept)          -1.87      1.59
    2 SD[id:marker] time        time        (Intercept)          -0.039     0.494
    3 SD[id:marker] (Intercept) (Intercept) latent SD (lambda)    0.74      0.636
    4 SD[id:marker] time        time        latent SD (lambda)    0.514     0.403
        Q2.5 Q97.5  Rhat ess_bulk ess_tail
       <dbl> <dbl> <dbl>    <dbl>    <dbl>
    1 -5.42  0.458 1.00      282.     310.
    2 -0.998 0.962 0.999     159.     350.
    3  0.036 2.52  1.00      443.     316.
    4  0.017 1.49  1.00      323.     292.

``` r

fixef(fit)
```

          component event assoc_term        term Estimate Est.Error   Q2.5 Q97.5
    1  longitudinal  <NA>       <NA> (Intercept)    1.811     0.668  0.408 2.837
    2  longitudinal  <NA>       <NA>        time    0.188     0.286 -0.369 0.739
    3  longitudinal  <NA>       <NA>          x1    0.810     0.529 -0.367 1.728
    4         event event       <NA>          x1   -0.052     0.829 -1.593 1.489
    5         event event       <NA>          x2   -0.205     0.940 -2.094 1.531
    6 marker_weight  <NA>     shared mean weight   -1.605     1.396 -4.899 0.365
       Rhat ess_bulk ess_tail
    1 0.998 355.9659 188.4075
    2 1.002 286.0901 271.6760
    3 1.005 401.7544 320.4063
    4 1.022 855.3959 314.6455
    5 1.001 436.1793 347.6246
    6 1.000 268.2684 128.6597

``` r

coef(fit)
```

    $formulaLong
    $formulaLong$id
       id        term Estimate Est.Error   Q2.5 Q97.5  Rhat ess_bulk ess_tail
    1   1 (Intercept)    1.857     0.700  0.451 2.954 0.999 230.9204 193.2691
    2   1        time    0.404     0.494 -0.476 1.361 1.005 292.2982 372.0344
    3  10 (Intercept)    1.769     0.726  0.275 2.917 0.999 215.7018 169.7796
    4  10        time   -0.542     0.587 -1.724 0.672 1.006 240.9057 280.6671
    5   2 (Intercept)    1.767     0.705 -0.084 2.825 1.000 238.7136 208.2397
    6   2        time    0.085     0.443 -0.842 0.895 1.008 276.1880 247.7924
    7   3 (Intercept)    1.911     0.739  0.434 3.257 0.999 198.6921 164.1739
    8   3        time    0.591     0.434 -0.298 1.379 1.015 314.9444 356.0969
    9   4 (Intercept)    1.831     0.756  0.179 3.153 1.004 222.8757 180.6438
    10  4        time    0.808     0.486 -0.130 1.777 0.999 257.5783 329.3857
    11  5 (Intercept)    1.787     0.746  0.155 2.896 0.998 212.2547 216.9369
    12  5        time    0.918     0.622 -0.316 2.187 0.998 273.2440 350.3450
    13  6 (Intercept)    1.664     0.723  0.155 2.772 1.000 178.8404 178.2656
    14  6        time    0.172     0.410 -0.639 0.959 1.004 338.2438 356.0969
    15  7 (Intercept)    1.930     0.717  0.458 3.220 1.000 195.4408 193.4196
    16  7        time   -0.459     0.564 -1.504 0.653 1.005 210.6477 370.1080
    17  8 (Intercept)    1.972     0.717  0.384 3.210 1.000 215.0235 186.6017
    18  8        time   -0.413     0.454 -1.230 0.569 1.011 244.7537 326.7148
    19  9 (Intercept)    1.726     0.733  0.314 2.870 1.000 191.3919 191.9124
    20  9        time    0.564     0.488 -0.377 1.610 1.000 306.1700 457.7611

    $formulaLong$marker
      marker        term Estimate Est.Error   Q2.5 Q97.5  Rhat ess_bulk ess_tail
    1     m1 (Intercept)    2.001     0.415  1.156 2.788 0.998 619.3644 361.2798
    2     m1          x1    0.961     0.522 -0.171 2.023 0.999 606.0756 265.3052
    3     m2 (Intercept)    1.909     0.394  1.128 2.641 1.000 659.7722 390.3540
    4     m2          x1    0.944     0.454  0.136 1.810 1.002 544.6139 445.2679

    $formulaLong$marker_by_id
       id marker        term Estimate Est.Error   Q2.5 Q97.5  Rhat ess_bulk
    1   1     m1 (Intercept)    1.836     0.738  0.247 3.041 1.002 195.3303
    2   1     m1        time   -0.147     0.514 -1.125 0.805 1.004 353.9033
    3   1     m2 (Intercept)    1.832     0.761  0.416 3.163 1.000 199.2666
    4   1     m2        time    0.832     0.527 -0.205 1.881 1.019 387.2622
    5  10     m1 (Intercept)    1.858     0.749  0.253 3.074 0.999 208.4218
    6  10     m1        time    0.087     0.634 -1.201 1.315 1.001 271.2302
    7  10     m2 (Intercept)    1.723     0.737  0.203 2.893 0.999 235.9054
    8  10     m2        time   -0.565     0.672 -2.004 0.662 1.003 247.3512
    9   2     m1 (Intercept)    1.696     0.723  0.109 2.835 0.998 223.4939
    10  2     m1        time   -0.208     0.514 -1.257 0.759 1.004 305.0752
    11  2     m2 (Intercept)    1.916     0.772  0.276 3.338 1.001 193.2944
    12  2     m2        time    0.483     0.527 -0.580 1.552 1.003 366.9827
    13  3     m1 (Intercept)    1.853     0.800  0.128 3.416 0.999 172.5930
    14  3     m1        time    0.620     0.480 -0.345 1.546 1.000 420.5334
    15  3     m2 (Intercept)    1.826     0.737  0.216 2.992 0.999 186.0429
    16  3     m2        time    0.153     0.463 -0.767 1.064 0.999 436.6906
    17  4     m1 (Intercept)    1.927     0.855  0.194 3.570 1.004 179.7751
    18  4     m1        time    0.288     0.528 -0.746 1.261 1.003 282.3756
    19  4     m2 (Intercept)    1.740     0.785  0.053 3.020 0.998 229.9107
    20  4     m2        time    0.642     0.537 -0.382 1.632 1.007 304.1326
    21  5     m1 (Intercept)    1.869     0.762  0.380 3.204 1.002 220.2330
    22  5     m1        time    1.627     0.665  0.422 2.909 0.999 258.0237
    23  5     m2 (Intercept)    1.705     0.720  0.004 2.858 0.998 169.2496
    24  5     m2        time    0.079     0.641 -1.110 1.256 0.999 317.9872
    25  6     m1 (Intercept)    1.765     0.768 -0.044 3.151 1.001 175.0023
    26  6     m1        time   -0.100     0.485 -1.052 0.863 1.010 399.8949
    27  6     m2 (Intercept)    1.717     0.703  0.056 2.717 1.001 168.4039
    28  6     m2        time    0.468     0.468 -0.425 1.402 0.999 408.3611
    29  7     m1 (Intercept)    1.851     0.782  0.138 3.276 0.999 208.5633
    30  7     m1        time   -0.870     0.639 -2.165 0.384 1.000 249.6671
    31  7     m2 (Intercept)    1.919     0.785  0.494 3.331 0.998 189.3267
    32  7     m2        time    0.229     0.605 -0.941 1.485 1.002 237.9042
    33  8     m1 (Intercept)    1.893     0.740  0.371 3.231 0.999 198.1855
    34  8     m1        time    0.073     0.466 -0.865 1.046 1.002 252.4152
    35  8     m2 (Intercept)    1.927     0.756  0.448 3.300 1.001 187.8563
    36  8     m2        time   -0.266     0.483 -1.226 0.607 1.000 265.3006
    37  9     m1 (Intercept)    1.695     0.760 -0.124 2.804 0.999 191.4044
    38  9     m1        time    0.017     0.508 -1.068 0.920 0.999 322.6817
    39  9     m2 (Intercept)    1.836     0.778  0.294 3.228 0.999 203.2314
    40  9     m2        time    0.727     0.586 -0.371 1.843 1.001 286.5473
       ess_tail
    1  150.1085
    2  411.9569
    3  177.2075
    4  358.8392
    5  180.9275
    6  229.4178
    7  254.4247
    8  376.8660
    9  225.1493
    10 259.3056
    11 170.5551
    12 311.3205
    13 175.5882
    14 343.4290
    15 213.3680
    16 452.7353
    17 194.5557
    18 298.8751
    19 228.3575
    20 307.6760
    21 158.1732
    22 412.2037
    23 174.5643
    24 461.2318
    25 174.5643
    26 370.7469
    27 204.5691
    28 360.4133
    29 180.9275
    30 295.5926
    31 179.4603
    32 279.1936
    33 150.4984
    34 327.3183
    35 180.9275
    36 261.4892
    37 219.6419
    38 212.3503
    39 170.3290
    40 408.3394


    $assoc
      assoc_term term Estimate Est.Error   Q2.5  Q97.5  Rhat ess_bulk ess_tail
    1     shared   m1   -1.179     1.281 -4.604  0.191 0.999 308.4560 158.4132
    2     shared   m2   -2.507     1.337 -5.835 -0.653 1.004 364.2494 161.3811
              group
    1 marker_weight
    2 marker_weight

    $formulaEvent
      event term Estimate Est.Error   Q2.5 Q97.5  Rhat ess_bulk ess_tail
    1 event   x1   -0.052     0.829 -1.593 1.489 1.022 814.5392 314.6455
    2 event   x2   -0.205     0.940 -2.094 1.531 1.001 426.9095 347.6246
      Hazard.Ratio HR.Q2.5 HR.Q97.5
    1        0.949   0.203    4.433
    2        0.815   0.123    4.623

    $formulaDist
    list()

    $formulaVCov
    $formulaVCov$population
    NULL

    $formulaVCov$id
       id         block row col        term Estimate Est.Error   Q2.5 Q97.5  Rhat
    1   1 SD[id:marker]   1   1 (Intercept)   -0.123     0.966 -2.330 1.551 1.000
    2   1 SD[id:marker]   2   2 (Intercept)    0.038     0.624 -1.023 1.737 1.018
    3  10 SD[id:marker]   1   1 (Intercept)    0.009     0.921 -1.888 2.004 1.000
    4  10 SD[id:marker]   2   2 (Intercept)    0.025     0.617 -1.484 1.359 1.003
    5   2 SD[id:marker]   1   1 (Intercept)   -0.033     0.991 -2.463 2.105 1.005
    6   2 SD[id:marker]   2   2 (Intercept)   -0.039     0.552 -1.253 1.157 1.002
    7   3 SD[id:marker]   1   1 (Intercept)   -0.091     1.029 -2.835 2.084 1.003
    8   3 SD[id:marker]   2   2 (Intercept)   -0.117     0.578 -1.729 0.966 0.999
    9   4 SD[id:marker]   1   1 (Intercept)   -0.007     0.879 -1.925 1.991 1.003
    10  4 SD[id:marker]   2   2 (Intercept)   -0.132     0.626 -1.585 1.103 1.003
    11  5 SD[id:marker]   1   1 (Intercept)    0.072     0.904 -1.971 2.018 1.010
    12  5 SD[id:marker]   2   2 (Intercept)    0.368     0.589 -0.456 1.764 1.002
    13  6 SD[id:marker]   1   1 (Intercept)   -0.051     0.865 -1.972 1.648 0.999
    14  6 SD[id:marker]   2   2 (Intercept)   -0.150     0.641 -1.624 0.983 1.000
    15  7 SD[id:marker]   1   1 (Intercept)   -0.068     0.943 -2.163 1.624 1.003
    16  7 SD[id:marker]   2   2 (Intercept)    0.229     0.598 -0.855 1.755 0.998
    17  8 SD[id:marker]   1   1 (Intercept)   -0.037     0.868 -2.153 1.452 1.013
    18  8 SD[id:marker]   2   2 (Intercept)   -0.187     0.687 -1.968 0.956 1.001
    19  9 SD[id:marker]   1   1 (Intercept)   -0.073     0.987 -2.623 1.949 1.004
    20  9 SD[id:marker]   2   2 (Intercept)   -0.001     0.578 -1.301 1.231 1.005
       ess_bulk ess_tail
    1  445.9349 424.3182
    2  395.1673 357.6137
    3  357.3181 234.7502
    4  444.8508 330.4114
    5  573.4763 465.4078
    6  279.0113 292.2525
    7  619.5966 427.9046
    8  441.7741 451.1299
    9  425.2130 255.6396
    10 436.2357 334.5920
    11 469.7286 312.9874
    12 240.9191 245.8280
    13 416.2875 311.3205
    14 386.6729 266.8336
    15 519.7571 398.0777
    16 292.6329 188.5822
    17 585.0811 473.4051
    18 377.5392 450.4360
    19 286.4082 199.5485
    20 414.0917 431.8954

``` r

posterior_assoc(fit, summary = TRUE)
```


    Posterior association effects

    Association term: cv_total
    # A tibble: 2 × 8
      term  Estimate Est.Error  Q2.5  Q97.5  Rhat ess_bulk ess_tail
      <chr>    <dbl>     <dbl> <dbl>  <dbl> <dbl>    <dbl>    <dbl>
    1 m1      -0.434     0.388 -1.47  0.125 1.00      395.     193.
    2 m2      -1.11      0.708 -2.66 -0.302 0.999     366.     332.

`posterior_summary(what = "model")` returns the same model summary as
[`summary()`](https://rdrr.io/r/base/summary.html). Its `fixef`,
`ranef`, and `coef` views form the common summary layer used by those
three conventional methods.
[`posterior_assoc()`](https://trinhdhk.github.io/joinme/reference/assoc.md)
returns association-specific posterior output in a format that stays
close to the main summary tables.

For covariance-style channels (`corr`, `vcov`), the printed and tabular
output uses readable matrix labels instead of raw numeric indices.
[`coef()`](https://rdrr.io/r/stats/coef.html) returns coefficients on
the model-matrix scale after adding the population-level contribution
and the matching group-level deviation. Use `summary = FALSE` to
retrieve posterior draws instead of posterior summaries.

## Dynamic prediction

We condition on the observed history up to the last measurement for
subject 1 and obtain predictive summaries beyond that time.

``` r

ndL <- sim$dataLong[sim$dataLong$id == 1, ]
ndE <- sim$dataEvent[sim$dataEvent$id == 1, ]
time_start <- max(ndL$time)

pred <- posterior_predict(
  fit,
  newdataLong = ndL,
  newdataEvent = ndE,
  time_start = time_start,
  times = seq(time_start, time_start + 1, length.out = 20),
  control = list(
    n_samples = 50, 
    n_pred_draws = 1000,
    iter_warmup = 200,
    chains = 1
    ),
  seed = 2026
)
```


      |
      |                                                            |   0%


      |
      |============================================================| 100%

``` r

# control$n_samples: extracted posterior parameter draws
# control$n_pred_draws: dynpred output draw count

# use predict(..., scale = c("epred", "linpred", "predict"))
# to return multiple longitudinal scales in one call (default = all scales)
```

``` r

plot(pred, type = c("longitudinal", "survival"), combined = TRUE)
```

![](joinme-intro_files/figure-html/intro-plot-1.png)

The heatmap is an alternative to the usual separated marker curves. It
orders markers by the magnitude of posterior mean change and makes tiles
transparent when the posterior sign certainty does not exceed
`1 - threshold`.

`plot(fit, ...)` uses fitted posterior samples directly and has no
conditioning argument. Conditional covariate profiles are evaluated with
[`conditional_effects()`](https://paulbuerkner.com/brms/reference/conditional_effects.brmsfit.html),
while [`predict()`](https://rdrr.io/r/stats/predict.html) remains
reserved for dynamic subject-specific forecasts given observed
longitudinal history.

``` r

profiles <- make_conditions(sim$dataEvent, vars = "x1")

effects_data <- conditional_effects(
  fit,
  effects = list(longitudinal = "time", event = "x1"),
  conditions = profiles,
  process = c("longitudinal", "event"),
  longitudinal_estimand = "population",
  plot = FALSE
)

plot(effects_data, ask = FALSE)
```

For the longitudinal process, set `longitudinal_estimand = "marker"` to
add the fitted marker-level deviation, or use `"marginal_marker"` to
average marker-specific predictions within each posterior draw. This
latter estimand applies each marker’s inverse link before averaging, so
it represents the average expected marker response rather than the
response of an average marker.

[`conditional_contrast()`](https://trinhdhk.github.io/joinme/reference/conditional_contrast.md)
compares two complete or partial covariate profiles. The two predictions
are evaluated against the same posterior draw and are contrasted before
the posterior centre and interval are calculated. Hence the reported
uncertainty retains posterior dependence between the profiles. The
longitudinal result is group A minus group B on the chosen response or
linear-predictor scale. The event result is the hazard ratio of A
relative to B, or its logarithm when `event_scale = "log_hazard_ratio"`.

``` r

contrast_conditions <- data.frame(
  time = seq(0, max(sim$dataLong$time), length.out = 50),
  x2 = 0,
  cond__ = paste0("time=", round(seq(0, max(sim$dataLong$time), length.out = 50), 2))
)

contrast_data <- conditional_contrast(
  fit,
  groupA = c(x1 = 1),
  groupB = c(x1 = -1),
  conditions = contrast_conditions,
  process = c("longitudinal", "event"),
  method = "posterior_epred",
  longitudinal_estimand = "marginal_marker",
  event_scale = "hazard_ratio",
  plot = FALSE
)

plot(contrast_data, condition_variable = "time", ask = FALSE)
```

Each row of `conditions` supplies covariates common to the two groups.
Undeclared predictors are represented by their mean in the fitting data
for a numeric predictor and by their first fitted level for a
categorical predictor. For the event process, the contrast concerns the
direct covariate component of `formulaEvent`; the baseline hazard and
longitudinal association contribution are common to the two profiles.

### Step 8: Diagnostics and extraction helpers

``` r

diagnosis(fit)
```


    Diagnostics Summary
    # A tibble: 11 × 2
       metric               value
       <chr>                <chr>
     1 draws                500
     2 divergences          11
     3 treedepth_hits       0
     4 ebfmi_min            0.865
     5 max_rhat             1.027
     6 min_ess_bulk         131.299
     7 min_ess_tail         118.267
     8 n_terms_total        28
     9 n_terms_bad_rhat     2
    10 n_terms_low_ess_bulk 0
    11 n_terms_low_ess_tail 0

    Per-Parameter Diagnostics
    # A tibble: 20 × 13
       family           parameter             term                  Estimate
       <chr>            <chr>                 <chr>                    <dbl>
     1 fixef            (Intercept)           (Intercept)              1.81
     2 fixef            time                  time                     0.188
     3 fixef            x1                    x1                       0.81
     4 baseline_hazard  basis_1               basis_1                 -5.79
     5 baseline_hazard  basis_2               basis_2                 -0.588
     6 baseline_hazard  basis_3               basis_3                 -0.303
     7 baseline_hazard  basis_4               basis_4                 -0.039
     8 baseline_hazard  basis_5               basis_5                  0.195
     9 baseline_hazard  basis_6               basis_6                  0.416
    10 survival_process x1                    x1                      -0.052
    11 survival_process x2                    x2                      -0.205
    12 assoc            cv_total              cv_total                 1.02
    13 marker_weights   mean weight           mean weight             -1.60
    14 marker_weights   SD weight             SD weight                1.03
    15 distributional   sigma[family=student] sigma[family=student]    1.48
    16 distributional   nu[family=student]    nu[family=student]       2.61
    17 corr             id                    <NA>                     0.190
    18 corr             id                    <NA>                     0
    19 corr             id                    <NA>                     0
    20 corr             id                    <NA>                     0.624
       Est.Error       Q2.5 Q97.5   Rhat ess_bulk ess_tail block row
           <dbl>      <dbl> <dbl>  <dbl>    <dbl>    <dbl> <chr> <chr>
     1     0.668   0.408    2.84   0.998     356.     188. <NA>  <NA>
     2     0.286  -0.369    0.739  1.00      286.     272. <NA>  <NA>
     3     0.529  -0.367    1.73   1.00      402.     320. <NA>  <NA>
     4     3.48  -13.2      0.776  1.00      530.     339. <NA>  <NA>
     5     0.915  -2.63     1.02   1.00      267.     154. <NA>  <NA>
     6     0.631  -1.66     0.797  1.00      256.     118. <NA>  <NA>
     7     0.537  -1.01     1.04   1.01      296.     143. <NA>  <NA>
     8     0.657  -0.983    1.43   0.999     385.     306. <NA>  <NA>
     9     0.97   -1.36     2.51   1.01      395.     247. <NA>  <NA>
    10     0.829  -1.59     1.49   1.02      855.     315. <NA>  <NA>
    11     0.94   -2.09     1.53   1.00      436.     348. <NA>  <NA>
    12     0.574   0.191    2.51   1.01      227.     193. <NA>  <NA>
    13     1.40   -4.90     0.365  1         268.     129. <NA>  <NA>
    14     0.536   0.263    2.26   0.998     258.     248. <NA>  <NA>
    15     0.153   1.20     1.78   1.00      333.     472. <NA>  <NA>
    16     0.503   2.02     3.77   1.01      256.     172. <NA>  <NA>
    17     0.312   0.000233 1.10   0.999     349.     333. id    (Intercept)
    18     0       0        0     NA          NA       NA  id    (Intercept)
    19     0       0        0     NA          NA       NA  id    time
    20     0.505   0.0479   1.97   1.00      243.     162. id    time
       col
       <chr>
     1 <NA>
     2 <NA>
     3 <NA>
     4 <NA>
     5 <NA>
     6 <NA>
     7 <NA>
     8 <NA>
     9 <NA>
    10 <NA>
    11 <NA>
    12 <NA>
    13 <NA>
    14 <NA>
    15 <NA>
    16 <NA>
    17 (Intercept)
    18 time
    19 (Intercept)
    20 time
    ... truncated to 20 rows; inspect $by_parameter for the full table.

``` r

diagnosis(pred)
```

                     metric       value
    1                 draws 1000.000000
    2           divergences    0.000000
    3        treedepth_hits    0.000000
    4             ebfmi_min    1.450653
    5              max_rhat    1.051430
    6          min_ess_bulk   16.105405
    7          min_ess_tail   11.647129
    8         n_terms_total    7.000000
    9      n_terms_bad_rhat    4.000000
    10 n_terms_low_ess_bulk    2.000000
    11 n_terms_low_ess_tail    4.000000

`ranef(fit)`, `coef(fit)`, and `vcov(fit)` are nested by `formulaLong`
and `formulaDist`. Set `summary = FALSE` in
[`fixef()`](https://rdrr.io/pkg/nlme/man/fixed.effects.html),
[`ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html), or
[`coef()`](https://rdrr.io/r/stats/coef.html) for draw-level output. The
same views are available one layer lower through
`posterior_summary(what = ..., summary = FALSE)` and through the
corresponding structured selectors in
[`extract()`](https://trinhdhk.github.io/joinme/reference/extract.md).
For prediction objects, `ranef(pred)` and `vcov(pred)` are available
only when marker covariance depends on id.
