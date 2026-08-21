# Survival submodel

This directory tells the event-time part of the joint model.

Each cause-specific hazard combines a centred spline representation of the
baseline hazard, event covariates and the active longitudinal association
predictor. Exact, right, left and interval censoring share an interval-row
representation. Gauss--Kronrod quadrature evaluates the full event-time
likelihood rather than a Cox partial likelihood.

R scales only the integration coordinates and interval widths. It evaluates
the baseline-hazard basis and the ordinary Cox model matrix on original study
time before supplying them to Stan. `W` is the event-endpoint Cox design used
by an exact failure. `W_gk` is the corresponding design at every quadrature
ordinate, so a baseline covariate may interact with a time spline without
freezing that spline at the observed endpoint. Dynamic prediction receives
the same original-time node designs.

For a proper interval-censored observation \((L,R]\), the R preparation gives
Stan two adjacent risk rows. The first contributes survival from zero to
\(L\), and the second contributes the conditional probability of failure
between \(L\) and \(R\). Their product is
\(S(L)\{1-\exp[-(H(R)-H(L))]\}=S(L)-S(R)\). This distinction prevents the
lower inspection time from being mistaken for delayed entry.

## Read the blocks in order

1. [data](data/README.md) — censoring intervals, event covariates,
   baseline-hazard bases and quadrature designs.
2. [transformed data](transformed_data/README.md) — why no further
   parameter-free event preparation is needed.
3. [parameters](parameters/README.md) — cause-specific baseline-hazard and
   event-covariate coefficients.
4. [transformed parameters](transformed_parameters/README.md) — how retained
   event quantities enter the joint hazard calculation.
5. [model](model/README.md) — baseline-hazard priors and curvature penalties.
6. [generated quantities](generated_quantities/README.md) — event log
   likelihood, cumulative hazard and conditional survival.

The likelihood and future survival calculation remain at full-model level
because the event predictor consumes longitudinal effects, marker weights,
association coefficients and functional transformations.

[Back to the master submodel map](../README.md).
