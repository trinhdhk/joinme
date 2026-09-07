# The JoiNMe Framework: Methodology and Statistical Theory

## Abstract

The `joinme` package provides a Bayesian framework for joint analysis of
multivariate longitudinal biomarkers and time-to-event outcomes. The
longitudinal process is expressed through a layered random-effects
decomposition that separates subject-level signals, marker-level
deviations, and marker-by-subject departures.

The survival component links risk to summaries of the longitudinal
trajectories, most notably the current value (CV), current slope (CS),
and a variance summary of the marker-by-subject covariance structure
through correlation-style (`corr`) and covariance-style (`vcov`)
channels. Association summaries can be transformed by user-defined
mathematical expressions, monotone splines, or ordered piecewise-linear
rules. Bayesian estimation propagates uncertainty jointly across the
longitudinal, event and association components.

In the public summaries and plotting output, covariance-style
associations are reported on the effective hazard scale used by the
fitted model. This keeps the reported `corr[...]` and `vcov[...]`
coefficients aligned with the printed matrix labels, association plots,
and posterior extraction helpers.

For monotone spline transforms, `JoiNMe` supports both raw-scale
I-splines and bounded-domain I-splines evaluated on `plogis(x)`. The
bounded-domain variant keeps the spline input inside (0, 1), which often
yields more stable knot placement and smoother optimisation when raw
association features are wide or heavy-tailed.

## Introduction

Joint modelling of longitudinal and survival data is a standard approach
for analysing time-dependent covariates measured with error \[1\]. In
many clinical and observational settings, subject status is monitored
via multiple biomarkers (for example, blood pressure, bilirubin,
cognitive scores) alongside a terminal event (for example, death,
dropout). Standard approaches often face two limitations:

1.  **Computational Scalability**: As the number of markers D increases,
    the dimension of the random effects vector grows linearly or
    quadratically, creating a identifiablity issue for numerical
    integration or MCMC sampling.
2.  **Rigidity**: Models often assume Gaussian error terms and linear
    associations, limiting their applicability to skewed data or
    non-linear biological mechanisms.

`JoiNMe` addresses these challenges through a structural decomposition
of the linear predictor and a flexible distributional regression
framework. This vignette outlines the statistical theory governing the
package.

## The Multivariate Longitudinal Submodel

Let \mathcal{D}\_i = \\y\_{id}(t\_{ij}) : j = 1, \dots, n\_{id}, d = 1,
\dots, D\\ denote the observed longitudinal data for subject i, where
y\_{id}(t) is the response for marker d at time t. We assume generalised
linear mixed models (GLMMs) for each marker, linked via a shared latent
structure.

### Structural Decomposition of the Linear Predictor

The linear predictor \eta\_{id}(t) for subject i and marker d at time t
is decomposed into fixed effects, shared random effects, and
marker-specific deviations:

\eta\_{id}(t) = \mathbf{x}\_{id}(t)^\top \boldsymbol{\beta} +
\mathbf{z}^{(id)}\_{i}(t)^\top \mathbf{u}\_i +
\mathbf{z}^{(mk)}\_{d}(t)^\top \mathbf{v}\_d +
\mathbf{z}^{(idm)}\_{id}(t)^\top \mathbf{w}\_{id}

where:

- **Fixed Effects** (\boldsymbol{\beta}): Capture population-level
  trends across covariates \mathbf{x}.
- **Shared Subject Effects** (\mathbf{u}\_i \sim N(\mathbf{0},
  \Sigma_u)): Represents the underlying “frailty” or shared process
  common to all markers for subject i (e.g., general health status).
- **Marker-Specific Random Effects** (\mathbf{v}\_d \sim N(\mathbf{0},
  \Sigma_v)): Captures heterogeneity in marker trajectories common
  across subjects (e.g., marker-specific baseline levels or rates of
  change).
- **Marker-by-Subject Interactions** (\mathbf{w}\_{id} \sim
  N(\mathbf{0}, \Sigma\_{w,i})): Captures the specific deviation of
  marker d for subject i from the shared trends.

This hierarchical structure allows information sharing across markers
(shrinking estimates towards the shared process \mathbf{u}\_i) while
maintaining flexibility for individual marker dynamics.

### Distributional Assumptions

The observed response y\_{id}(t) is generated from a distribution
\mathcal{F} with location parameter linked to \eta\_{id}(t) and
potentially other ancillary parameters \boldsymbol{\psi}\_{d}:

y\_{id}(t) \sim \mathcal{F}\_d(\cdot \| \eta\_{id}(t),
\boldsymbol{\psi}\_{d})

`JoiNMe` supports a wide range of families \mathcal{F}\_d, including:

- **Gaussian / Student’s t**: For continuous, potentially heavy-tailed
  data.
- **Beta**: For bounded (0,1) proportions.
- **Bernoulli / Binomial**: For binary or count-of-success data.
- **Poisson / Negative Binomial**: For count data.
- **Double Exponential / Skew Double Exponential**: For continuous
  outcomes with Laplace tails and, when required, asymmetric conditional
  quantiles.
- **Cumulative Logit**: For ordinal outcomes.

Furthermore, distributional parameters (for example, residual scale
\sigma) can themselves be modelled as functions of covariates
(distributional regression): \log(\sigma\_{id}(t)) =
\mathbf{x}\_{\sigma, id}(t)^\top \boldsymbol{\beta}\_\sigma

`joinem` also supports **family-scoped** distributional regressions
using bracket syntax on the left-hand side, e.g.

- `sigma[family=student_t] ~ 1 + time`
- `sigma[family=gaussian] ~ 1 + x1`

These scopes are estimated as separate coefficient blocks and only
contribute to rows of the scoped family. Without brackets (for example,
`sigma ~ ...`), the formula applies to all rows where the parameter
exists.

### Fixed skew-Laplace quantiles

For the skew double exponential distribution, JoiNMe uses a quantile
parameterisation. Given location \mu, scale \sigma\>0, and 0\<\tau\<1,
its density is

f(y\mid\mu,\sigma,\tau) = \frac{2\tau(1-\tau)}{\sigma} \begin{cases}
\exp\\-2(1-\tau)\|y-\mu\|/\sigma\\, & y\<\mu,\\
\exp\\-2\tau\|y-\mu\|/\sigma\\, & y\geq\mu. \end{cases}

Thus \mu is the conditional \tau-quantile, and \tau=0.5 recovers vthe
symmetric Laplace distribution.

A marker-specific fixed value is declared as
`jm_family("skew_laplace", tau = 0.8)`. If it is omitted, `tau` is
estimated through a `tau ~ ...` distributional regression when present,
or as a family-level constant otherwise. Fixed and estimated
skew-Laplace markers may coexist, and posterior prediction selects the
appropriate value by marker.

## Covariance regression

JoiNMe models the marker-by-subject covariance factor as

L_i=\operatorname{diag}(\boldsymbol\sigma_i)K_i,

where \boldsymbol\sigma_i contains marginal standard deviations and K_i
is a lower-triangular Cholesky correlation factor with unit row norms.
The two parts may depend on different subject-level predictors:

``` r

formulaVCov = list(
  sd = ~ treatment + age,
  corr = ~ baseline_score
)
```

A single formula remains shorthand for applying the same regression to
both parts. Component intercepts are represented explicitly, so formula
intercept columns are removed from both model matrices.

For standard-deviation coordinate r,

\eta^{sd}\_{ir} =\alpha^{sd}\_r +\boldsymbol
x\_{sd,i}^\top\boldsymbol\beta^{sd}\_r +\lambda\_{m(r)}z\_{L,i,m(r)},
\qquad \sigma\_{ir}=g\_{sd}(\eta^{sd}\_{ir}),

where g\_{sd} is softplus by default and may instead be exponential. For
an off-diagonal coordinate h=(r,c), c\<r,

\eta^{K}\_{ih} =\alpha^{K}\_h +\boldsymbol
x\_{corr,i}^\top\boldsymbol\beta^{K}\_h +\lambda\_{m(h)}z\_{L,i,m(h)},
\qquad p\_{i,rc}=\tanh(\eta^{K}\_{ih}).

The partial correlations are converted sequentially to the Cholesky row,
K\_{i,rc}=p\_{i,rc}\prod\_{h\<c}\sqrt{1-p\_{i,rh}^2} for c\<r and
K\_{i,rr}=\prod\_{h\<r}\sqrt{1-p\_{i,rh}^2}. Hence SD and dependence
predictors are estimated in separate coefficient blocks and together
determine
\Sigma_i=\operatorname{diag}(\boldsymbol\sigma_i)K_iK_i^\top\operatorname{diag}(\boldsymbol\sigma_i).

The lower-triangle coordinate m retains one non-negative scalar loading
\lambda_m. Thus
\boldsymbol\lambda_L=(\lambda_1,\ldots,\lambda\_{M\_{\mathrm{cov}}})^\top
is a vector and its action is
\operatorname{diag}(\boldsymbol\lambda_L)\boldsymbol z\_{L,i}, not a
dense loading matrix. In the ordinary model, \lambda_m^2 is the
conditional residual variance of the corresponding *untransformed*
covariance predictor. It is not directly a standard deviation,
correlation, or covariance in L_iL_i^\top.

`jm_priors(vcov = list(sd = ..., corr = ...))` assigns independent
intercept and slope prior families to the two blocks. Their formula
dimensions and packed orders are shared by fitting, dynamic prediction,
and both simulation entry points. The dedicated covariance-regression
vignette gives the complete recursion and recovery-study syntax.

## The Survival Submodel

The risk of the terminal event is modelled using a proportional hazards
structure. The observation may identify an exact event time, establish
survival beyond a right-censoring time, establish failure before a
left-censoring time, or locate the event in an inspection interval.

The hazard function \lambda_i(t) at time t is given by:

\lambda_i(t) = \lambda_0(t) \exp\left( \mathbf{w}\_{surv, i}^\top
\boldsymbol{\gamma} + \sum\_{k=1}^K f_k(\mathcal{H}\_{i}(t),
\boldsymbol{\alpha}\_k) \right)

### Baseline Hazard \lambda_0(t)

The baseline hazard is modelled flexibly using B-splines (or P-splines):
\log \lambda_0(t) = \sum\_{l=1}^L \phi_l B_l(t) where B_l(t) are
B-spline basis functions.

### Censoring contributions

Let H_i(t)=\int_0^t\lambda_i(u)\\du and S_i(t)=\exp\\-H_i(t)\\.
Conditional on the longitudinal random effects, the event likelihood
contribution is

\mathcal L\_{S,i}= \begin{cases} S_i(L_i), & T_i\>L_i \quad\text{(right
censored)},\\ \lambda_i(t_i)S_i(t_i), & T_i=t_i \quad\text{(exact)},\\
1-S_i(R_i), & T_i\le R_i \quad\text{(left censored)},\\
S_i(L_i)-S_i(R_i), & L_i\<T_i\le R_i \quad\text{(interval censored)}.
\end{cases}

Thus JoiNMe evaluates a full survival likelihood, rather than a Cox
partial likelihood. For `Surv(lower, upper, type = "interval2")`, the
event data contain one row per subject. A proper interval (L_i,R_i\]
contributes known survival to L_i followed by failure within (L_i,R_i\].
The two log contributions are

-H_i(L_i) +\log\left\[1-\exp\\-\[H_i(R_i)-H_i(L_i)\]\\\right\]
=\log\\S_i(L_i)-S_i(R_i)\\.

The lower inspection time is therefore not treated as delayed entry.
Event covariates supplied on the subject’s interval2 row are held over
the represented risk time, whilst the latent longitudinal association
remains time-varying. Left and interval censoring presently describe one
event type; an unobserved cause cannot be assigned a cause-specific
failure density.

### Association Structure

The term \mathcal{H}\_{i}(t) = \\ \eta\_{id}(s) : s \le t, d=1\dots D \\
represents the history of the longitudinal processes. `JoiNMe`
summarises this history using a small set of interpretable features:

- **Current value (CV)**: the value of the latent longitudinal linear
  predictor at time t.
- **Current slope (CS)**: the derivative of the latent longitudinal
  linear predictor at time t.
- **Correlation association (CORR)**: a subject-specific, time-constant
  summar of the off-diagonal entries of the Cholesky-correlation factor
  `K_i` from the marker-by-subject covariance model.
- **Variance-covariance association (VCOV)**: a
  subject-specific,time-constant summary of those same off-diagonal
  `K_i` entries together with the subject-specific standard deviations.
  This channel therefore avoids double-counting information that would
  appear again in `L_i L_i^\top`.

CV is the fitted latent trajectory at the event time and CS is its
derivative per unit of original study time. These summaries can be
transformed before entering the hazard, as described below.

For covariance-style channels (`corr`, `vcov`), the declared
`assoc$slope` prior applies directly to the effective coefficient
entering the log hazard. The coefficient is not multiplied by a second
random association scale. Normal, fixed-degrees-of-freedom Student-t,
Laplace, and regularised horseshoe priors are selected independently
through `jm_priors(assoc = list(slope = ...))`.

The covariance-regression layer itself uses an identified
parameterisation for the subject-specific Cholesky factor: each
lower-triangular entry has its own latent perturbation z\_{im} and a
non-negative loading \lambda_m, so \eta\_{im} = \alpha_m + x_i^\top
\beta_m + \lambda_m z\_{im} with \lambda_m \ge 0 and z\_{im} \sim
\mathcal{N}(0, 1).

#### Marker-weighted contributions (CV and CS)

When multiple markers are present, `JoiNMe` forms marker-average
association components using weights \omega_d (one per marker). These
weights apply **both** to the current value (CV) and current slope (CS)
summaries. Let \eta\_{id}^{(mk)}(t) denote the marker-specific deviation
component of the linear predictor for marker d. Then the marker-weighted
summaries are:

\eta\_{i}^{(mk)}(t) = \frac{1}{D}\sum\_{d=1}^{D} \omega_d
\eta\_{id}^{(mk)}(t), \quad \frac{d}{dt}\eta\_{i}^{(mk)}(t) =
\frac{1}{D}\sum\_{d=1}^{D} \omega_d \frac{d}{dt}\eta\_{id}^{(mk)}(t).

These enter the hazard through marker-specific association terms (for
example, `cv_marker`, `cs_marker`) with transform-before-averaging
semantics for marker-aggregated terms.

To make the weighting explicit, let \mathbf{v}\_d denote marker-level
random effects and \mathbf{w}\_{id} denote marker-by-subject random
effects. The marker-weighted summaries are formed as weighted means
across markers with weights \omega_d. The prior declaration supplies
known offsets directly.

For a stochastic family, let w^{raw}\_{sd} = w^{(0)}\_{sd} +
\mu\_{\omega,s} + \delta\_{sd}, where the common location borrows
information across markers and \delta\_{sd} is a centred, unit-scale
departure. With `family = "constant"` or `"none"`, both fitted terms are
absent and w^{raw}\_{sd} = w^{(0)}\_{sd}. Set \omega\_{sd} =
w^{raw}\_{sd}. The orresponding weighted means are computed as

\bar{\mathbf{v}} = \sum\_{d=1}^D \omega_d \mathbf{v}\_d, \qquad
\bar{\mathbf{w}}\_i = \sum\_{d=1}^D \omega_d \mathbf{w}\_{id}.

These weighted means are then inserted into the marker components of the
current value and slope calculations.

#### Priors on marker contributions

Marker weights and CV/CS association coefficients are regularized to
avoid overfitting when D is large.

- **Marker weights**: when estimated, `joinme` uses
  w\_{sd}=w^{(0)}\_{sd}+\mu\_{\omega,s}+z\_{sd}. The fitted common
  location uses `marker_weights$intercept`; centred unit-scale
  coordinates use `marker_weights$family`. Student-t degrees of freedom
  are shared across all weight sets, and no additional scale multiplies
  the marker departure. We may add some known contribution w^{(0)}\_{sd}
  by declaring `marker_weights$offset`. To disable the fitting process
  and only use these known weights, we can set `family = "constant"` or
  `family = "none"`.
- **Caveat: CV/CS total and marker coefficients**: role-specific priors
  are transformed from standardised raw parameters, with a single sign
  anchor for identifiability when marker-weight channels are active.
  Specifically, only the first active latent in the order `cv_total`,
  `cs_total`, `cv_marker`, `cs_marker` is constrained non-negative (via
  absolute-value mapping); remaining active channels are unconstrained.

This keeps marker effects signed (positive or negative) while improving
sampling geometry for hierarchical association parameters.

`joinme` also supports non-linear transformations of these terms via the
`transforms` argument. Available transformation modes include identity,
a user-defined mathematical expression, monotone I-splines, and ordered
piecewise-linear rules whose knot ordinates are learned from simplex
increments.

Within a functional expression, the standard Normal distribution
function \Phi(x) and its inverse \Phi^{-1}(p) are not interchangeable. A
probit longitudinal link uses the latter as its forward link and the
former as its inverse link.

For `corr` and `vcov`, the transformed summary is time-constant for a
fixed subject because it depends only on that subject’s
covariance-regression Cholesky factor. The implemented likelihood
therefore computes the transformed covariance-style scalar once per
subject, then reuses that same value in both:

- the exact event-time log-hazard \log \lambda_i(T_i),
- the Gauss-Kronrod node evaluations used inside \int_0^{T_i}
  \lambda_i(s) ds.

## Bayesian Inference and Estimation

The joint likelihood is the product of the longitudinal likelihood
\mathcal{L}\_{Long} and the censoring-appropriate survival likelihood
\mathcal{L}\_{Surv}:

p(\theta, \mathbf{u}, \mathbf{v}, \mathbf{w} \| \mathbf{y}, \mathcal
O_T) \propto \left\[ \prod\_{i=1}^N \prod\_{d=1}^D
\prod\_{j=1}^{n\_{id}} p(y\_{idj} \| \mathbf{u}\_i, \mathbf{v}\_d,
\mathbf{w}\_{id}) \right\] \times \left\[ \prod\_{i=1}^N \mathcal
L\_{S,i}(\mathcal O\_{T_i}) \right\] \times p(\text{priors})

### Priors

Default priors are weakly informative:

- **Ordinary regressions**: global intercepts default to centred
  Student-t_6 priors with scale 2, while slopes default to centred
  Student-t_6 priors with scale 1. `longitudinal`, `survival`, `vcov`,
  `functional`, and each `formulaDist` left-hand side may replace
  intercept and slope roles independently. A distributional prior may
  repeat a family scope from `formulaDist`, for example
  `` `sigma[family='student']` ``. A marker scope is also accepted when
  the marker uniquely identifies that family-scoped coefficient vector;
  shared family coefficients must use the family name.
- **Association coefficients**: CV/CS total, mean, and marker
  coefficients use the independently declared `assoc$slope`
  transformation. For marker-weight-involving channels (`cv_total`,
  `cs_total`, `cv_marker`, `cs_marker`), only the first active latent is
  sign-anchored non-negative to improve identifiability; other channels
  remain free to be positive or negative.
- **Marker weights**: when estimated, every weight set has a fitted
  common location and marker-specific additive departures around
  user-supplied offsets. The common location uses the
  `marker_weights$intercept` prior and borrows information from all
  markers in its set. The standardised departures use the family named
  by `marker_weights$family`; the family name `"student_t"` additionally
  fits one shared degrees-of-freedom value with
  \nu-2\sim\operatorname{Gamma}(2,0.1) using departures from every set.
- **Random-effect scales**: subject and marker standard deviations use
  exponential priors. Marker-by-subject covariance scales arise from
  their separate SD regression and non-negative residual loadings.
- **Random-effect correlations**: correlation matrices use an
  LKJ-Cholesky prior with concentration parameter \eta, allowing mild
  regularisation toward identity while remaining flexible.
- **Baseline hazard**: the baseline level has an intercept prior. Spline
  coefficients use a second-difference prior to favour a smooth baseline
  hazard without penalising its overall level.
- **Penalised transformations**: when `type = "ispline_penalised"`, the
  I-spline coefficients use a smoothness penalty (second differences) to
  stabilise nonlinear association transforms in sparse time regions.
- **Family-specific parameters**: scale parameters use exponential
  priors, Student-t degrees of freedom use gamma priors, and skewness
  parameters on (0,1) use beta priors, matching the parameter support.

## Competing risks and dynamic prediction

`joinme` supports competing risks through the event-type indicator
(`event_type`) or, when available, a multi-state pair of variables
(`stage_from`, `stage_to`) that are collapsed into a single event-type
factor. However, it does not support directly muti-state model like
`JMbayes2` or `INLAjoint`. The competing-risk likelihood is a product of
cause-specific hazards, and the default prediction methods return
overall survival across all causes. Cause-specific prediction requires a
separate post-processing step.

The baseline hazard and hazard covariates are then indexed by event
type, so th instantaneous risk is a cause-specific hazard.

In the likelihood, the cumulative hazard integrates the sum of
cause-specific hazards, yielding an overall survival function. In
prediction, the dynamic survival and cumulative hazard outputs
correspond to the overall survival across all causes, conditional on
being event-free up to the landmark time. Cause- specific prediction
would require a separate post-processing step that is not returned by
the default prediction methods.

## References

1.  Rizopoulos D. Joint Models for Longitudinal and Time-to-Event Data:
    With Applications in R. Chapman and Hall/CRC; 2012.
2.  Carpenter B, Gelman A, Hoffman MD, et al. Stan: A probabilistic
    programming language. J Stat Softw. 2017;76(1).
3.  Brilleman SL, Elci EM, Novik JB, Wolfe R. Joint longitudinal and
    time-to-event models for health data. Stat Methods Med Res. 2018.
