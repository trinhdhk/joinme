# Technical Specification and Implementation Details

## 1 Overview

This document provides the formal technical specification for the
`joinme` package model structure. It details the precise matrix
formulations, transformation logic, and the exact algorithms used by the
underlying Stan engine. Ideally, this document serves as a reference for
advanced methodologists or developers seeking to understand or extend
the package internals.

## 2 Mathematical formulation

The joint model consists of two linked submodels: a multivariate
longitudinal mixed-effects model and a survival model with time-varying
covariates defined by the longitudinal process.

### 2.1 Notation

- $`N`$: Number of subjects ($`i = 1, \dots, N`$).
- $`D`$: Number of longitudinal markers ($`d = 1, \dots, D`$).
- $`n_{id}`$: Number of observations for subject $`i`$ marker $`d`$.
- $`T_i`$: Observed event or censoring time for subject $`i`$.
- $`\delta_i`$: Event indicator (1 if event, 0 if censored).
- $`y_{ijk}`$: The $`k`$-th observation of marker $`d`$ for subject
  $`i`$ at time $`t_{ijk}`$.

### 2.2 Longitudinal Submodel

We assume the observed outcomes $`y_{ijk}`$ follow a distribution from
the exponential family (e.g., Gaussian, Student-t) centred at a linear
predictor $`\eta_{ijk}(t)`$.

``` math
y_{ijk} \sim \mathcal{D}(\eta_{id}(t_{ijk}), \sigma_y, \nu)
```

The linear predictor $`\eta_{id}(t)`$ is decomposed into fixed effects
and random effects components. A key feature of `JoiNMe` is the flexible
random effects structure:

``` math
\eta_{id}(t) = \mathbf{x}_{id}(t)^\top \boldsymbol{\beta} + \mathbf{z}_{id, id}(t)^\top \mathbf{b}_i + \mathbf{z}_{id, mk}(t)^\top \mathbf{u}_d + \mathbf{z}_{id, idm}(t)^\top \mathbf{w}_{id}
```

where:

- $`\mathbf{x}_{id}(t)`$: Fixed effects design vector (time-varying).
- $`\mathbf{b}_i \sim N(\mathbf{0}, \Sigma_b)`$: Subject-specific random
  effects (shared across markers or independent).
- $`\mathbf{u}_d \sim N(\mathbf{0}, \Sigma_u)`$: Marker-specific random
  effects (e.g., marker-specific average intercept).
- $`\mathbf{w}_{id} \sim N(\mathbf{0}, \Sigma_w)`$:
  Subject-marker-specific random effects (the interaction level).

**Optional marker-by-id block.** If the marker block omits the inner
`( ... | id )` term, then $`\mathbf{w}_{id}`$ is removed and the
marker-by-id block has dimension $`Q_{idm} = 0`$. In this case,
covariance-style associations (`corr`, `vcov`) are not available,
because there is no subject-specific marker covariance to model.

This decomposition allows capturing correlations within subjects across
markers, within markers across subjects, and the specific deviation of a
subject’s trajectory for a specific marker.

### 2.3 Survival Submodel

The hazard function $`h_i(t)`$ for subject $`i`$ at time $`t`$ is
modelled as:

``` math
h_i(t) = h_0(t) \exp\left( \mathbf{w}_i^\top \boldsymbol{\gamma} + \sum_{k=1}^K g_k(\mathcal{M}_i(t), \boldsymbol{\alpha}_k) \right)
```

where:

- $`h_0(t)`$: Baseline hazard function modelled via spline basis or a
  user-specified baseline hazard formula.
- $`\mathbf{w}_i`$: Baseline covariates (baseline age, sex, treatment).
- $`\mathcal{M}_i(t) = \{ \eta_{i1}(t), \dots, \eta_{iD}(t), \eta'_{i1}(t), \dots \}`$:
  The collection of longitudinal trajectories and their derivatives
  (slopes) at time $`t`$.
- $`g_k(\cdot)`$: A functional transformation linking the longitudinal
  history to the hazard (the “association structure”).

Common forms for $`g_k(\cdot)`$ in `JoiNMe` are constructed from the
following association summaries, each optionally transformed:

- **Current value (CV)**: the value of the longitudinal linear predictor
  at time $`t`$.
- **Current slope (CS)**: the derivative of the linear predictor at time
  $`t`$.
- **Correlation association (CORR)**: a subject-specific, time-constant
  summary of off-diagonal marker-by-id correlations.
- **Variance-covariance association (VCOV)**: a subject-specific,
  time-constant summary of lower-triangular marker-by-id Cholesky-factor
  features, taken directly from `L`.

In implementation, CV and CS are evaluated at a default set of 15
Gauss-Kronrod nodes on the scaled time interval (configurable via
`quadrature_nodes`), and CS is approximated by a forward finite
difference with step size $`\varepsilon`$.

When multiple markers are present, `JoiNMe` forms marker-average CV and
CS components using weights $`\omega_d`$ (one per marker). These weights
apply to marker summaries after transformation and before association
scaling. When marker-weight estimation is enabled, latent perturbations
are added to base weights and the resulting weights are used directly
before aggregation. Marker weights and association coefficients are
regularized with standard Student-$`t_6`$, Laplace, or Normal priors for
`shrinkage = 0, 1, 2`, respectively, using a positive non-centred
parameterisation $`\alpha = z_{\alpha} \cdot sd_{\alpha}`$.

For covariance-style channels (`corr`, `vcov`), the fitted association
coefficients use the same positive-scale parameterisation: a global
scale (`s_corr` or `s_vcov`) multiplies a latent coefficient vector, and
summaries report the resulting effective coefficients on the hazard
scale.

For covariance regression, each lower-triangular entry has its own
latent perturbation $`z_{im}`$ and a non-negative loading $`\lambda_m`$,
so the identified form is
$`\eta_{im} = \alpha_m + x_i^\top \beta_m + \lambda_m z_{im}`$ with
$`\lambda_m \ge 0`$ and $`z_{im} \sim \mathcal{N}(0, 1)`$.

## 3 Implementation

The model is implemented in Stan, utilising Hamiltonian Monte Carlo
(HMC) for efficient posterior sampling.

### 3.1 Time Scaling

To check numerical stability and ensure spline bases are well-defined,
all internal time calculations are performed on a scaled domain
$`t' \in [0, 1]`$:

``` math
t' = \frac{t}{t_{max}}, \quad t_{max} = \max_i(T_i)
```

Coefficients reported to the user are automatically rescaled to the
original time scale.

### 3.2 Integration Strategy (Gauss-Kronrod)

The survival likelihood requires evaluating the cumulative hazard
$`H_i(T_i) = \int_0^{T_i} h_i(s) ds`$. Since the hazard contains
time-varying random effects $`\eta_{id}(s)`$, this integral has no
closed form.

`JoiNMe` employs **Gauss-Kronrod quadrature** (default 15-point,
configurable) to approximate this integral numerically within the Stan
`log_prob` block. The Stan program uses hardcoded nodes and weights for
the supported node counts, so only the node count is passed in the data.

``` math
\int_0^{T_i} h_i(s) ds \approx \frac{T_i}{2} \sum_{j=1}^{K} w_j h_i \left( \frac{T_i}{2}(x_j + 1) \right)
```

This method provides high accuracy with relatively few function
evaluations compared to Simpson’s rule or trapezoidal integration.

### 3.3 Arbitrary Transformation

`JoiNMe` supports user-specified association transforms in four modes:

- **Identity**: no transformation.
- **Functional bytecode**: a data-defined Reverse Polish Notation (RPN)
  program.
- **Monotone I-spline**: smooth monotone mapping with user-supplied
  knots and coefficients.
- **Piecewise linear**: linear interpolation between user-supplied
  knots.

The functional bytecode evaluator implements the following operations:

- 0: PUSH_X
- 1: PUSH_CONST
- 2: ADD
- 3: SUB
- 4: MUL
- 5: DIV
- 6: LOG
- 7: EXP
- 8: SQRT
- 9: INV_LOGIT
- 10: LOGIT
- 11: RECIPROCAL
- 12: POW
- 13: SIN
- 14: COS
- 15: TAN
- 16: ABS
- 17: SQUARE
- 18: SINH
- 19: COSH
- 20: TANH
- 21: ASINH
- 22: ACOSH
- 23: ATANH
- 24: SOFTPLUS
- 25: CBRT

### 3.4 Parallelisation

The model supports `reduce_sum` for within-chain parallelisation. The
likelihood evaluation is sharded by subject $`i`$.

- **Longitudinal Likelihood**: Calculated independently for each block
  of subjects.
- **Survival Likelihood**: Calculated independently for each block of
  subjects.

The user controls the `grainsize` to balance overhead vs. parallelism.

## 4 Stan Data Block Specification

The generated Stan data block organizes the inputs efficiently:

| Variable | Description |
|:---|:---|
| `n_id` | number of subjects |
| `N` | total longitudinal observations |
| `id`, `marker` | subject and marker indices for each observation |
| `y_real`, `y_int`, `trials` | continuous and discrete outcomes plus binomial trials |
| `X_obs`, `Z_id_obs`, `Z_mk_obs`, `Z_idm_obs` | fixed and random-effect design matrices |
| `S_event`, `d_event`, `event_type` | scaled event times, event indicator, cause type |
| `Bs_event_c`, `Bs_gk_c` | spline bases at event times and GK nodes |
| `X_gk_*`, `Z_*_gk_*` | design matrices at GK nodes (now/fwd for slopes) |
| `tf_mode_*`, `functional_ops_*`, `const_data_*` | association transform configuration |
| `knots_*`, `coeff_*`, `spline_degree_*` | spline-based transform configuration |

## 5 Prediction

Dynamic predictions $`\pi_i(u|t)`$ are efficiently computed using the
same Stan code structure (via `generated quantities` or a separate
`standalone generated quantities` run) to ensure consistency.

``` math
S_i(u|t) = \frac{\exp(-\Lambda_i(u))}{\exp(-\Lambda_i(t))}
```

where $`\Lambda_i(t)`$ uses the same Gauss-Kronrod integrator and
association bytecode as the likelihood.
