# JoiNMe longitudinal family specification

Creates a marker-specific family specification for
`joinme(..., families = ...)`. A specification may define a custom
longitudinal link (or inverse link) and, for the skew-Laplace family, a
fixed quantile parameter.

Character values supplied through `link` name the *forward* link
\\g(\mu)\\. Formula values supplied through `link` are inverted
symbolically and then compiled by the same functional bytecode evaluator
used for custom association transformations. In contrast, `inv_link`
directly specifies \\g^{-1}(\eta)\\.

The named forward links and their inverse links are:

- `"identity"`:

  \\g(\mu)=\mu\\ and \\g^{-1}(\eta)=\eta\\.

- `"log"`:

  \\g(\mu)=\log(\mu)\\ and \\g^{-1}(\eta)=\exp(\eta)\\.

- `"logit"`:

  \\g(\mu)=\operatorname{logit}(\mu)\\ and
  \\g^{-1}(\eta)=\operatorname{logit}^{-1}(\eta)\\.

- `"probit"`:

  \\g(\mu)=\Phi^{-1}(\mu)\\ and \\g^{-1}(\eta)=\Phi(\eta)\\.

- `"exp"`:

  \\g(\mu)=\exp(\mu)\\ and \\g^{-1}(\eta)=\log(\eta)\\.

## Usage

``` r
joinme_family(name, link = NULL, inv_link = NULL, tau = NULL)

jm_family(name, link = NULL, inv_link = NULL, tau = NULL)
```

## Arguments

- name:

  Character scalar family name, e.g. `"student_t"`, `"bernoulli"`.

- link:

  Optional character scalar naming a forward link, or a one-sided
  formula defining a forward link in `x`, for example `~ log(x + 1)`.
  Formula links must be one-to-one compositions that can be inverted
  algebraically by
  [`invert_transform_expr()`](https://trinhdhk.github.io/joinme/reference/invert_transform_expr.md).

- inv_link:

  Optional one-sided formula `~ ...` or expression string using `x`
  (e.g. `~ exp(x)`, `~ inv_logit(x)`, or `~ Phi(x)`).

- tau:

  Optional fixed quantile/asymmetry parameter for `"skew_laplace"` or
  `"skew_double_exponential"`. It must be a finite scalar strictly
  between zero and one. Under Stan's quantile parameterisation,
  `tau = 0.5` is the symmetric Laplace distribution. Values zero and one
  are excluded because they yield a degenerate, non-normalisable
  limiting distribution. When omitted, `tau` is estimated from a
  distributional regression or as a family-level parameter.

## Value

Object of class `"JoiNMe_family_spec"`.

## Examples

``` r
jm_family("poisson", link = "log")
#> $family
#> [1] "poisson"
#> 
#> $link
#> [1] "log"
#> 
#> $inv_link
#> $inv_link$bytecode
#> [1] 0 7
#> 
#> $inv_link$const_data
#> numeric(0)
#> 
#> $inv_link$n_ops
#> [1] 2
#> 
#> $inv_link$n_const
#> [1] 0
#> 
#> 
#> $tau
#> [1] NA
#> 
#> attr(,"class")
#> [1] "JoiNMe_family_spec"
jm_family("poisson", link = ~ log(x))
#> $family
#> [1] "poisson"
#> 
#> $link
#> [1] "log"
#> 
#> $inv_link
#> $inv_link$bytecode
#> [1] 0 7
#> 
#> $inv_link$const_data
#> numeric(0)
#> 
#> $inv_link$op_iota_intercept_idx
#> [1] 0 0
#> 
#> $inv_link$op_iota_slope_idx
#> [1] 0 0
#> 
#> $inv_link$n_ops
#> [1] 2
#> 
#> $inv_link$n_bytecode
#> [1] 2
#> 
#> $inv_link$n_const
#> [1] 0
#> 
#> $inv_link$n_iota_intercept
#> [1] 0
#> 
#> $inv_link$n_iota_slope
#> [1] 0
#> 
#> 
#> $tau
#> [1] NA
#> 
#> attr(,"class")
#> [1] "JoiNMe_family_spec"
jm_family("bernoulli", inv_link = ~ inv_logit(x))
#> $family
#> [1] "bernoulli"
#> 
#> $link
#> [1] "logit"
#> 
#> $inv_link
#> $inv_link$bytecode
#> [1] 0 9
#> 
#> $inv_link$const_data
#> numeric(0)
#> 
#> $inv_link$op_iota_intercept_idx
#> [1] 0 0
#> 
#> $inv_link$op_iota_slope_idx
#> [1] 0 0
#> 
#> $inv_link$n_ops
#> [1] 2
#> 
#> $inv_link$n_bytecode
#> [1] 2
#> 
#> $inv_link$n_const
#> [1] 0
#> 
#> $inv_link$n_iota_intercept
#> [1] 0
#> 
#> $inv_link$n_iota_slope
#> [1] 0
#> 
#> 
#> $tau
#> [1] NA
#> 
#> attr(,"class")
#> [1] "JoiNMe_family_spec"
jm_family("bernoulli", link = ~ inv_Phi(x))
#> $family
#> [1] "bernoulli"
#> 
#> $link
#> [1] "probit"
#> 
#> $inv_link
#> $inv_link$bytecode
#> [1]  0 26
#> 
#> $inv_link$const_data
#> numeric(0)
#> 
#> $inv_link$op_iota_intercept_idx
#> [1] 0 0
#> 
#> $inv_link$op_iota_slope_idx
#> [1] 0 0
#> 
#> $inv_link$n_ops
#> [1] 2
#> 
#> $inv_link$n_bytecode
#> [1] 2
#> 
#> $inv_link$n_const
#> [1] 0
#> 
#> $inv_link$n_iota_intercept
#> [1] 0
#> 
#> $inv_link$n_iota_slope
#> [1] 0
#> 
#> 
#> $tau
#> [1] NA
#> 
#> attr(,"class")
#> [1] "JoiNMe_family_spec"
jm_family("bernoulli", inv_link = ~ Phi(x))
#> $family
#> [1] "bernoulli"
#> 
#> $link
#> [1] "probit"
#> 
#> $inv_link
#> $inv_link$bytecode
#> [1]  0 26
#> 
#> $inv_link$const_data
#> numeric(0)
#> 
#> $inv_link$op_iota_intercept_idx
#> [1] 0 0
#> 
#> $inv_link$op_iota_slope_idx
#> [1] 0 0
#> 
#> $inv_link$n_ops
#> [1] 2
#> 
#> $inv_link$n_bytecode
#> [1] 2
#> 
#> $inv_link$n_const
#> [1] 0
#> 
#> $inv_link$n_iota_intercept
#> [1] 0
#> 
#> $inv_link$n_iota_slope
#> [1] 0
#> 
#> 
#> $tau
#> [1] NA
#> 
#> attr(,"class")
#> [1] "JoiNMe_family_spec"
jm_family("skew_laplace", tau = 0.8)
#> $family
#> [1] "skew_double_exponential"
#> 
#> $link
#> [1] "identity"
#> 
#> $inv_link
#> $inv_link$bytecode
#> [1] 0
#> 
#> $inv_link$const_data
#> numeric(0)
#> 
#> $inv_link$n_ops
#> [1] 1
#> 
#> $inv_link$n_const
#> [1] 0
#> 
#> 
#> $tau
#> [1] 0.8
#> 
#> attr(,"class")
#> [1] "JoiNMe_family_spec"
```
