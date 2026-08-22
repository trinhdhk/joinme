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
#> Joint Nested Mixed-effects Longitudinal Family Specification:
#>   Family: poisson
#>   Link: log
#>   Inverse link: exp(x)
jm_family("poisson", link = ~ log(x))
#> Joint Nested Mixed-effects Longitudinal Family Specification:
#>   Family: poisson
#>   Link: log
#>   Inverse link: exp(x)
jm_family("bernoulli", inv_link = ~ inv_logit(x))
#> Joint Nested Mixed-effects Longitudinal Family Specification:
#>   Family: bernoulli
#>   Link: logit
#>   Inverse link: inv_logit(x)
jm_family("bernoulli", link = ~ inv_Phi(x))
#> Joint Nested Mixed-effects Longitudinal Family Specification:
#>   Family: bernoulli
#>   Link: probit
#>   Inverse link: Phi(x)
jm_family("bernoulli", inv_link = ~ Phi(x))
#> Joint Nested Mixed-effects Longitudinal Family Specification:
#>   Family: bernoulli
#>   Link: probit
#>   Inverse link: Phi(x)
jm_family("skew_laplace", tau = 0.8)
#> Joint Nested Mixed-effects Longitudinal Family Specification:
#>   Family: skew_double_exponential
#>   Link: identity
#>   Inverse link: x
```
