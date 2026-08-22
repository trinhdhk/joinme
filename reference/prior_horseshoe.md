# Declare a regularised horseshoe prior

Declares the regularised horseshoe hierarchy used by brms and rstanarm.
The coefficient is conditionally Normal with local and global
half-Student-t scales and a finite Student-t slab. The slab regularises
very large signals, whilst the local scales permit coefficient-specific
escape from regularisation.

## Usage

``` r
prior_horseshoe(
  df = 1,
  global_df = 1,
  global_scale = 1,
  slab_df = 4,
  slab_scale = 2
)
```

## Arguments

- df:

  Positive fixed degrees of freedom for local regularisation scales.

- global_df:

  Positive fixed degrees of freedom for the global scale.

- global_scale:

  Positive global scale. Smaller values express stronger prior sparsity.

- slab_df:

  Positive fixed degrees of freedom for the regularising slab.

- slab_scale:

  Positive scale of the regularising slab.

## Value

A `joinme_prior_spec` object for use inside
[`jm_priors()`](https://trinhdhk.github.io/joinme/reference/joinme_priors.md).
