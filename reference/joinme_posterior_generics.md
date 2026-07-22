# Posterior and log-likelihood generics for JoiNMe objects

Package-local generics for posterior predictive extraction and pointwise
log-likelihood evaluation.

These generics mirror the common Bayesian-modeling API used by packages
such as brms and rstantools, but they are defined inside JoiNMe so the
package no longer requires rstantools at load time.

## Usage

``` r
posterior_predict(object, ...)

posterior_epred(object, ...)

posterior_linpred(object, ...)

log_lik(object, ...)
```

## Arguments

- object:

  An object supporting the requested posterior or log-likelihood
  extraction.

- ...:

  Additional arguments passed to the class-specific method.

## Value

Class-specific posterior prediction or log-likelihood output.
