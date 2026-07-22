% Manual page maintained directly until roxygen is regenerated.
\name{posterior_predict}
\alias{posterior_predict}
\alias{posterior_epred}
\alias{posterior_linpred}
\alias{log_lik}
\title{Posterior and log-likelihood generics for JoiNMe objects}
\usage{
posterior_predict(object, ...)

posterior_epred(object, ...)

posterior_linpred(object, ...)

log_lik(object, ...)
}
\arguments{
\item{object}{An object supporting the requested posterior or log-likelihood extraction.}

\item{...}{Additional arguments passed to the class-specific method.}
}
\value{
Class-specific posterior prediction or log-likelihood output.
}
\description{
Package-local generics for posterior predictive extraction and pointwise
log-likelihood evaluation.

These generics mirror the common Bayesian-modeling API used by packages such as
\pkg{brms} and \pkg{rstantools}, but they are defined inside \pkg{JoiNMe} so the
package no longer requires \pkg{rstantools} at load time.
}
