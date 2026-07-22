% Manual page maintained directly until roxygen is regenerated.
\name{longitudinal_plot}
\alias{longitudinal_plot}
\alias{survival_plot}
\alias{cumhaz_plot}
\alias{association_plot}
\alias{diagnostic_plot}
\alias{mcmc_plot}
\title{Convenience plot helpers for JoiNMe objects}
\usage{
longitudinal_plot(object, ...)

survival_plot(object, ...)

cumhaz_plot(object, ...)

association_plot(object, ...)

diagnostic_plot(
  object,
  type = c("rhat", "ess_bulk", "ess_tail", "mcse_mean", "mcse_sd", "running_mean",
    "running_quantile"),
  ...
)

mcmc_plot(
  object,
  pars = NA,
  type = c("intervals", "areas", "dens", "dens_overlay", "hist", "trace", "violin",
    "acf", "rhat", "neff"),
  variable = NULL,
  regex = FALSE,
  fixed = FALSE,
  draws = NULL,
  seed = 1,
  ...
)
}
\arguments{
\item{object}{A \code{JoiNMeFit}, \code{JoiNMeFit}, \code{JoiNMeDynPred}, or
\code{PredJoiNMeFit} object as appropriate for the requested helper.}

\item{type}{The requested diagnostic or bayesplot geometry.}

\item{pars}{Deprecated alias of \code{variable} for \code{mcmc_plot()}.}

\item{variable}{Optional renamed posterior variable names for
\code{mcmc_plot()}.}

\item{regex}{Logical; whether \code{variable} should be treated as a regular
expression.}

\item{fixed}{Deprecated logical alias controlling exact versus regular-expression
matching when \code{pars} is supplied.}

\item{draws}{Optional number of posterior draws to keep.}

\item{seed}{Integer seed used when subsetting posterior draws.}

\item{...}{Additional arguments forwarded to the corresponding \code{plot()} or
\code{bayesplot::mcmc_*()} implementation.}
}
\value{
For the plot helpers, a \code{ggplot} object, a combined plot object, or a
named list of plots depending on the requested object and plot family.
}
\description{
Convenience wrappers around the main \code{plot()} methods for JoiNMe fit and
dynamic prediction objects.

\code{longitudinal_plot()}, \code{survival_plot()}, and \code{cumhaz_plot()}
route to the corresponding fitted or predicted trajectory plot.

\code{association_plot()} exposes the fitted association-curve display for
\code{JoiNMeFit} objects.

\code{diagnostic_plot()} exposes the scalar and running sampler diagnostics for
\code{JoiNMeFit} objects.

\code{mcmc_plot()} applies bayesplot MCMC geometries to the renamed posterior
draws returned by \code{draws()}, which means the displayed parameter labels are
the same user-facing names used by summaries and diagnostics.
}
