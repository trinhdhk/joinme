% Manual page maintained directly until roxygen is regenerated.
\name{draws}
\alias{draws}
\title{Extract renamed posterior draws from JoiNMe objects}
\usage{
draws(object, ...)
}
\arguments{
\item{object}{A \code{JoiNMeFit}, \code{JoiNMeFit}, \code{JoiNMeDynPred}, or
\code{PredJoiNMeFit} object.}

\item{...}{Method-specific arguments. For fit and dynamic-prediction objects,
the supported arguments are \code{variables}, \code{regex}, \code{draws},
\code{seed}, and \code{format}.}
}
\value{
A posterior draw object in the requested format. Supported formats are
\code{"draws_array"}, \code{"draws_matrix"}, and \code{"draws_df"}.
}
\description{
Extract posterior draws using user-facing parameter names rather than raw Stan
variable names.

Compared with \code{extract()}, \code{draws()} is the convenience layer for
downstream posterior workflows. It returns a single \pkg{posterior}-compatible
object that is ready for variable filtering, summary, and visualisation.

For \code{JoiNMeFit} objects, \code{draws()} relabels fixed effects,
survival coefficients, association parameters, and supported distributional
parameters using the same naming logic consumed by summaries and diagnostics.

For \code{JoiNMeDynPred} objects, \code{draws()} flattens the stored prediction
draw blocks into posterior-compatible formats so the results can be used with
packages such as \pkg{posterior} and \pkg{bayesplot}.

Use \code{draws()} when you want one renamed posterior object. Use
\code{extract()} when you want component-specific extraction together with the
associated mapping metadata.
}
