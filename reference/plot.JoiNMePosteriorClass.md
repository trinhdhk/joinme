# Plot posterior class probabilities

Displays posterior mean class probabilities and their central 95%
intervals with a common colour for each latent class. Subject and marker
allocation domains are returned separately because they have different
natural observational units.

## Usage

``` r
# S3 method for class 'JoiNMePosteriorClass'
plot(x, domain = NULL, ...)
```

## Arguments

- x:

  A summary-form `JoiNMePosteriorClass` object.

- domain:

  Optional allocation domain, either `"subject"` or `"marker"`.

- ...:

  Unused.

## Value

A ggplot when one domain is requested or available; otherwise a named
list of ggplots.
