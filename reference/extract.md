# Extract posterior draws from JoiNMe objects

S3 generic to extract component-specific posterior payloads from
`JoiNMeFit` and `JoiNMeDynPred` objects.

Compared with
[`draws()`](https://trinhdhk.github.io/joinme/reference/draws.md),
`extract()` provides lower-level interface. It works component by
component and returns the metadata needed to understand how a requested
summary term maps back to the stored Stan variables or prediction draw
blocks.

Use `extract()` when you need a specific model component, the
corresponding `term_map`, or a specialised payload such as
`what = "association_plot"`. Use
[`draws()`](https://trinhdhk.github.io/joinme/reference/draws.md) when
you want a single posterior object ready for `posterior` or `bayesplot`
workflows.

## Usage

``` r
extract(object, ...)
```

## Arguments

- object:

  A supported JoiNMe object.

- ...:

  Additional method-specific arguments.

## See also

[`draws()`](https://trinhdhk.github.io/joinme/reference/draws.md)
