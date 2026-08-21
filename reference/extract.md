# Extract posterior draws from JoiNMe objects

S3 generic to extract component-specific posterior results from
`JoiNMeFit` and `JoiNMeDynPred` objects.

Together with
[`posterior_draws()`](https://trinhdhk.github.io/joinme/reference/posterior_draws.md),
`extract()` is the lowest public posterior interface. It works component
by component and returns the metadata needed to understand how a
requested summary term maps back to the stored Stan variables or
prediction draw blocks.

The coefficient hierarchy is deliberately one-directional. `extract()`
owns the draw-level representations selected by
`what = "fixed_effects"`, `"random_effects"`, and `"coefficients"`;
[`posterior_summary()`](https://paulbuerkner.com/brms/reference/posterior_summary.html)
summarises those representations; and
[`fixef()`](https://rdrr.io/pkg/nlme/man/fixed.effects.html),
[`ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html), and
[`coef()`](https://rdrr.io/r/stats/coef.html) are the conventional
high-level entry points. This single path avoids parallel coefficient
APIs with competing semantics.

Use `extract()` when you need a specific model component, the
corresponding `term_map`, or a specialised result such as
`what = "association_plot"`. Use
[`posterior_draws()`](https://trinhdhk.github.io/joinme/reference/posterior_draws.md)
when you want a single posterior object ready for `posterior` or
`bayesplot` workflows.

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

[`posterior_draws()`](https://trinhdhk.github.io/joinme/reference/posterior_draws.md)
