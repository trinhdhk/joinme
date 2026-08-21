# Extract stored posterior draw components from dynamic prediction objects

Extracts stored prediction draw blocks without flattening them first.

This is the structured companion to
[`posterior_draws.JoiNMeDynPred()`](https://trinhdhk.github.io/joinme/reference/posterior_draws.md).
Use
[`extract()`](https://trinhdhk.github.io/joinme/reference/extract.md)
when you want to keep the original prediction block semantics
(`longitudinal`, `survival`, `cumhaz`, random effects, and scale/id
filters). Use
[`posterior_draws()`](https://trinhdhk.github.io/joinme/reference/posterior_draws.md)
when you want those blocks flattened into one `posterior`-compatible
draw object with composite variable labels.

## Usage

``` r
# S3 method for class 'JoiNMeDynPred'
extract(
  object,
  what = c("longitudinal", "longitudinal_fitted", "survival", "cumhaz",
    "random_effects_id", "random_effects_marker_id"),
  id = NULL,
  scale = NULL,
  ...
)
```

## Arguments

- object:

  A `JoiNMeDynPred` object.

- what:

  Draw block selector: `"longitudinal"`, `"longitudinal_fitted"`,
  `"survival"`, `"cumhaz"`, `"random_effects_id"`,
  `"random_effects_marker_id"`.

- id:

  Optional character/integer id filter.

- scale:

  Optional scale filter for longitudinal blocks (`epred`, `linpred`,
  `predict`).

## Value

A list with fields:

- `posterior_draws`: numeric matrix or list of matrices

- `meta`: extraction metadata
