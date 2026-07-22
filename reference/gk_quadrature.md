# Gauss-Kronrod quadrature grid helper

Build a quadrature grid on `[0, 1]` for survival integration.

## Usage

``` r
gk_quadrature(nodes = 15L)
```

## Arguments

- nodes:

  Positive integer node count. Allowed values are exactly `7`, `15`,
  `31`, `41`, `51`, and `61`. Defaults to `15`.

## Value

A named list with components: \describe \itemn_gkInteger total node
count. \itemnodesNumeric vector of nodes on `[0, 1]`.
\itemweightsNumeric vector of normalised weights summing to `1`.
\itempanelsAlways `1L` (single fixed rule). \itemruleCharacter rule
label (`"gk7"`, `"gk15"`, `"gk31"`, `"gk41"`, `"gk51"`, or `"gk61"`).

## Details

Only fixed single-panel rules are supported.
