# Precompile CmdStanR models

Manually precompiles the packaged Stan models using CmdStanR. This is
useful when you want to force compilation ahead of time, rather than
waiting for the first call with engine = "cmdstanr".

## Usage

``` r
precompile_cmdstanr_models(force_recompile = TRUE, cleanup = TRUE)
```

## Arguments

- force_recompile:

  Logical; recompile even if cached.

- cleanup:

  Logical; whether to clean up old cached executables before compiling.

## Value

Invisibly returns a named list of compiled CmdStanR models.
