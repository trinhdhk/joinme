# Resolve recursive Stan include dependencies

Resolve recursive Stan include dependencies

## Usage

``` r
.stan_dependency_files(stan_file, visited = character())
```

## Arguments

- stan_file:

  Path to root Stan file.

- visited:

  Internal recursion guard.

## Value

Character vector of dependency file paths (including the root file).
