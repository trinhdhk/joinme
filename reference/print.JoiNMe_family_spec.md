# Method for printing a `JoiNMe_family_spec` object

Prints the family, forward link, and inverse-link expression. The
inverse link is reconstructed from its stored functional bytecode so
that a family specification displays mathematical R syntax rather than
the underlying instruction and constant vectors.

## Usage

``` r
# S3 method for class 'JoiNMe_family_spec'
print(x, ...)
```

## Arguments

- x:

  An object of class `JoiNMe_family_spec`.

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns the input object `x`.
