# JoiNMe longitudinal family specification

Creates a family specification object for `joinme(..., families = ...)`
with an optional custom longitudinal link (or inverse-link) per marker
family.

This helper is intentionally small and explicit: it stores only
canonical family name plus a supported inverse-link choice that can be
passed to Stan as compact integer link codes.

Supported links / inverse-links:

- `"identity"` : \\g^{-1}(x)=x\\

- `"log"` : \\g^{-1}(x)=\log(x)\\

- `"logit"` : \\g^{-1}(x)=\operatorname{logit}^{-1}(x)\\

- `"probit"` : \\g^{-1}(x)=\Phi(x)\\

- `"exp"` : \\g^{-1}(x)=\exp(x)\\

## Usage

``` r
joinme_family(name, link = NULL, inv_link = NULL)

jm_family(name, link = NULL, inv_link = NULL)
```

## Arguments

- name:

  Character scalar family name, e.g. `"student_t"`, `"bernoulli"`.

- link:

  Optional character scalar naming the inverse-link.

- inv_link:

  Optional one-sided formula `~ ...` or expression string using `x`
  (e.g. `~ exp(x)`, `~ inv_logit(x)`, `~ probit(x)`).

## Value

Object of class `"JoiNMe_family_spec"`.
