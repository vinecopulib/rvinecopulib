# Margin-family fitting protocol

Margin families describe how a candidate is fitted. Implementations
provide `fit_margin()` and
[`margin_info()`](https://vinecopulib.github.io/rvinecopulib/reference/margin_protocol.md)
methods. Their
[`margin_info()`](https://vinecopulib.github.io/rvinecopulib/reference/margin_protocol.md)
result must contain `family_name` and `types`. Every fitter receives the
observations, observation weights, and requested variable type. It must
either use the weights or reject them.

## Usage

``` r
fit_margin(family, x, weights = numeric(), type = "c")
```

## Arguments

- family:

  a margin-family specification.

- x:

  observed values.

- weights:

  observation weights, or
  [`numeric()`](https://rdrr.io/r/base/numeric.html).

- type:

  requested variable type.

## Value

`fit_margin()` returns a fitted margin.
