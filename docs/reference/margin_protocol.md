# Fitted marginal distribution protocol

Fitted margins provide distribution evaluation and model metadata
through a small set of S3 generics. Methods for `dmargin()`,
`pmargin()`, and `qmargin()` dispatch on `margin`, their second
argument. `margin_info()` dispatches on its only argument.

## Usage

``` r
dmargin(x, margin)

pmargin(x, margin)

qmargin(p, margin)

margin_info(object)
```

## Arguments

  - x:
    
    vector of evaluation points.

  - margin:
    
    a fitted marginal distribution object.

  - p:
    
    vector of probabilities.

  - object:
    
    a fitted margin or margin-family specification.

## Value

`dmargin()`, `pmargin()`, and `qmargin()` return numeric vectors.
`margin_info()` returns a list of model information.

## Details

A fitted margin class must implement all four generics documented here.
Evaluation methods must be vectorized. `margin_info()` returns a list
with entries `family_name`, `type`, `support`, `npars`, and `loglik`.
`type` is `"c"` for a continuous distribution, `"d"` for an
integer-supported distribution, or `"zi"` for a continuous distribution
with an atom at zero. `support` contains the mathematical lower and
upper bounds. `npars` must be non-negative or `NA` when unavailable;
`loglik` may be `NA` for a fixed margin.

The type determines left limits throughout rvinecopulib. They are `F(x)`
for continuous margins, `F(x - 1)` for integer-supported margins, and
`F(0) - P(X = 0)` at zero for zero-inflated margins. Consequently,
`dmargin(0, margin)` must return the atom's mass for a zero-inflated
margin.
