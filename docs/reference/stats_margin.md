# Create a fixed margin from a stats distribution

`stats_margin()` adapts the univariate distributions documented in
[stats::Distributions](https://rdrr.io/r/stats/Distributions.html) to
the fitted-margin protocol. Distribution parameters are supplied through
`...`; the resulting margin reports the distribution's parameter count
but has no fitted log-likelihood.

## Usage

``` r
stats_margin(family, ...)
```

## Arguments

  - family:
    
    distribution suffix such as `"norm"`, `"pois"`, or `"binom"`.

  - ...:
    
    parameters passed to the distribution's density, CDF, and quantile
    functions.

## Value

A fixed margin implementing the [margin
protocol](https://vinecopulib.github.io/rvinecopulib/reference/margin_protocol.md).

## Examples

``` r
normal <- stats_margin("norm", mean = 1, sd = 2)
poisson <- stats_margin("pois", lambda = 3)
margin_info(normal)
#> $family_name
#> [1] "norm"
#> 
#> $type
#> [1] "c"
#> 
#> $support
#> [1] -Inf  Inf
#> 
#> $npars
#> [1] 2
#> 
#> $loglik
#> [1] NA
#> 
margin_info(poisson)
#> $family_name
#> [1] "pois"
#> 
#> $type
#> [1] "d"
#> 
#> $support
#> [1]   0 Inf
#> 
#> $npars
#> [1] 1
#> 
#> $loglik
#> [1] NA
#> 
```
