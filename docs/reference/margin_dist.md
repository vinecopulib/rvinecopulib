# Create a fitted or fixed custom margin

`margin_dist()` creates a fitted or fixed margin from evaluation
functions and explicit metadata. The functions must each accept one
vector argument and return a numeric vector of the same length.

## Usage

``` r
margin_dist(
  d,
  p,
  q,
  family = "custom",
  type = "c",
  support = c(-Inf, Inf),
  npars = 0,
  loglik = NA_real_
)
```

## Arguments

  - d:
    
    density or probability mass function.

  - p:
    
    distribution function.

  - q:
    
    quantile function.

  - family:
    
    descriptive family name.

  - type:
    
    margin type: `"c"`, `"d"`, or `"zi"`.

  - support:
    
    mathematical support as its lower and upper bounds.

  - npars:
    
    effective number of parameters.

  - loglik:
    
    maximized marginal log-likelihood, or `NA` for a fixed margin.

## Value

A fitted margin implementing the [margin
protocol](https://vinecopulib.github.io/rvinecopulib/reference/margin_protocol.md).

## Examples

``` r
normal_margin <- margin_dist(
  d = function(x) dnorm(x, mean = 1, sd = 2),
  p = function(x) pnorm(x, mean = 1, sd = 2),
  q = function(p) qnorm(p, mean = 1, sd = 2),
  family = "norm",
  support = c(-Inf, Inf),
  npars = 2
)
dmargin(1, normal_margin)
#> [1] 0.1994711
pmargin(1, normal_margin)
#> [1] 0.5
qmargin(0.5, normal_margin)
#> [1] 1
```
