# Define a custom margin family

`margin_family()` defines a candidate family for `vine()`. Its fitting
function must accept `x`, `weights`, and `type`, and return an object
implementing the [fitted-margin
protocol](https://vinecopulib.github.io/rvinecopulib/reference/margin_protocol.md).
Additional named fitting arguments can be stored in `fit_args`.

## Usage

``` r
margin_family(fit, family_name = "custom", types = "c", fit_args = list())
```

## Arguments

  - fit:
    
    function accepting `x`, `weights`, and `type`.

  - family\_name:
    
    descriptive family name.

  - types:
    
    supported variable types.

  - fit\_args:
    
    named list of additional arguments passed to `fit`.

## Value

A margin-family specification.

## Examples

``` r
normal_family <- margin_family(
  fit = function(x, weights, type) {
    if (length(weights)) {
      location <- weighted.mean(x, weights)
      scale <- sqrt(weighted.mean((x - location)^2, weights))
    } else {
      location <- mean(x)
      scale <- sd(x)
    }
    margin_dist(
      d = function(y) dnorm(y, location, scale),
      p = function(y) pnorm(y, location, scale),
      q = function(p) qnorm(p, location, scale),
      family = "normal",
      type = type,
      npars = 2,
      loglik = sum(dnorm(x, location, scale, log = TRUE))
    )
  },
  family_name = "normal"
)
```
