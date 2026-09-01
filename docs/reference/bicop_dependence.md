# Dependence measures of a bivariate copula

Computes the four corner tail-dependence coefficients or Blomqvist's
beta for a bivariate copula distribution.

## Usage

``` r
tail_dep(object)

blomqvist_beta(object)
```

## Arguments

  - object:
    
    a
    [bicop\_dist](https://vinecopulib.github.io/rvinecopulib/reference/bicop_dist.md)
    or fitted
    [bicop](https://vinecopulib.github.io/rvinecopulib/reference/bicop.md)
    object.

## Value

`tail_dep()` returns a 2 by 2 matrix whose rows refer to the lower and
upper tail of the first variable and whose columns refer to the lower
and upper tail of the second variable. `blomqvist_beta()` returns a
numeric scalar.

## Examples

``` r
cop <- bicop_dist("clayton", 0, 2)
tail_dep(cop)
#>          variable2
#> variable1     lower upper
#>     lower 0.7071068     0
#>     upper 0.0000000     0
blomqvist_beta(cop)
#> [1] 0.5118579
```
