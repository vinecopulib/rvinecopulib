# Conversion between Kendall's tau and parameters

Conversion between Kendall's tau and parameters

## Usage

``` r
par_to_ktau(family, rotation, parameters)

ktau_to_par(family, tau)
```

## Arguments

  - family:
    
    a copula family (see `bicop_dist()`) or a
    [bicop\_dist](https://vinecopulib.github.io/rvinecopulib/reference/bicop_dist.md)
    object.

  - rotation:
    
    the rotation of the copula, one of `0`, `90`, `180`, `270`.

  - parameters:
    
    vector or matrix of copula parameters, not used when `family` is a
    `bicop_dist` object.

  - tau:
    
    Kendall's \\(\\tau\\).

## Examples

``` r
# the following are equivalent
par_to_ktau(bicop_dist("clayton", 0, 3))
#> [1] 0.6
par_to_ktau("clayton", 0, 3)
#> [1] 0.6

ktau_to_par("clayton", 0.5)
#>      [,1]
#> [1,]    2
ktau_to_par(bicop_dist("clayton", 0, 3), 0.5)
#>      [,1]
#> [1,]    2
```
