# Convert list to bicop object

Convert list to bicop object

## Usage

``` r
as.bicop(object, check = TRUE)
```

## Arguments

  - object:
    
    a list containing entries for `"family"`, `"rotation"`,
    `"parameters"`, and `"npars"`.

  - check:
    
    whether to check for validity of the family/parameter specification.

## Value

A bicop object corresponding to the specification in `object`.

## Examples

``` r
as.bicop(list(family = "gumbel", rotation = 90, parameters = 2, npars = 1))
#> Bivariate copula fit ('bicop'): family = gumbel, rotation = 90, parameters = 2, var_types = c,c
```
