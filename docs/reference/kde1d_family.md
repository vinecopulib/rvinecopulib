# Define a kde1d margin family

The arguments configure `kde1d::kde1d()` and travel with the family
object, rather than being interpreted by `vine()`.

## Usage

``` r
kde1d_family(xmin = NaN, xmax = NaN, mult = 1, bw = NA, deg = 2)
```

## Arguments

  - xmin, xmax, mult, bw, deg:
    
    arguments passed to `kde1d::kde1d()`.

## Value

A margin-family specification supporting all rvinecopulib variable
types.
