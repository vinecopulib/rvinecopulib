# Plotting R-vine structures

Plot one or all trees of an R-vine structure.

## Usage

``` r
# S3 method for class 'rvine_structure'
plot(x, ...)

# S3 method for class 'rvine_matrix'
plot(x, ...)
```

## Arguments

  - x:
    
    an `rvine_structure` or `rvine_matrix` object.

  - ...:
    
    passed to `plot.vinecop_dist()`.

## Examples

``` r
plot(cvine_structure(1:5))

plot(rvine_structure_sim(5))

mat <- rbind(c(1, 1, 1), c(2, 2, 0), c(3, 0, 0))
plot(rvine_matrix(mat))

plot(rvine_matrix_sim(5))
```
