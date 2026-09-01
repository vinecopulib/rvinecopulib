# Truncate a vine copula model

Extracts a truncated sub-vine based on a truncation level supplied by
user.

## Usage

``` r
truncate_model(object, ...)

# S3 method for class 'rvine_structure'
truncate_model(object, trunc_lvl, ...)

# S3 method for class 'rvine_matrix'
truncate_model(object, trunc_lvl, ...)

# S3 method for class 'vinecop_dist'
truncate_model(object, trunc_lvl, ...)

# S3 method for class 'vine_dist'
truncate_model(object, trunc_lvl, ...)
```

## Arguments

  - object:
    
    a model object.

  - ...:
    
    further arguments passed to specific methods.

  - trunc\_lvl:
    
    tree level after which the vine copula should be truncated.

## Details

While a vine model for a `d` dimensional random vector contains at most
`d-1` nested trees, this function extracts a sub-model based on a given
truncation level.

For instance, `truncate_model(object, 1)` results in a 1-truncated vine
(i.e., a vine with a single tree). Similarly `truncate_model(object, 2)`
results in a 2-truncated vine (i.e., a vine with two trees). Note that
`truncate_model(truncate_model(object, 1), 2)` returns a 1-truncated
vine.

## Examples

``` r
# specify pair-copulas
bicop <- bicop_dist("bb1", 90, c(3, 2))
pcs <- list(
  list(bicop, bicop), # pair-copulas in first tree
  list(bicop) # pair-copulas in second tree
)

# specify R-vine matrix
mat <- matrix(c(1, 2, 3, 1, 2, 0, 1, 0, 0), 3, 3)

# set up vine structure
structure <- as_rvine_structure(mat)

# truncate the model
truncate_model(structure, 1)
#> 3-dimensional R-vine structure ('rvine_structure'), 1-truncated
#> 1 1 1
#>   2  
#> 3    

# set up vine copula model
vc <- vinecop_dist(pcs, mat)

# truncate the model
truncate_model(vc, 1)
#> 3-dimensional vine copula model ('vinecop_dist'), 1-truncated
```
