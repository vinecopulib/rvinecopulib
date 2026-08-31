# Coerce various kind of objects to R-vine structures and matrices

`as_rvine_structure` and `as_rvine_matrix` are new S3 generics allowing
to coerce objects into R-vine structures and matrices (see
`rvine_structure()` and `rvine_matrix()`).

## Usage

``` r
as_rvine_structure(x, ...)

as_rvine_matrix(x, ...)

# S3 method for class 'rvine_structure'
as_rvine_structure(x, ..., validate = FALSE)

# S3 method for class 'rvine_structure'
as_rvine_matrix(x, ..., validate = FALSE)

# S3 method for class 'list'
as_rvine_structure(x, ..., is_natural_order = FALSE)

# S3 method for class 'list'
as_rvine_matrix(x, ..., is_natural_order = FALSE)

# S3 method for class 'rvine_matrix'
as_rvine_structure(x, ..., validate = FALSE)

# S3 method for class 'rvine_matrix'
as_rvine_matrix(x, ..., validate = FALSE)

# S3 method for class 'matrix'
as_rvine_structure(x, ..., validate = TRUE)

# S3 method for class 'matrix'
as_rvine_matrix(x, ..., validate = TRUE)
```

## Arguments

  - x:
    
    An object of class `rvine_structure`, `rvine_matrix`, `matrix` or
    `list` that can be coerced into an R-vine structure or R-vine matrix
    (see *Details*).

  - ...:
    
    Other arguments passed on to individual methods.

  - validate:
    
    When \`TRUE“, verifies that the input is a valid rvine-structure
    (see *Details*). You may want to suppress this when you know that
    you already have a valid structure and you want to save some time,
    or to explicitly enable it if you have a structure that you want to
    re-check.

  - is\_natural\_order:
    
    A flag indicating whether the `struct_array` element of `x` is
    assumed to be provided in natural order already (a structure is in
    natural order if the anti-diagonal is 1, .., d from bottom left to
    top right).

## Value

Either an object of class `rvine_structure` or of class `rvine_matrix`
(see `rvine_structure()` or `rvine_matrix()`).

## Details

The coercion to `rvine_structure` and `rvine_matrix` can be applied to
different kind of objects Currently, `rvine_structure`, `rvine_matrix`,
`matrix` and `list` are supported.

For `as_rvine_structure`:

  - `rvine_structure` : the main use case is to re-check an object via
    `validate = TRUE`.

  - `rvine_matrix` and `matrix` : allow to coerce matrices into R-vine
    structures (see `rvine_structure()` for more details). The main
    difference between `rvine_matrix` and `matrix` is the nature of the
    validity checks.

  - `list` : must contain named elements `order` and `struct_array` to
    be coerced into an R-vine structure (see `rvine_structure()` for
    more details).

For `as_rvine_matrix`:

  - `rvine_structure` : allow to coerce an `rvine_structure` into an
    R-vine matrix (useful e.g. for printing).

  - `rvine_matrix`: similar to `as_rvine_structure` for
    `rvine_structure`, the main use case is to re-check an object via
    `validate = TRUE`.

  - `matrix` : allow to coerce matrices into R-vine matrices (mainly by
    checking that the matrix defines a valid R-vine, see
    `rvine_matrix()` for more details).

  - `list` : must contain named elements `order` and `struct_array` to
    be coerced into an R-vine matrix (see `rvine_structure()` for more
    details).

## See also

rvine\_structure rvine\_matrix

## Examples

``` r
# R-vine structures can be constructed from the order vector and struct_array
rvine_structure(order = 1:4, struct_array = list(
  c(4, 4, 4),
  c(3, 3),
  2
))
#> 4-dimensional R-vine structure ('rvine_structure')
#> 4 4 4 4
#> 3 3 3  
#> 2 2    
#> 1      

# ... or a similar list can be coerced into an R-vine structure
as_rvine_structure(list(order = 1:4, struct_array = list(
  c(4, 4, 4),
  c(3, 3),
  2
)))
#> 4-dimensional R-vine structure ('rvine_structure')
#> 4 4 4 4
#> 3 3 3  
#> 2 2    
#> 1      

# similarly, standard matrices can be coerced into R-vine structures
mat <- matrix(c(4, 3, 2, 1, 4, 3, 2, 0, 4, 3, 0, 0, 4, 0, 0, 0), 4, 4)
as_rvine_structure(mat)
#> 4-dimensional R-vine structure ('rvine_structure')
#> 4 4 4 4
#> 3 3 3  
#> 2 2    
#> 1      

# or truncate and construct the structure
mat[3, 1] <- 0
as_rvine_structure(mat)
#> 4-dimensional R-vine structure ('rvine_structure'), 2-truncated
#> 4 4 4 4
#> 3 3 3  
#>   2    
#> 1      

# throws an error
mat[3, 1] <- 5
try(as_rvine_structure(mat))
#> Error : not a valid R-vine array: the upper left triangle can only contain numbers between 1 and d (number of variables).
```
