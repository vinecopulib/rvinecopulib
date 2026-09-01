# Extracts components of `bicop_dist` and `vinecop_dist` objects

Extracts either the structure matrix (for `vinecop_dist` only), or
pair-copulas, their parameters, Kendall's taus, or families (for
`bicop_dist` and `vinecop_dist`).

## Usage

``` r
get_structure(object)

get_pair_copula(object, tree = NA, edge = NA)

get_parameters(object, tree = NA, edge = NA)

get_ktau(object, tree = NA, edge = NA)

get_family(object, tree = NA, edge = NA)

get_all_pair_copulas(object, trees = NA)

get_all_parameters(object, trees = NA)

get_all_ktaus(object, trees = NA)

get_all_families(object, trees = NA)
```

## Arguments

  - object:
    
    a `bicop_dist`, `vinecop_dist` or `vine_dist` object.

  - tree:
    
    tree index (not required if `object` is of class `bicop_dist`).

  - edge:
    
    edge index (not required if `object` is of class `bicop_dist`).

  - trees:
    
    the trees to extract from `object` (`trees = NA` extracts all
    trees).

## Value

The structure matrix, or pair-copulas, their parameters, Kendall's taus,
or families.

## Details

\#' The get\_structure method (for `vinecop_dist` or `vine_dist` objects
only) extracts the structure (see
[rvine\_structure](https://vinecopulib.github.io/rvinecopulib/reference/rvine_structure.md)
for more details).

The get\_matrix method (for `vinecop_dist` or `vine_dist` objects only)
extracts the structure matrix (see
[rvine\_structure](https://vinecopulib.github.io/rvinecopulib/reference/rvine_structure.md)
for more details).

The other `get_xyz` methods for `vinecop_dist` or `vine_dist` objects
return the entries corresponding to the pair-copula indexed by its
`tree` and `edge`. When `object` is of class `bicop_dist`, `tree` and
`edge` are not required.  

  - `get_pair_copula()` = the pair-copula itself (see
    [bicop](https://vinecopulib.github.io/rvinecopulib/reference/bicop.md)).

  - `get_parameters()` = the parameters of the pair-copula (i.e., a
    `numeric` scalar, vector, or matrix).

  - `get_family()` = a character for the family (see
    [bicop](https://vinecopulib.github.io/rvinecopulib/reference/bicop.md)
    for implemented families).

  - `get_ktau()` = a `numeric` scalar with the pair-copula Kendall's
    tau.

The `get_all_xyz` methods (for `vinecop_dist` or `vine_dist` objects
only) return lists of lists, with each element corresponding to a tree
in `trees`, and then elements of the sublists correspond to edges. The
returned lists have two additional attributes:

  - `"d"` = the dimension of the model.

  - `"trees"` = the extracted trees.

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

# set up vine copula model
vc <- vinecop_dist(pcs, mat)

# get the structure
get_structure(vc)
#> 3-dimensional R-vine structure ('rvine_structure')
#> 1 1 1
#> 2 2  
#> 3    
all(get_matrix(vc) == mat)
#> [1] TRUE

# get pair-copulas
get_pair_copula(vc, 1, 1)
#> Bivariate copula ('bicop_dist'): family = bb1, rotation = 90, parameters = 3, 2, var_types = c,c
get_all_pair_copulas(vc)
#> Nested list of lists for the pair-copulas of a 3 dimensional vine with all trees: 
#> - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#> x[[1]] -> a list with the 2 pair-copulas of tree 1. 
#> x[[2]] -> a list with the 1 pair-copulas of tree 2. 
all.equal(get_all_pair_copulas(vc), pcs,
          check.attributes = FALSE, check.environment = FALSE)
#> [1] TRUE
```
