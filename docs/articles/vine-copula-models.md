# Vine copula models and structures

A vine copula decomposes a multivariate dependence model into a sequence
of bivariate copulas. The structure determines which conditioned and
conditioning variables belong to each edge; one `bicop_dist` object
supplies the model for that edge.

## Read an R-vine structure

Consider the triangular array

    4 4 4 4
    3 3 3
    2 2
    1

It represents these pair copulas:

| Tree | Edge | Conditioned pair | Conditioning set |
| ---: | ---: | ---------------- | ---------------- |
|    1 |    1 | 1, 4             | —                |
|    1 |    2 | 2, 4             | —                |
|    1 |    3 | 3, 4             | —                |
|    2 |    1 | 1, 3             | 4                |
|    2 |    2 | 2, 3             | 4                |
|    3 |    1 | 1, 2             | 3, 4             |

In R, the structure can be constructed from its order and off-diagonal
rows:

``` r
structure <- rvine_structure(
  order = 1:4,
  struct_array = list(c(4, 4, 4), c(3, 3), 2)
)
structure
#> 4-dimensional R-vine structure ('rvine_structure')
#> 4 4 4 4
#> 3 3 3  
#> 2 2    
#> 1
as_rvine_matrix(structure)
#> 4-dimensional R-vine matrix ('rvine_matrix')
#> 4 4 4 4
#> 3 3 3  
#> 2 2    
#> 1
```

The square R-vine matrix contains the same information with zeros
filling the unused triangle. `as_rvine_structure()` converts it back to
the compressed representation.

`cvine_structure()` and `dvine_structure()` construct the two common
order-determined special cases. Random valid structures are available
through `rvine_structure_sim()`.

``` r
cvine_structure(c(4, 1, 3, 2))
#> 4-dimensional R-vine structure ('rvine_structure')
#> 2 2 2 2
#> 3 3 3  
#> 1 1    
#> 4
dvine_structure(c(4, 1, 3, 2))
#> 4-dimensional R-vine structure ('rvine_structure')
#> 1 3 2 2
#> 3 2 3  
#> 2 1    
#> 4
rvine_structure_sim(4)
#> 4-dimensional R-vine structure ('rvine_structure')
#> 3 3 3 3
#> 1 1 1  
#> 4 4    
#> 2
```

See `?rvine_structure` for the complete validity and proximity
conditions.

## Construct a vine copula from components

For a three-dimensional model, the first tree has two pair copulas and
the second has one. The nested list follows
`pair_copulas[[tree]][[edge]]`.

``` r
pair_copulas <- list(
  list(
    bicop_dist("gaussian", parameters = 0.6),
    bicop_dist("clayton", parameters = 1.5)
  ),
  list(bicop_dist("frank", parameters = 2))
)
model <- vinecop_dist(
  pair_copulas,
  structure = dvine_structure(1:3)
)
model
#> 3-dimensional vine copula model ('vinecop_dist')
summary(model)
#> # A data.frame: 3 x 10 
#>  tree edge conditioned conditioning var_types   family rotation parameters df
#>     1    1        1, 2                    c,c gaussian        0        0.6  1
#>     1    2        2, 3                    c,c  clayton        0        1.5  1
#>     2    1        1, 3            2       c,c    frank        0          2  1
#>   tau
#>  0.41
#>  0.43
#>  0.21
```

The model can immediately be evaluated and simulated:

``` r
u <- rvinecop(5, model)
dvinecop(u, model)
#> [1] 0.3040095 0.3516670 1.8689282 1.3118066 1.7333509
pvinecop(u[1:2, ], model, n_mc = 1000)
#> [1] 0.143 0.204
```

## Fit and select a model

With no structure supplied, `vinecop()` selects trees sequentially using
the Dissmann algorithm and fits a pair copula on every selected edge.

``` r
n <- 150
z <- rnorm(n)
x <- cbind(
  z + rnorm(n, sd = 0.7),
  0.6 * z + rnorm(n),
  -0.4 * z + rnorm(n),
  rnorm(n)
)
u <- pseudo_obs(x)

fit <- vinecop(u, family_set = "onepar", selcrit = "bic")
fit
#> 4-dimensional vine copula fit ('vinecop')
#> nobs = 150   logLik = 17.21   npars = 6   AIC = -22.42   BIC = -4.36
summary(fit)
#> # A data.frame: 6 x 11 
#>  tree edge conditioned conditioning var_types   family rotation parameters df
#>     1    1        4, 2                    c,c      joe        0        1.1  1
#>     1    2        3, 1                    c,c    frank        0       -1.3  1
#>     1    3        2, 1                    c,c gaussian        0       0.33  1
#>     2    1        4, 1            2       c,c      joe       90        1.1  1
#>     2    2        3, 2            1       c,c  clayton      270       0.15  1
#>     3    1        4, 3         1, 2       c,c      joe       90        1.2  1
#>     tau loglik
#>   0.060    1.4
#>  -0.143    3.3
#>   0.217    8.1
#>  -0.043    1.3
#>  -0.068    1.2
#>  -0.087    2.0
```

The most important pair-copula controls (`family_set`, `par_method`,
`nonpar_method`, `mult`, `selcrit`, `weights`, `presel`, and
`allow_rotations`) have the same meaning as in `bicop()`. The [bivariate
copula
guide](https://vinecopulib.github.io/rvinecopulib/articles/bivariate-copulas.md)
describes the available families, rotations, and estimation methods.

## Inspect the selected trees

Getters work on both manually constructed and fitted models. Tree and
edge indices are one-based in the R interface.

``` r
get_structure(fit)
#> 4-dimensional R-vine structure ('rvine_structure')
#> 2 1 1 1
#> 1 2 2  
#> 3 3    
#> 4
get_matrix(fit)
#> 4-dimensional R-vine matrix ('rvine_matrix')
#> 2 1 1 1
#> 1 2 2  
#> 3 3    
#> 4
get_pair_copula(fit, tree = 1, edge = 1)
#> Bivariate copula fit ('bicop'): family = joe, rotation = 0, parameters = 1.11, var_types = c,c
get_family(fit, tree = 1, edge = 1)
#> [1] "joe"
get_parameters(fit, tree = 1, edge = 1)
#>          [,1]
#> [1,] 1.111375
get_ktau(fit, tree = 1, edge = 1)
#> [1] 0.06016693
get_all_families(fit)
#> Nested list of lists for the copula families of a 4 dimensional vine with all trees: 
#> - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#> x[[1]] -> a list with the 3 copula families of tree 1. 
#> x[[2]] -> a list with the 2 copula families of tree 2. 
#> x[[3]] -> a list with the 1 copula families of tree 3.
```

`summary(fit)` returns the edge table invisibly, which is convenient for
programmatic inspection.

## Control structure selection

The edge weight used to build each tree is selected by `tree_crit`:

  - `"tau"` for absolute Kendall’s tau;
  - `"rho"` for absolute Spearman’s rho;
  - `"hoeffd"` for Hoeffding’s D;
  - `"mcor"` for maximum correlation;
  - `"joe"` for a Gaussian-copula mutual-information criterion;
  - a custom function of `data` and `weights`.

`tree_algorithm = "mst_prim"` and `"mst_kruskal"` select a maximum
spanning tree. `"random_weighted"` and `"random_unweighted"` draw random
spanning trees, which can be useful in ensemble or exploratory
workflows.

``` r
custom_strength <- function(data, weights) {
  if (length(weights)) {
    abs(weighted.mean(data[, 1] * data[, 2], weights))
  } else {
    abs(mean(data[, 1] * data[, 2]))
  }
}

custom_fit <- vinecop(u, tree_crit = custom_strength)
random_fit <- vinecop(u, tree_algorithm = "random_weighted")
```

Custom criteria run on the calling thread. Pair-copula fitting may still
use multiple `cores`.

## Fix all or part of a structure

Supply a full structure to keep it fixed while selecting pair-copula
families. A truncated structure fixes its existing trees and asks
`vinecop()` to select the remaining trees up to `trunc_lvl`.

``` r
fixed_structure_fit <- vinecop(
  u,
  structure = dvine_structure(1:4),
  family_set = c("gaussian", "clayton")
)

first_tree <- rvine_structure(order = 1:4, struct_array = list(rep(4, 3)))
partial_fit <- vinecop(
  u,
  structure = first_tree,
  family_set = "onepar"
)
```

## Truncation and thresholding

`trunc_lvl` limits the number of fitted trees. Pair copulas above that
level are independence copulas. `threshold` replaces edges whose
selection strength is below the threshold by independence copulas.

``` r
truncated <- vinecop(u, family_set = "onepar", trunc_lvl = 1)
thresholded <- vinecop(u, family_set = "onepar", threshold = 0.1)
truncate_model(fit, trunc_lvl = 1)
#> 4-dimensional vine copula fit ('vinecop'), 1-truncated
#> nobs = 150   logLik = 12.74   npars = 3   AIC = -19.49   BIC = -10.46
```

Set `trunc_lvl = NA` or `threshold = NA` to select the value
automatically with mBICV. The fitted threshold and truncation level are
visible in the model summary.

## Refit an existing model

Pass a fitted model through `vinecop_object` to update parameters while
keeping its structure and families fixed. This is useful for rolling
windows or bootstrap samples.

``` r
refitted <- vinecop(u, vinecop_object = fit)
stopifnot(identical(get_all_families(refitted), get_all_families(fit)))
```

For a manually constructed distribution, use `vinecop()` to select from
data or construct a new `vinecop_dist()` after changing its components.

## Related documentation

  - [Bivariate copula
    models](https://vinecopulib.github.io/rvinecopulib/articles/bivariate-copulas.md)
    provides the details of the pair-copula building blocks.
  - [Discrete, mixed, and zero-inflated
    data](https://vinecopulib.github.io/rvinecopulib/articles/discrete-data.md)
    explains the copula-scale representation required when some
    variables have atoms.
  - [Conditional simulation and Rosenblatt
    transforms](https://vinecopulib.github.io/rvinecopulib/articles/conditional-simulation.md)
    describes conditioning-aware structures and sampling orders.
  - [Scores, Hessians, and varying
    parameters](https://vinecopulib.github.io/rvinecopulib/articles/likelihood-inference.md)
    explains edge parameter ordering and stepwise versus full
    derivatives.
  - The [`vinecop()`
    reference](https://vinecopulib.github.io/rvinecopulib/reference/vinecop.md)
    lists all fitting and selection arguments; the [structure
    reference](https://vinecopulib.github.io/rvinecopulib/reference/rvine_structure.md)
    gives the formal matrix validity conditions.
