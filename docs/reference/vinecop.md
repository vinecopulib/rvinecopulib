# Fitting vine copula models

Automated fitting and model selection for vine copula models with
continuous or discrete data. Selection of the structure is performed
using the algorithm of Dissmann et al. (2013).

## Usage

``` r
vinecop(
  data,
  var_types = rep("c", NCOL(data)),
  family_set = "all",
  structure = NA,
  par_method = "mle",
  nonpar_method = "constant",
  mult = 1,
  selcrit = "aic",
  weights = numeric(),
  psi0 = 0.9,
  presel = TRUE,
  allow_rotations = TRUE,
  trunc_lvl = Inf,
  tree_crit = "tau",
  threshold = 0,
  keep_data = FALSE,
  vinecop_object = NULL,
  show_trace = FALSE,
  cores = 1,
  tree_algorithm = "mst_prim",
  conditioning_set = NULL
)
```

## Arguments

  - data:
    
    a matrix or data.frame with at least two columns, containing the
    (pseudo-)observations for the two variables (copula data should have
    approximately uniform margins). More columns are required for
    discrete models, see *Details*.

  - var\_types:
    
    variable types, a length d vector; e.g., `c("c", "c")` for two
    continuous variables, or `c("c", "d")` for first variable continuous
    and second discrete.

  - family\_set:
    
    a character vector of families; see `bicop()` for additional
    options.

  - structure:
    
    an `rvine_structure` object, namely a compressed representation of
    the vine structure, or an object that can be coerced into one (see
    `rvine_structure()` and `as_rvine_structure()`). The dimension must
    be `length(pair_copulas[[1]]) + 1`; `structure = NA` performs
    automatic selection based on Dissman's algorithm. See *Details* for
    partial selection of the structure.

  - par\_method:
    
    the estimation method for parametric models, either `"mle"` for
    maximum likelihood or `"itau"` for inversion of Kendall's tau (only
    available for one-parameter families and `"t"`).

  - nonpar\_method:
    
    the estimation method for nonparametric models, either `"constant"`
    for the standard transformation estimator, or
    `"linear"`/`"quadratic"` for the local-likelihood approximations of
    order one/two.

  - mult:
    
    multiplier for the smoothing parameters of nonparametric families.
    Values larger than 1 make the estimate more smooth, values less than
    1 less smooth.

  - selcrit:
    
    criterion for family selection, either `"loglik"`, `"aic"`, `"bic"`,
    `"mbic"`. For `vinecop()` there is the additional option `"mbicv"`.

  - weights:
    
    optional vector of weights for each observation.

  - psi0:
    
    prior probability of a non-independence copula (only used for
    `selcrit = "mbic"` and `selcrit = "mbicv"`).

  - presel:
    
    whether the family set should be thinned out according to symmetry
    characteristics of the data.

  - allow\_rotations:
    
    whether to allow rotations of the copula.

  - trunc\_lvl:
    
    the truncation level of the vine copula; `Inf` means no truncation,
    `NA` indicates that the truncation level should be selected
    automatically by `mBICV()`.

  - tree\_crit:
    
    the criterion for tree selection, one of `"tau"`, `"rho"`,
    `"hoeffd"`, `"mcor"`, or `"joe"` for Kendall's \\(\\tau\\),
    Spearman's \\(\\rho\\), Hoeffding's \\(D\\), maximum correlation, or
    logarithm of the partial correlation, respectively. Alternatively, a
    function with arguments `data` and `weights` may be supplied. `data`
    is a two-column matrix of pair-copula pseudo-observations. `weights`
    contains the corresponding observation weights, standardized by the
    backend to sum to the original number of observations, or
    `numeric(0)` when no weights were supplied. Rows containing missing
    values or zero weights are removed before the function is called.
    The function must return one finite numeric value; its absolute
    value is used as the edge strength. Custom functions always run
    serially on the calling thread; pair-copula fitting may still use
    `cores` threads.

  - threshold:
    
    for thresholded vine copulas; `NA` indicates that the threshold
    should be selected automatically by `mBICV()`.

  - keep\_data:
    
    whether the data should be stored (necessary for using `fitted()`).

  - vinecop\_object:
    
    a `vinecop` object to be updated; if provided, only the parameters
    are fit; structure and families are kept the same.

  - show\_trace:
    
    logical; whether a trace of the fitting progress should be printed.

  - cores:
    
    number of cores to use; if more than 1, estimation of pair copulas
    within a tree is done in parallel.

  - tree\_algorithm:
    
    The algorithm for building the spanning tree (`"mst_prim"`,
    `"mst_kruskal"`, `"random_weighted"`, or `"random_unweighted"`)
    during the tree-wise structure selection. `"mst_prim"` and
    `"mst_kruskal"` use Prim's and Kruskal's algorithms respectively to
    select the maximum spanning tree, maximizing the sum of the edge
    weights (i.e., `tree_crit`). `"random_weighted"` and
    `"random_unweighted"` use Wilson's algorithm to generate a random
    spanning tree, either with probability proportional to the product
    of the edge weights (weighted) or uniformly (unweighted).

  - conditioning\_set:
    
    variable indices or names to place at the end of the vine order. The
    resulting model can be sampled conditionally with `rvinecop()` by
    supplying `u_cond`. Conditioning-aware selection supports fixed or
    automatically selected truncation levels and requires an MST tree
    algorithm (`"mst_prim"` or `"mst_kruskal"`).

## Value

Objects inheriting from `vinecop` and `vinecop_dist` for `vinecop()`. In
addition to the entries provided by `vinecop_dist()`, there are:

  - `threshold`, the (set or estimated) threshold used for thresholding
    the vine.

  - `data` (optionally, if `keep_data = TRUE` was used), the dataset
    that was passed to `vinecop()`.

  - `controls`, a `list` with fit controls that was passed to
    `vinecop()`.

  - `nobs`, the number of observations that were used to fit the model.

## Details

### Missing data

If there are missing data (i.e., `NA` entries), incomplete observations
are discarded before fitting a pair-copula. This is done on a
pair-by-pair basis so that the maximal available information is used.

### Discrete variables

The dependence measures used to select trees (default: Kendall's tau)
are corrected for ties (see
[wdm::wdm](https://tnagler.github.io/wdm-r/reference/wdm.html)).

Let `n` be the number of observations and `d` the number of variables.
When at least one variable is discrete, two types of "observations" are
required in `data`: the first `n x d` block contains realizations of
\\(F\_{X\_j}(X\_j)\\). The second `n x d` block contains realizations of
\\(F\_{X\_j}(X\_j^-)\\). The minus indicates a left-sided limit of the
cdf. For, e.g., an integer-valued variable, it holds
\\(F\_{X\_j}(X\_j^-) = F\_{X\_j}(X\_j - 1)\\). For continuous variables
the left limit and the cdf itself coincide. Respective columns can be
omitted in the second block.

### Structure selection

Selection of the structure is performed using the algorithm of Dissmann,
J. F., E. C. Brechmann, C. Czado, and D. Kurowicka (2013). *Selecting
and estimating regular vine copulae and application to financial
returns.* Computational Statistics & Data Analysis, 59 (1), 52-69. The
dependence measure used to select trees (default: Kendall's tau) is
corrected for ties and can be changed using the `tree_crit` argument,
which can be set to `"tau"`, `"rho"` or `"hoeffd"`. Both Prim's
(default: `"mst_prim"`) and Kruskal's (`"mst_kruskal"`) algorithms are
available through `tree_algorithm` to set the maximum spanning tree
selection algorithm. An alternative to the maximum spanning tree
selection is to use random spanning trees, which can be selected using
`tree_algorithm` and come in two flavors, both using Wilson's algorithm
loop erased random walks:

  - "random\_weighted"\` generates a random spanning tree with
    probability proportional to the product of the weights (i.e., the
    dependence) of the edges in the tree.

  - "random\_unweighted"\` generates a random spanning tree uniformly
    over all spanning trees satisfying the proximity condition.

### Partial structure selection

It is possible to fix the vine structure only in the first trees and
select the remaining ones automatically. To specify only the first `k`
trees, supply a `k`-truncated `rvine_structure()` or `rvine_matrix()`.
All trees up to `trunc_lvl` will then be selected automatically.

## References

Dissmann, J. F., E. C. Brechmann, C. Czado, and D. Kurowicka (2013).
*Selecting and estimating regular vine copulae and application to
financial returns.* Computational Statistics & Data Analysis, 59 (1),
52-69.

## See also

`vinecop()`, `dvinecop()`, `pvinecop()`, `rvinecop()`, `plot.vinecop()`,
`contour.vinecop()`

## Examples

``` r
## simulate dummy data
x <- rnorm(30) * matrix(1, 30, 5) + 0.5 * matrix(rnorm(30 * 5), 30, 5)
u <- pseudo_obs(x)

## fit and select the model structure, family and parameters
fit <- vinecop(u)
summary(fit)
#> # A data.frame: 10 x 11 
#>  tree edge conditioned conditioning var_types  family rotation parameters df
#>     1    1        4, 1                    c,c     bb7        0   2.8, 3.8  2
#>     1    2        5, 3                    c,c clayton        0        3.7  1
#>     1    3        3, 1                    c,c     bb7        0   2.6, 8.9  2
#>     1    4        2, 1                    c,c     bb7        0   3.6, 2.9  2
#>     2    1        4, 3            1       c,c   indep        0             0
#>     2    2        5, 1            3       c,c   indep        0             0
#>     2    3        3, 2            1       c,c clayton      180       0.32  1
#>     3    1        4, 5         3, 1       c,c     joe        0        1.4  1
#>     3    2        5, 2         1, 3       c,c     joe        0        1.5  1
#>     4    1        4, 2      5, 3, 1       c,c     joe      270        1.3  1
#>    tau loglik
#>   0.70   23.4
#>   0.65   21.2
#>   0.78   34.7
#>   0.71   23.1
#>   0.00    0.0
#>   0.00    0.0
#>   0.14    1.3
#>   0.19    2.2
#>   0.21    3.5
#>  -0.15    1.6
plot(fit)

contour(fit)


## select by BIC from one-parameter families
fit <- vinecop(u, family_set = "onepar", selcrit = "bic")
summary(fit)
#> # A data.frame: 10 x 11 
#>  tree edge conditioned conditioning var_types   family rotation parameters df
#>     1    1        4, 1                    c,c gaussian        0       0.89  1
#>     1    2        5, 3                    c,c  clayton        0        3.7  1
#>     1    3        3, 1                    c,c  clayton        0        6.5  1
#>     1    4        2, 1                    c,c   gumbel        0        3.2  1
#>     2    1        4, 3            1       c,c  clayton        0       0.42  1
#>     2    2        5, 1            3       c,c    frank        0        1.1  1
#>     2    3        3, 2            1       c,c  clayton      180       0.33  1
#>     3    1        4, 5         3, 1       c,c      joe        0        1.4  1
#>     3    2        5, 2         1, 3       c,c      joe        0        1.5  1
#>     4    1        4, 2      5, 3, 1       c,c  clayton       90       0.32  1
#>    tau loglik
#>   0.70  21.50
#>   0.65  21.16
#>   0.76  32.36
#>   0.68  20.58
#>   0.17   1.37
#>   0.13   0.55
#>   0.14   1.14
#>   0.19   2.31
#>   0.21   3.70
#>  -0.14   0.96

## 1-truncated, Gaussian D-vine
fit <- vinecop(
  u,
  structure = dvine_structure(1:5),
  family_set = "gauss",
  trunc_lvl = 1
)
plot(fit)

contour(fit)


## Partial structure selection with only first tree specified
structure <- rvine_structure(order = 1:5, list(rep(5, 4)))
structure
#> 5-dimensional R-vine structure ('rvine_structure'), 1-truncated
#> 5 5 5 5 5
#>       4  
#>     3    
#>   2      
#> 1        
fit <- vinecop(u, structure = structure, family_set = "gauss")
plot(fit)


## Model for discrete data
x <- qpois(u, 1)  # transform to Poisson margins
# we require two types of observations (see Details)
u_disc <- cbind(ppois(x, 1), ppois(x - 1, 1))
fit <- vinecop(u_disc, var_types = rep("d", 5))

## Model for mixed data
x <- qpois(u[, 1], 1)  # transform first variable to Poisson margin
# we require two types of observations (see Details)
u_disc <- cbind(ppois(x, 1), u[, 2:5], ppois(x - 1, 1))
fit <- vinecop(u_disc, var_types = c("d", rep("c", 4)))
```
