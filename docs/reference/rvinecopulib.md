# High Performance Algorithms for Vine Copula Modeling

rvinecopulib provides high-performance tools for constructing, fitting,
selecting, simulating, and visualizing bivariate and vine copula models.
It also fits complete multivariate distributions by combining marginal
models with a vine copula. Continuous, integer-valued discrete, mixed,
and zero-inflated variables are supported.

## Details

The package is the R interface to the header-only 'vinecopulib' C++
library, which is bundled with rvinecopulib. Users do not need to
install the C++ library separately. 'vinecopulib' is licensed under the
MIT License and rvinecopulib under the GNU GPL version 3.

## Main capabilities

rvinecopulib provides:

  - parametric, nonparametric, rotated, and extreme-value pair copulas;

  - automatic pair-copula family, vine structure, truncation, and
    threshold selection;

  - complete multivariate distributions with nonparametric, parametric,
    or custom margins;

  - continuous, integer-valued discrete, mixed, and zero-inflated data;

  - observation weights and parallel fitting;

  - conditional simulation and Rosenblatt transforms;

  - likelihood scores, Hessians, and observation-specific parameters;

  - custom marginal families and custom tree-selection criteria.

## Modeling interfaces

The framework is exposed at three levels. `vine()` models observations
on their original scale by fitting the marginal distributions and
dependence model together. `vinecop()` models uniform
pseudo-observations when only the dependence model is required.
`bicop()` provides the corresponding two-variable workflow.

The constructors `vine_dist()`, `vinecop_dist()`, and `bicop_dist()`
create models from explicitly specified components.

## Learning more

Start with `vignette("getting-started")`. Dedicated articles cover
bivariate copulas, vine structures and selection, marginal models,
discrete data, conditional simulation, and likelihood derivatives. The
[package website](https://vinecopulib.github.io/rvinecopulib/) contains
the rendered articles and complete function reference.

## See also

Useful links:

  - <https://vinecopulib.github.io/rvinecopulib/>

  - Report bugs at <https://github.com/vinecopulib/rvinecopulib/issues>

## Author

Thomas Nagler, Thibault Vatter

## Examples

``` r
## complete distribution on the original data scale
x <- data.frame(a = rnorm(50), b = rexp(50), n = rpois(50, 2))
fit <- vine(
  x,
  var_types = c("c", "c", "d"),
  copula_controls = list(family_set = "onepar")
)
rvine(3, fit)
#>               a         b n
#> [1,] -0.5230271 0.2626702 2
#> [2,] -0.9188961 0.2766614 3
#> [3,] -1.4100087 0.2144410 2

## dependence model on the copula scale
u <- pseudo_obs(as.matrix(USArrests))
copula_fit <- vinecop(u, family_set = "onepar")
summary(copula_fit)
#> # A data.frame: 6 x 11 
#>  tree edge conditioned conditioning var_types   family rotation parameters df
#>     1    1        3, 4                    c,c   gumbel      180        1.5  1
#>     1    2        4, 2                    c,c   gumbel      180        2.1  1
#>     1    3        2, 1                    c,c   gumbel      180        2.5  1
#>     2    1        3, 2            4       c,c gaussian        0     -0.034  1
#>     2    2        4, 1            2       c,c      joe      180        1.4  1
#>     3    1        3, 1         2, 4       c,c  clayton      270       0.43  1
#>     tau loglik
#>   0.318  5.944
#>   0.525 17.882
#>   0.594 24.983
#>  -0.022  0.027
#>   0.192  4.465
#>  -0.177  2.318
#> 
#> ---
#> 1 <-> Murder,   2 <-> Assault,   3 <-> UrbanPop,   4 <-> Rape 

## bivariate model
uv <- rbicop(100, "gaussian", 0, 0.5)
bicop(uv, family_set = "par")
#> Bivariate copula fit ('bicop'): family = gaussian, rotation = 0, parameters = 0.53, var_types = c,c
```
