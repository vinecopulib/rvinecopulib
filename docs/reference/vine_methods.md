# Vine based distributions

Density, distribution function and random generation for the vine based
distribution.

## Usage

``` r
dvine(x, vine, cores = 1)

pvine(x, vine, n_mc = 10^4, cores = 1)

rvine(n, vine, qrng = FALSE, cores = 1, x_cond = NULL, conditioning_set = NULL)
```

## Arguments

  - x:
    
    evaluation points, either a length d vector or a d-column matrix,
    where d is the number of variables in the vine.

  - vine:
    
    an object of class `"vine_dist"`.

  - cores:
    
    number of cores to use; if larger than one, computations are done in
    parallel on `cores` batches .

  - n\_mc:
    
    number of samples used for quasi Monte Carlo integration.

  - n:
    
    number of observations.

  - qrng:
    
    if `TRUE`, generates quasi-random numbers using the multivariate
    Generalized Halton sequence up to dimension 300 and the Generalized
    Sobol sequence in higher dimensions (default `qrng = FALSE`).

  - x\_cond:
    
    optional conditioning values for `rvine()` on the original data
    scale. A vector or one-row object is repeated `n` times;
    alternatively, supply an `n`-row matrix or data frame for
    observation-specific conditioning values. If `NULL`, `rvine()`
    performs unconditional simulation.

  - conditioning\_set:
    
    variable indices or names corresponding to the columns of `x_cond`.
    When `NULL`, the columns correspond to the last variables of the
    current copula order. Discrete left limits are computed internally
    from the fitted margins.

## Value

`dvine()` gives the density, `pvine()` gives the distribution function,
and `rvine()` generates unconditional or conditional random deviates.

The length of the result is determined by `n` for `rvine()`, and the
number of rows in `u` for the other functions.

The `vine` object is recycled to the length of the result.

## Details

See [vine](https://vinecopulib.github.io/rvinecopulib/reference/vine.md)
for the estimation and construction of vine models. Here, the density,
distribution function and random generation for the vine distributions
are standard.

The functions are based on `dvinecop()`, `pvinecop()` and `rvinecop()`
for
[vinecop](https://vinecopulib.github.io/rvinecopulib/reference/vinecop.md)
objects. Margins are evaluated through `dmargin()`, `pmargin()`, and
`qmargin()`. Methods are provided for margins fitted by `vine()` and for
the fixed
[stats::Distributions](https://rdrr.io/r/stats/Distributions.html)
specifications accepted by `vine_dist()`.

## Examples

``` r
# specify pair-copulas
bicop <- bicop_dist("bb1", 90, c(3, 2))
pcs <- list(
  list(bicop, bicop), # pair-copulas in first tree
  list(bicop) # pair-copulas in second tree
)

# set up vine copula model
mat <- rvine_matrix_sim(3)
vc <- vine_dist(list(stats_margin("norm")), pcs, mat)

# simulate from the model
x <- rvine(200, vc)
pairs(x)


# evaluate the density and cdf
dvine(x[1, ], vc)
#> [1] 0.3210905
pvine(x[1, ], vc)
#> [1] 1e-04
```
