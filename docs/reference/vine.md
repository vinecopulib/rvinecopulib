# Vine copula models

Automated fitting or creation of custom vine copula models

## Usage

``` r
vine(
  data,
  margins_controls = list(),
  copula_controls = list(family_set = "all", structure = NA, par_method = "mle",
    nonpar_method = "constant", mult = 1, selcrit = "aic", psi0 = 0.9, presel = TRUE,
    allow_rotations = TRUE, trunc_lvl = Inf, tree_crit = "tau", threshold = 0, keep_data
    = FALSE, show_trace = FALSE, cores = 1, tree_algorithm = "mst_prim", conditioning_set
    = NULL),
  weights = numeric(),
  keep_data = FALSE,
  cores = 1,
  var_types = NULL
)

vine_dist(margins, pair_copulas, structure)
```

## Arguments

  - data:
    
    a matrix or data.frame. Discrete variables have to be declared as
    `ordered()`.

  - margins\_controls:
    
    a list with arguments to be passed to marginal fitting. The
    supported entries are
    
      - `family_set` candidate families used for every variable, or a
        list with one entry per variable. The default `"kde1d"`
        preserves nonparametric margin fitting. Other character names
        refer to `univariateML` families; configured built-in and custom
        candidates are created with `kde1d_family()`,
        `univariateML_family()`, and `margin_family()`.
    
      - `selcrit` selection criterion, one of `"loglik"`, `"aic"`, or
        `"bic"`.
    
      - `cores` number of cores used for marginal fitting. `NULL`
        inherits the top-level `cores` argument. Use one core for
        stochastic custom fitters when results must be invariant to the
        number of cores.

  - copula\_controls:
    
    a list with arguments to be passed to `vinecop()`.

  - weights:
    
    optional numeric vector of weights for each observation, used for
    both marginal and copula estimation. Every margin family receives
    these weights and is responsible for using or explicitly rejecting
    them.

  - keep\_data:
    
    whether the original data should be stored; if you want to store the
    pseudo-observations used for fitting the copula, use the
    `copula_controls` argument.

  - cores:
    
    the number of cores to use for parallel computations. Unless
    overridden by `margins_controls$cores`, margins are fitted in forked
    processes when `cores > 1` on non-Windows systems and serially on
    Windows. Stochastic custom margin fits may depend on the number of
    forked processes.

  - var\_types:
    
    optional variable types, one for each variable after factor
    expansion: `"c"` for continuous, `"d"` for integer-valued discrete,
    or `"zi"` for continuous with an atom at zero. Types are inferred
    from `ordered()` and `zero_inflated()` columns when omitted.

  - margins:
    
    a list containing one marginal distribution per variable. Each
    margin can be a `margin_dist()` object, another fitted object
    implementing the [margin
    protocol](https://vinecopulib.github.io/rvinecopulib/reference/margin_protocol.md),
    a `kde1d` object, or a `stats_margin()` object. Legacy `list(distr =
    ...)` specifications are converted through `stats_margin()`. If
    `margins` has length one, it is recycled for every component.

  - pair\_copulas:
    
    A nested list of 'bicop\_dist' objects, where
    `pair_copulas[[t]][[e]]` corresponds to the pair-copula at edge `e`
    in tree `t`.

  - structure:
    
    an `rvine_structure` object, namely a compressed representation of
    the vine structure, or an object that can be coerced into one (see
    `rvine_structure()` and `as_rvine_structure()`). The dimension must
    be `length(pair_copulas[[1]]) + 1`.

## Value

Objects inheriting from `vine_dist` for `vine_dist()`, and `vine` and
`vine_dist` for `vine()`.

Objects from the `vine_dist` class are lists containing:

  - `margins`, a list of marginals (see below).

  - `copula`, an object of the class `vinecop_dist`, see
    `vinecop_dist()`.

For objects from the `vine` class, `copula` is also an object of the
class `vine`, see `vinecop()`. Additionally, objects from the `vine`
class contain:

  - `margins_controls`, a `list` with the controls used to fit and
    select the marginal families.

  - `copula_controls`, a `list` with the set of fit controls that was
    passed to `vinecop()` when estimating the copula.

  - `data` (optionally, if `keep_data = TRUE` was used), the dataset
    that was passed to `vine()`.

  - `nobs`, an `integer` containing the number of observations that was
    used to fit the model.

Concerning `margins`:

  - For objects created with `vine_dist()`, it simply corresponds to the
    `margins` argument.

  - For objects created with `vine()`, it is a list of fitted objects
    implementing the [margin
    protocol](https://vinecopulib.github.io/rvinecopulib/reference/margin_protocol.md).

## Details

`vine_dist()` creates a vine copula by specifying the margins, a nested
list of `bicop_dist` objects and a quadratic structure matrix.

`vine()` provides automated fitting for vine copula models.
Margin-specific fitting options are properties of the family
specification, not controls interpreted by `vine()`. For example, use
`kde1d_family(xmin = 0, mult = 2)` in `family_set`. `copula_controls` is
a list with the same parameters as `vinecop()` (except for `data`).

A character `margins_controls$family_set` is used for every variable.
Its entries can be `"kde1d"` or names of models available in
univariateML. A list supplies one candidate set per variable, for
example `list(c("norm", "cauchy"), c("pois", "geom"))`. Each entry may
itself be a character vector, a margin-family object, or a list
combining both. An unnamed list is matched by position; a named list
must contain every expanded variable name exactly once and is reordered
to match the data.

The aliases `"par"` and `"parametric"` expand to all `univariateML`
families, `"nonpar"` and `"nonparametric"` select `"kde1d"`, and `"all"`
combines both. Parametric names and aliases require the suggested
`univariateML` package; the default `"kde1d"`, `kde1d_family()`, and
custom `margin_family()` objects do not. Identical candidate
specifications are fitted only once.

Candidates incompatible with a variable's type are removed before
fitting. rvinecopulib fits every remaining candidate and selects the
first minimum of negative log-likelihood, AIC, or BIC according to
`selcrit`. Competing candidates therefore need a finite
`margin_info()$loglik` value, and AIC or BIC additionally requires a
finite `margin_info()$npars`. A sole candidate without a parameter count
is retained with a warning, but model AIC and BIC are then unavailable.
An unsupported candidate may fail while another candidate succeeds; such
failures and all candidate warnings are reported. Protocol violations
always stop fitting. The built-in univariateML family explicitly rejects
observation weights.

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

# set up vine copula model with Gaussian margins
vc <- vine_dist(list(stats_margin("norm")), pcs, mat)

# show model
summary(vc)
#> $margins
#> # A data.frame: 3 x 2 
#>  margin distr
#>       1  norm
#>       2  norm
#>       3  norm
#> 
#> $copula
#> # A data.frame: 3 x 10 
#>  tree edge conditioned conditioning var_types family rotation parameters df
#>     1    1        3, 1                    c,c    bb1       90       3, 2  2
#>     1    2        2, 1                    c,c    bb1       90       3, 2  2
#>     2    1        3, 2            1       c,c    bb1       90       3, 2  2
#>   tau
#>  -0.8
#>  -0.8
#>  -0.8
#> 

# simulate some data
x <- rvine(50, vc)

# estimate a vine copula model
fit <- vine(x, copula_controls = list(family_set = "par"))
summary(fit)
#> $margins
#> # A data.frame: 3 x 8 
#>  margin name family type xmin xmax npars loglik
#>       1   V1  kde1d    c -Inf  Inf   4.6    -59
#>       2   V2  kde1d    c -Inf  Inf   4.2    -60
#>       3   V3  kde1d    c -Inf  Inf   5.7    -57
#> 
#> $copula
#> # A data.frame: 3 x 11 
#>  tree edge conditioned conditioning var_types family rotation parameters df
#>     1    1        3, 1                    c,c gumbel      270        4.2  1
#>     1    2        2, 1                    c,c gumbel      270          4  1
#>     2    1        3, 2            1       c,c gumbel      270        4.1  1
#>    tau loglik
#>  -0.76     52
#>  -0.75     49
#>  -0.75     50
#> 

# parametric margin selection (requires the suggested univariateML package)
if (requireNamespace("univariateML", quietly = TRUE)) {
  fit_par_margins <- vine(
    x,
    margins_controls = list(family_set = c("norm", "cauchy"))
  )
}
#> Loading required namespace: intervals

## model for discrete data
x <- as.data.frame(x)
x[, 1] <- ordered(round(x[, 1]), levels = seq.int(-5, 5))
fit_disc <- vine(x, copula_controls = list(family_set = "par"))
summary(fit_disc)
#> $margins
#> # A data.frame: 3 x 8 
#>  margin name family type xmin xmax npars loglik
#>       1   V1  kde1d    d    0   10   5.2    -69
#>       2   V2  kde1d    c -Inf  Inf   4.2    -60
#>       3   V3  kde1d    c -Inf  Inf   5.7    -57
#> 
#> $copula
#> # A data.frame: 3 x 11 
#>  tree edge conditioned conditioning var_types family rotation parameters df
#>     1    1        3, 1                    c,d  frank        0        -16  1
#>     1    2        2, 1                    c,d  frank        0        -10  1
#>     2    1        3, 2            1       c,c      t        0 0.25, 2.00  2
#>    tau loglik
#>  -0.78   31.4
#>  -0.67   21.3
#>   0.16    3.4
#> 
```
