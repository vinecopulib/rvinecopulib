# Getting started with rvinecopulib

rvinecopulib has three main modeling interfaces. The right one depends
on the scale of the data and whether the marginal distributions are part
of the model.

| Starting point              | Function    | Result                    |
| --------------------------- | ----------- | ------------------------- |
| Original observations       | `vine()`    | Margins and a vine copula |
| Uniform pseudo-observations | `vinecop()` | A vine copula             |
| Two uniform variables       | `bicop()`   | A bivariate copula        |

The fitted classes add fit information to distribution classes. A `vine`
object is also a `vine_dist`, a `vinecop` is also a `vinecop_dist`, and
a `bicop` is also a `bicop_dist`. Simulation and evaluation therefore
use the same functions for fitted and manually constructed models.

## A complete model on the data scale

Use `vine()` when the observations are on their original scale. It
estimates a marginal distribution for every variable, transforms the
data to the copula scale, and fits the dependence model.

``` r
n <- 120
z <- rnorm(n)
x <- data.frame(
  amount = exp(z + rnorm(n, sd = 0.5)),
  duration = exp(0.5 * z + rnorm(n, sd = 0.8)),
  events = rpois(n, exp(0.3 + 0.25 * z))
)

fit <- vine(
  x,
  var_types = c("c", "c", "d"),
  copula_controls = list(family_set = "onepar"),
  keep_data = TRUE
)
fit
#> 3-dimensional vine distribution fit ('vine')
#> nobs = 120   logLik = -484.16   npars = 24.21   AIC = 1016.74   BIC = 1084.22
summary(fit)
#> $margins
#> # A data.frame: 3 x 8 
#>  margin     name family type xmin xmax npars loglik
#>       1   amount  kde1d    c -Inf  Inf   7.9   -186
#>       2 duration  kde1d    c -Inf  Inf   9.2   -141
#>       3   events  kde1d    d -Inf  Inf   4.1   -175
#> 
#> $copula
#> # A data.frame: 3 x 11 
#>  tree edge conditioned conditioning var_types   family rotation parameters df
#>     1    1        3, 2                    d,c    frank        0        1.2  1
#>     1    2        2, 1                    c,c gaussian        0       0.46  1
#>     2    1        3, 1            2       d,c  clayton        0       0.05  1
#>    tau loglik
#>  0.133  2.472
#>  0.307 14.530
#>  0.024  0.093
#> 
#> ---
#> 1 <-> amount,   2 <-> duration,   3 <-> events
```

Integer-valued discrete columns are declared with `var_types = "d"`.
Ordered factors are detected automatically. See [Marginal
models](https://vinecopulib.github.io/rvinecopulib/articles/marginal-models.md)
for the default KDE margins, parametric selection, custom margin
families, zero inflation, and observation weights.

The fitted object represents the joint distribution. Its density or mass
can be evaluated on the original data scale, and new observations can be
drawn with their original column names and types.

``` r
dvine(x[1:4, ], fit)
#> [1] 0.04863737 0.06268189 0.07279592 0.03497973
predict(fit, x[1:4, ], what = "pdf")
#> [1] 0.04863737 0.06268189 0.07279592 0.03497973
logLik(fit)
#> 'log Lik.' -484.1607 (df=24.20907)
AIC(fit)
#> [1] 1016.739

simulated <- rvine(5, fit)
simulated
#>         amount  duration events
#> [1,] 0.8176936 0.5882763      1
#> [2,] 0.2327406 0.7762025      1
#> [3,] 1.3789696 2.3025974      0
#> [4,] 8.2736166 2.1162245      2
#> [5,] 1.2403309 0.7031729      0
```

`pvine()` and `predict(..., what = "cdf")` approximate a multivariate
CDF by Monte Carlo. Increase `n_mc` when accurate tail probabilities are
important.

## A dependence model on the copula scale

Use `vinecop()` when the margins are not part of the model. Its input
consists of approximately uniform pseudo-observations. For continuous
data, `pseudo_obs()` applies scaled ranks column by column.

``` r
u <- pseudo_obs(as.matrix(USArrests))
copula_fit <- vinecop(
  u,
  family_set = c("gaussian", "clayton", "gumbel", "frank")
)
summary(copula_fit)
#> # A data.frame: 6 x 11 
#>  tree edge conditioned conditioning var_types   family rotation parameters df
#>     1    1        3, 4                    c,c   gumbel      180        1.5  1
#>     1    2        4, 2                    c,c   gumbel      180        2.1  1
#>     1    3        2, 1                    c,c   gumbel      180        2.5  1
#>     2    1        3, 2            4       c,c gaussian        0     -0.034  1
#>     2    2        4, 1            2       c,c  clayton        0       0.52  1
#>     3    1        3, 1         2, 4       c,c  clayton      270       0.45  1
#>     tau loglik
#>   0.318  5.944
#>   0.525 17.882
#>   0.594 24.983
#>  -0.022  0.027
#>   0.206  4.308
#>  -0.183  2.504
#> 
#> ---
#> 1 <-> Murder,   2 <-> Assault,   3 <-> UrbanPop,   4 <-> Rape
```

The result can be evaluated and simulated directly on the unit
hypercube.

``` r
new_u <- rvinecop(6, copula_fit)
dvinecop(new_u, copula_fit)
#> [1] 6.4235884 1.2890814 0.7926246 6.6280552 1.0875679 2.3087096
pvinecop(new_u[1:2, ], copula_fit, n_mc = 1000)
#> [1] 0.070 0.409
```

Do not apply ordinary ranks blindly to discrete data. Discrete copula
models need both `F(x)` and its left limit `F(x-)`; the dedicated
discrete-data article explains the accepted layouts. `vine()` constructs
these quantities automatically from its fitted margins. See [Vine copula
models and
structures](https://vinecopulib.github.io/rvinecopulib/articles/vine-copula-models.md)
for structure selection, truncation, thresholding, and model inspection.

## A bivariate model

Use `bicop()` for two-dimensional copula data. It estimates each
compatible candidate and selects a family using AIC by default.

``` r
u2 <- rbicop(200, "clayton", 90, 2)
bivariate_fit <- bicop(u2, family_set = "par")
bivariate_fit
#> Bivariate copula fit ('bicop'): family = clayton, rotation = 90, parameters = 2.11, var_types = c,c

dbicop(u2[1:4, ], bivariate_fit)
#> [1]  1.329958 58.723070  2.498366  1.461224
tail_dep(bivariate_fit)
#>          variable2
#> variable1     lower upper
#>     lower 0.0000000     0
#>     upper 0.7195585     0
```

The [bivariate-model
guide](https://vinecopulib.github.io/rvinecopulib/articles/bivariate-copulas.md)
covers rotations, h-functions, dependence measures, family sets,
vectorized parameters, and nonparametric models.

## Fixed and fitted models

The `*_dist()` constructors are useful when the model components are
known. For example, a Gaussian copula with correlation 0.7 is a complete
bivariate distribution even though it has not been fitted:

``` r
fixed <- bicop_dist("gaussian", parameters = 0.7)
dbicop(c(0.25, 0.8), fixed)
#> [1] 0.3674023
rbicop(3, fixed)
#>            [,1]       [,2]
#> [1,] 0.04373282 0.07634724
#> [2,] 0.14755744 0.19933630
#> [3,] 0.26463486 0.46614350
```

Fitted objects additionally store a log-likelihood, observation count,
and fit controls. Set `keep_data = TRUE` only when methods such as
`fitted()` need the training data.

## Selection controls

The most frequently used controls are:

  - `family_set`: candidate pair-copula families or a predefined
    collection;
  - `selcrit`: `"loglik"`, `"aic"`, `"bic"`, `"mbic"`, or `"mbicv"`
    where applicable;
  - `par_method`: maximum likelihood (`"mle"`) or Kendall’s tau
    inversion (`"itau"`) for supported families;
  - `allow_rotations`: whether rotated Archimedean families may be
    selected;
  - `weights`: optional observation weights;
  - `cores`: parallel workers used by supported operations.

`vinecop()` additionally controls structure selection, truncation, and
thresholding. These options are developed in the vine-model article
rather than hidden in a single large example.

## Reproducibility and persistence

Ordinary model fitting is deterministic. Simulation uses R’s
random-number state, so `set.seed()` makes it reproducible. Quasi-random
simulation is available through `qrng = TRUE`.

All model objects can be stored with `saveRDS()` and restored with
`readRDS()`. Custom margin classes remain usable after restoration as
long as their S3 methods are available in the R session.

``` r
path <- tempfile(fileext = ".rds")
saveRDS(copula_fit, path)
restored <- readRDS(path)
unlink(path)

set.seed(102)
draws1 <- rvinecop(4, copula_fit)
set.seed(102)
draws2 <- rvinecop(4, restored)
stopifnot(isTRUE(all.equal(draws1, draws2)))
```

## Where to go next

  - [Marginal
    models](https://vinecopulib.github.io/rvinecopulib/articles/marginal-models.md)
    develops full distributions on the original data scale.
  - [Bivariate copula
    models](https://vinecopulib.github.io/rvinecopulib/articles/bivariate-copulas.md)
    documents the pair-copula building blocks and family-selection
    controls.
  - [Vine copula models and
    structures](https://vinecopulib.github.io/rvinecopulib/articles/vine-copula-models.md)
    covers high-dimensional dependence modeling.
  - [Discrete
    data](https://vinecopulib.github.io/rvinecopulib/articles/discrete-data.md),
    [conditional
    simulation](https://vinecopulib.github.io/rvinecopulib/articles/conditional-simulation.md),
    and [likelihood
    inference](https://vinecopulib.github.io/rvinecopulib/articles/likelihood-inference.md)
    introduce advanced workflows.
  - The [function
    reference](https://vinecopulib.github.io/rvinecopulib/reference/index.md)
    lists the complete API by modeling task.
