# Bivariate copula models

Bivariate copulas are the building blocks of vine models. rvinecopulib
can construct a known model, fit or select one from data, evaluate its
distribution functions, and summarize its dependence.

## Implemented families

| Type          | Family                | Identifier      |      Parameters |
| ------------- | --------------------- | --------------- | --------------: |
| Independence  | Independence          | `"indep"`       |               0 |
| Elliptical    | Gaussian              | `"gaussian"`    |               1 |
| Elliptical    | Student t             | `"t"`           |               2 |
| Archimedean   | Clayton               | `"clayton"`     |               1 |
| Archimedean   | Gumbel                | `"gumbel"`      |               1 |
| Archimedean   | Frank                 | `"frank"`       |               1 |
| Archimedean   | Joe                   | `"joe"`         |               1 |
| Archimedean   | BB1, BB6, BB7, BB8    | lower-case name |               2 |
| Extreme value | Tawn                  | `"tawn"`        |               3 |
| Nonparametric | Transformation kernel | `"tll"`         | effective count |

Unambiguous partial names such as `"gauss"` and collection aliases such
as `"par"` are accepted. For reusable code, full names make intent
clearer.

## Construct a known model

`bicop_dist()` specifies a family, rotation, parameters, and variable
types.

``` r
cop <- bicop_dist(
  family = "clayton",
  rotation = 90,
  parameters = 2
)
cop
#> Bivariate copula ('bicop_dist'): family = clayton, rotation = 90, parameters = 2, var_types = c,c
```

Rotations are `0`, `90`, `180`, or `270` degrees. They allow asymmetric
Archimedean copulas to represent different corners and negative
dependence. Elliptical, Frank, independence, and nonparametric copulas
already cover both signs without rotations.

## Evaluate and simulate

The density, CDF, and random generator follow R’s `d*`, `p*`, and `r*`
conventions.

``` r
points <- rbind(c(0.2, 0.7), c(0.8, 0.3))
dbicop(points, cop)
#> [1] 1.562211 1.901324
pbicop(points, cop)
#> [1] 0.08022147 0.13123681
rbicop(4, cop)
#>            [,1]       [,2]
#> [1,] 0.77133621 0.46791844
#> [2,] 0.79602646 0.08946854
#> [3,] 0.07922695 0.35240165
#> [4,] 0.96904452 0.02742706
```

The same functions also accept family, rotation, and parameters
directly:

``` r
dbicop(points, "gaussian", 0, 0.6)
#> [1] 0.6267684 0.6267684
```

## Conditional distributions and h-functions

H-functions are the conditional distributions used recursively in vine
copulas. If `cond_var = 1`, the first input variable is conditioned on
and the function returns the conditional CDF of the second. With
`inverse = TRUE`, it returns the corresponding conditional quantile.

``` r
h <- hbicop(points, cond_var = 1, family = cop)
h
#> [1] 0.4649857 0.6008183
hbicop(cbind(points[, 1], h), cond_var = 1, family = cop, inverse = TRUE)
#> [1] 0.7 0.3
```

This round trip recovers the second column up to numerical precision.
The same conditional-distribution operations drive vine recursion and
simulation; see [Conditional simulation and Rosenblatt
transforms](https://vinecopulib.github.io/rvinecopulib/articles/conditional-simulation.md)
for the multivariate workflow.

## Dependence measures

`par_to_ktau()` converts family parameters to Kendall’s tau, and
`ktau_to_par()` provides the inverse for families where tau determines
all parameters. Model objects expose Kendall’s tau through `get_ktau()`.

``` r
par_to_ktau("clayton", 90, 2)
#> [1] -0.5
get_ktau(cop)
#> [1] -0.5
blomqvist_beta(cop)
#> [1] -0.5118579
tail_dep(cop)
#>          variable2
#> variable1     lower upper
#>     lower 0.0000000     0
#>     upper 0.7071068     0
```

`tail_dep()` reports the four corner tail-dependence coefficients. This
is especially useful for rotated and asymmetric families, where a single
lower/upper pair is not enough to describe the model.

## Fit and select from data

`bicop()` expects two approximately uniform columns. It fits every
requested compatible family and selects the best according to `selcrit`.

``` r
u <- rbicop(300, "gumbel", 180, 2)
fit <- bicop(
  u,
  family_set = c("gaussian", "clayton", "gumbel", "frank"),
  selcrit = "bic"
)
fit
#> Bivariate copula fit ('bicop'): family = gumbel, rotation = 180, parameters = 2.16, var_types = c,c
summary(fit)
#> Bivariate copula fit ('bicop'): family = gumbel, rotation = 180, parameters = 2.16, var_types = c,c
#> Dependence: tau = 0.54; beta = 0.54; tail dependence: LL = 0.62
#> Fit: n = 300; logLik = 126.44; df = 1.00; AIC = -250.87; BIC = -247.17
```

Useful family collections include:

  - `"all"`, `"parametric"`, and `"nonparametric"`;
  - `"oneparametric"`, `"twoparametric"`, and `"threeparametric"`;
  - `"elliptical"`, `"archimedean"`, `"ev"`, and `"bbs"`;
  - `"itau"`, the families compatible with Kendall’s tau inversion.

Partial matching makes shorter forms such as `"onepar"` and `"par"`
available. Explicit character vectors are preferable when the candidate
set is part of a scientific specification.

Maximum likelihood is the default for parametric families. Set
`par_method = "itau"` for Kendall’s tau inversion. The
transformation-kernel model uses `nonpar_method` and `mult` to control
the local-likelihood order and smoothing.

## Observation weights and missing values

Pass one nonnegative weight per row with `weights`. The backend
standardizes weights internally. Missing or zero-weight observations do
not contribute to a pair fit.

``` r
weights <- seq_len(nrow(u))
weighted_fit <- bicop(
  u,
  family_set = c("gaussian", "gumbel"),
  weights = weights
)
weighted_fit
#> Bivariate copula fit ('bicop'): family = gumbel, rotation = 180, parameters = 2.1, var_types = c,c
```

## Discrete variables

A discrete copula observation needs both `F(x)` and `F(x-)`. For two
discrete variables, pass four columns: the two CDF values followed by
their two left limits. Redundant left-limit columns for continuous
variables may be omitted.

``` r
counts <- cbind(rpois(80, 2), rpois(80, 3))
u_discrete <- cbind(
  ppois(counts[, 1], 2),
  ppois(counts[, 2], 3),
  ppois(counts[, 1] - 1, 2),
  ppois(counts[, 2] - 1, 3)
)
fit_discrete <- bicop(
  u_discrete,
  var_types = c("d", "d"),
  family_set = "onepar"
)
fit_discrete
#> Bivariate copula fit ('bicop'): family = joe, rotation = 270, parameters = 1.05, var_types = d,d
```

The [discrete-data
guide](https://vinecopulib.github.io/rvinecopulib/articles/discrete-data.md)
derives the likelihood contribution and covers compact layouts, mixed
models, and atoms in detail.

## Observation-specific parameters

Evaluation and simulation can use one parameter set per observation
without constructing many model objects. A one-parameter family accepts
a vector; multi-parameter families use a matrix with one row per
observation.

``` r
rho <- seq(-0.7, 0.7, length.out = nrow(points))
dbicop(points, "gaussian", 0, parameters = rho)
#> [1] 1.6000932 0.4764093
rbicop(length(rho), "gaussian", 0, parameters = rho)
#>            [,1]      [,2]
#> [1,] 0.56361393 0.3312878
#> [2,] 0.01234443 0.1489253
```

This interface is intended for conditional or covariate-dependent copula
models. Parameters are not recycled, and vectorized evaluation is
limited to parametric families.

## Related documentation

  - [Vine copula models and
    structures](https://vinecopulib.github.io/rvinecopulib/articles/vine-copula-models.md)
    shows how bivariate copulas are assembled into a high-dimensional
    model.
  - [Scores, Hessians, and varying
    parameters](https://vinecopulib.github.io/rvinecopulib/articles/likelihood-inference.md)
    documents likelihood derivatives and parameter ordering.
  - The [`bicop()`
    reference](https://vinecopulib.github.io/rvinecopulib/reference/bicop.md)
    lists all fitting controls; [`bicop_dist()` and its distribution
    functions](https://vinecopulib.github.io/rvinecopulib/reference/bicop_dist.md)
    document construction and evaluation.
