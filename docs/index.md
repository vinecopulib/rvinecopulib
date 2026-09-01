# rvinecopulib

[![R-CMD-check](https://github.com/vinecopulib/rvinecopulib/actions/workflows/R-CMD-check.yaml/badge.svg?branch=main)](https://github.com/vinecopulib/rvinecopulib/actions/workflows/R-CMD-check.yaml)
[![Coverage
Status](https://img.shields.io/codecov/c/github/vinecopulib/rvinecopulib/main.svg)](https://app.codecov.io/github/vinecopulib/rvinecopulib?branch=main)
[![CRAN
version](https://www.r-pkg.org/badges/version/rvinecopulib)](https://cran.r-project.org/package=rvinecopulib)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/rvinecopulib)](https://cran.r-project.org/package=rvinecopulib)

rvinecopulib provides high-performance tools for bivariate and vine
copula models. It covers model construction, estimation and selection,
simulation, prediction, visualization, discrete and mixed data, and full
multivariate distributions with fitted margins. The package is the R
interface to the
[vinecopulib](https://github.com/vinecopulib/vinecopulib) C++ library.

## Main capabilities

  - fit, select, evaluate, and simulate [bivariate copula
    models](https://vinecopulib.github.io/rvinecopulib/articles/bivariate-copulas.html),
    including parametric, nonparametric, rotated, and extreme-value
    families;
  - fit and select high-dimensional [vine copula models and
    structures](https://vinecopulib.github.io/rvinecopulib/articles/vine-copula-models.html),
    with automatic family, structure, truncation, and threshold
    selection;
  - combine vine copulas with nonparametric, parametric, or custom
    [marginal
    models](https://vinecopulib.github.io/rvinecopulib/articles/marginal-models.html)
    to form complete multivariate distributions;
  - handle continuous, integer-valued discrete, mixed, and
    [zero-inflated
    data](https://vinecopulib.github.io/rvinecopulib/articles/discrete-data.html),
    with optional observation weights;
  - perform [conditional simulation and Rosenblatt
    transforms](https://vinecopulib.github.io/rvinecopulib/articles/conditional-simulation.html);
  - compute likelihood [scores and
    Hessians](https://vinecopulib.github.io/rvinecopulib/articles/likelihood-inference.html)
    and evaluate models with observation-specific parameters;
  - use parallel fitting, custom marginal families, and custom
    tree-selection criteria in extensible modeling workflows.

## Installation

Install the stable release from CRAN:

``` r
install.packages("rvinecopulib")
```

Install the development version from GitHub:

``` r
remotes::install_github("vinecopulib/rvinecopulib")
```

## Examples

rvinecopulib exposes the same modeling framework at three levels:

| Starting point                       | Main function | Model                  |
| ------------------------------------ | ------------- | ---------------------- |
| Two uniform variables                | `bicop()`     | One bivariate copula   |
| Uniform pseudo-observations          | `vinecop()`   | Dependence only        |
| Observations on their original scale | `vine()`      | Margins and dependence |

### A bivariate copula

`bicop()` fits and selects a copula model for two uniform variables.
Fixed models can be created with `bicop_dist()`.

``` r
u <- rbicop(200, family = "clayton", rotation = 90, parameters = 2)
bivariate_fit <- bicop(u, family_set = "par")
summary(bivariate_fit)

dbicop(u[1:5, ], bivariate_fit)
tail_dep(bivariate_fit)
```

See the [bivariate-copula
article](https://vinecopulib.github.io/rvinecopulib/articles/bivariate-copulas.html)
for implemented families, rotations, h-functions, dependence measures,
and selection controls.

### A vine copula on the copula scale

`pseudo_obs()` converts continuous observations to approximately uniform
scores. `vinecop()` then selects the vine structure, pair-copula
families, and parameters.

``` r
u <- pseudo_obs(as.matrix(USArrests))
copula_fit <- vinecop(u, family_set = "onepar")
summary(copula_fit)

simulated_u <- rvinecop(100, copula_fit)
dvinecop(simulated_u[1:5, ], copula_fit)
```

See the [vine-copula
article](https://vinecopulib.github.io/rvinecopulib/articles/vine-copula-models.html)
for structure construction, selection, truncation, and model inspection.

### A full distribution on the original scale

`vine()` fits one marginal distribution per variable and a vine copula
to their probability integral transforms. The default margins are
nonparametric; parametric and custom marginal families are also
supported.

``` r
n <- 150
latent <- rnorm(n)
x <- data.frame(
  amount = exp(latent + rnorm(n, sd = 0.5)),
  duration = exp(0.5 * latent + rnorm(n, sd = 0.7)),
  count = rpois(n, exp(0.2 + 0.3 * latent))
)

fit <- vine(
  x,
  var_types = c("c", "c", "d"),
  copula_controls = list(family_set = "onepar")
)
summary(fit)
rvine(5, fit)
```

See the [marginal-modeling
article](https://vinecopulib.github.io/rvinecopulib/articles/marginal-models.html)
for parametric selection, custom families, ordered variables, zero
inflation, and observation weights. The [getting-started
article](https://vinecopulib.github.io/rvinecopulib/articles/getting-started.html)
develops the complete workflow.

The constructors `bicop_dist()`, `vinecop_dist()`, and `vine_dist()`
create models from components specified directly.

The [complete API
reference](https://vinecopulib.github.io/rvinecopulib/reference/) and
all articles are available on the [package
website](https://vinecopulib.github.io/rvinecopulib/). Questions and bug
reports are welcome in the [GitHub issue
tracker](https://github.com/vinecopulib/rvinecopulib/issues).
