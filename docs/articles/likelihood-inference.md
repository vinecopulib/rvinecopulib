# Scores, Hessians, and varying parameters

rvinecopulib exposes observation-wise likelihood scores and average
Hessians for parametric bivariate and vine copula models. These are the
basic inputs for convergence checks, standard errors, sandwich
covariance estimates, and gradient-based extensions.

For a parameter vector `theta` and observation-wise log-likelihood
contributions `ell_i(theta)`, rvinecopulib reports
\[s_{i}(\theta) = \frac{\partial\ell_{i}(\theta)}{\partial\theta},\qquad\bar{H}(\theta) = \frac{1}{n}\sum\limits_{i = 1}^{n}\frac{\partial^{2}\ell_{i}(\theta)}{\partial\theta\partial\theta^{\top}}.\]
Thus `scores()` contains the row vectors `s_i`, while `hessian()`
returns `H-bar`.

## Bivariate scores and Hessians

Fit a parametric model and evaluate derivatives at its fitted
parameters:

``` r
u2 <- rbicop(250, "gaussian", 0, 0.6)
fit2 <- bicop(u2, family_set = "gaussian")

score2 <- scores(u2, fit2)
hessian2 <- hessian(u2, fit2)
dim(score2)
#> [1] 250   1
hessian2
#>           [,1]
#> [1,] -4.204222
colSums(score2)
#> [1] -0.0003204515
```

The score matrix has one row per observation and one column per model
parameter. At an interior maximum-likelihood estimate, its column sums
should be close to zero. The Hessian is averaged over observations;
multiply by the sample size when a summed log-likelihood Hessian is
required.

These derivatives are available for continuous parametric models.
Boundary estimates and weakly identified parameters need the same care
as in any likelihood analysis.

## Vine parameter ordering

For a vine, score columns follow `(tree, edge, parameter)` order: tree 1
from left to right, then tree 2, and so on; multi-parameter pair copulas
contribute their parameters in the order shown by `get_parameters()`.

``` r
n <- 220
z <- rnorm(n)
x <- cbind(
  z + rnorm(n),
  0.7 * z + rnorm(n),
  -0.5 * z + rnorm(n)
)
u <- pseudo_obs(x)
fit <- vinecop(
  u,
  structure = dvine_structure(1:3),
  family_set = "gaussian"
)
get_all_parameters(fit)
#> Nested list of lists for the copula parameters of a 3 dimensional vine with all trees: 
#> - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#> x[[1]] -> a list with the 2 copula parameters of tree 1. 
#> x[[2]] -> a list with the 1 copula parameters of tree 2.
```

The number and ordering of derivative columns can be checked against the
edge summary.

``` r
score <- scores(u, fit)
hess <- hessian(u, fit)
dim(score)
#> [1] 220   3
dim(hess)
#> [1] 3 3
summary(fit)
#> # A data.frame: 3 x 11 
#>  tree edge conditioned conditioning var_types   family rotation parameters df
#>     1    1        1, 2                    c,c gaussian        0        0.4  1
#>     1    2        2, 3                    c,c gaussian        0      -0.33  1
#>     2    1        1, 3            2       c,c gaussian        0      -0.22  1
#>    tau loglik
#>   0.26     18
#>  -0.22     12
#>  -0.14      5
```

## Stepwise and full likelihood derivatives

Vine models are usually fitted sequentially. With `step_wise = TRUE`
(the default), each pair copula treats the pseudo-observations passed
down from earlier trees as fixed. This is the objective optimized by the
sequential estimator.

With `step_wise = FALSE`, derivatives propagate through the complete
h-function cascade and refer to the full joint log-likelihood.

For a vine with edge set `E`, an observation contributes
\[\ell_{i}(\theta) = \sum\limits_{e \in E}\log c_{e}\{ u_{e,1,i}(\theta_{< e}),u_{e,2,i}(\theta_{< e});\theta_{e}\}.\]
Stepwise derivatives hold the conditional pseudo-observations `u_e`
fixed. Full derivatives also differentiate their dependence on
parameters in earlier trees, denoted by `theta_<e` above.

``` r
stepwise_gradient <- colSums(scores(u, fit, step_wise = TRUE))
full_gradient <- colSums(scores(u, fit, step_wise = FALSE))
rbind(stepwise = stepwise_gradient, full = full_gradient)
#>                   [,1]          [,2]         [,3]
#> stepwise -2.262873e-05 -2.705407e-06 6.958963e-05
#> full      3.299510e-01 -4.535448e-01 6.958963e-05
```

The stepwise gradient should be close to zero at a sequential fit. The
full gradient generally is not: a sequence of conditional maximizations
is not the same as a joint maximum. This distinction should be stated
explicitly whenever derivatives are used for inference.

## Reuse density intermediates

`dvinecop(..., keep_all = TRUE)` returns the density and the per-edge
values computed by the h-function cascade.

``` r
full <- dvinecop(u[1:5, ], fit, keep_all = TRUE)
names(full)
#> [1] "pdf"        "pdf_edges"  "hfunc1"     "hfunc2"     "hfunc1_sub"
#> [6] "hfunc2_sub"
full$pdf
#> [1] 3.1267573 1.4497370 1.0291963 1.1218710 0.3357038
```

The triangular `pdf_edges`, `hfunc1`, and `hfunc2` entries follow the
same tree and edge ordering as the fitted pair copulas. The `_sub`
h-functions contain left-limit evaluations for discrete models and are
empty for fully continuous models.

## Observation-specific bivariate parameters

Parametric copulas can be evaluated with a different parameter set in
every row. This is useful when another model maps covariates to valid
copula parameters.

``` r
points <- matrix(runif(200), ncol = 2)
rho <- seq(-0.8, 0.8, length.out = nrow(points))

density <- dbicop(points, "gaussian", 0, parameters = rho)
score_varying <- scores(points, bicop_dist("gaussian"), parameters = rho)
hessian_varying <- hessian(
  points,
  bicop_dist("gaussian"),
  parameters = rho
)
head(cbind(rho, density, score_varying))
#>             rho      density           
#> [1,] -0.8000000 0.3491747073  9.9373957
#> [2,] -0.7838384 1.4793496498 -1.4969973
#> [3,] -0.7676768 0.0000346894 58.7386032
#> [4,] -0.7515152 0.6343731400  2.9922591
#> [5,] -0.7353535 1.4060226704 -0.7165023
#> [6,] -0.7191919 0.7572296856  1.7017953
```

For a one-parameter family, `parameters` may be a vector. Otherwise it
must be a matrix with one row per observation and one column per family
parameter. Parameters are not recycled.

## Observation-specific vine parameters

The vine interface takes an `n` by `p` parameter matrix, where `p` is
the total number of pair-copula parameters and columns use the score
ordering.

``` r
base_parameters <- unlist(get_all_parameters(fit), use.names = FALSE)
parameter_matrix <- matrix(
  rep(base_parameters, each = nrow(u)),
  nrow = nrow(u)
)
parameter_matrix[, 1] <- seq(-0.6, 0.6, length.out = nrow(u))

density_varying <- dvinecop(u, fit, parameters = parameter_matrix)
score_varying <- scores(u, fit, parameters = parameter_matrix)
head(density_varying)
#> [1] 0.2906631 0.4834200 0.5210400 1.5026325 0.2022223 0.8928398
dim(score_varying)
#> [1] 220   3
```

Observation-specific vine parameters are supported for continuous,
all-parametric models. Each row must represent a valid complete
parameter vector. The model stored in `fit` is not modified.

## From derivatives to uncertainty estimates

The score and Hessian outputs deliberately provide ingredients rather
than a single universal covariance estimator. The correct assembly
depends on whether the estimator is stepwise or joint, whether weights
or clusters are present, and whether robust inference is desired.

For an ordinary interior bivariate MLE, the inverse negative summed
Hessian is the familiar model-based covariance estimate. More elaborate
vine inference should account for the sequential estimation scheme
instead of treating the stepwise estimate as a joint maximum.

In the ordinary independent-observation case, define
\[A = - \bar{H},\qquad B = \frac{1}{n}\sum\limits_{i = 1}^{n}s_{i}s_{i}^{\top}.\]
The model-based covariance estimate is `(-n H-bar)^(-1)`, while the
usual sandwich estimate is `A^(-1) B A^(-1) / n`. Weights, clusters, and
sequential estimation change how these pieces should be assembled.

## Related documentation

  - [Bivariate copula
    models](https://vinecopulib.github.io/rvinecopulib/articles/bivariate-copulas.md)
    and [vine copula models and
    structures](https://vinecopulib.github.io/rvinecopulib/articles/vine-copula-models.md)
    provide the model definitions used in the examples.
  - [Conditional simulation and Rosenblatt
    transforms](https://vinecopulib.github.io/rvinecopulib/articles/conditional-simulation.md)
    explains the h-function cascade that full derivatives propagate
    through.
  - The [bivariate distribution
    reference](https://vinecopulib.github.io/rvinecopulib/reference/bicop_methods.md)
    and [vine-copula distribution
    reference](https://vinecopulib.github.io/rvinecopulib/reference/vinecop_methods.md)
    specify derivative support and parameter shapes.
