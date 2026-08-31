# Discrete, mixed, and zero-inflated data

Copula likelihoods for variables with atoms require two probabilities at
every observation: the CDF value `F(x)` and its left limit `F(x-)`. For
an integer-valued variable, `F(x-) = F(x - 1)`. Continuous variables
satisfy `F(x-) = F(x)`.

The high-level `vine()` interface computes these values from fitted
margins. The copula-only `bicop()` and `vinecop()` interfaces accept
them explicitly.

## Discrete variables on the original scale

All ordinary discrete variables are assumed to have integer support.
Declare a numeric or integer column through `var_types = "d"`.

``` r
n <- 100
latent <- rnorm(n)
x <- data.frame(
  value = latent + rnorm(n),
  count = rpois(n, exp(0.5 + 0.3 * latent))
)

fit <- vine(
  x,
  var_types = c("c", "d"),
  copula_controls = list(family_set = "onepar")
)
summary(fit)$margins
#> # A data.frame: 2 x 8 
#>  margin  name family type xmin xmax npars loglik
#>       1 value  kde1d    c -Inf  Inf   5.6   -170
#>       2 count  kde1d    d -Inf  Inf   3.8   -156
dvine(x[1:5, ], fit)
#> [1] 0.086337749 0.070854020 0.005353429 0.024217577 0.096298044
```

Explicitly discrete columns are checked for finite integer-valued
observations. This catches accidental declarations of rounded-looking
continuous data.

## Ordered factors

Ordered factors are discrete variables whose levels are mapped
internally to the integer support `0, 1, ...`. Their class and levels
are restored after simulation.

``` r
survey <- data.frame(
  rating = ordered(
    sample(c("low", "middle", "high"), n, replace = TRUE),
    levels = c("low", "middle", "high")
  ),
  score = rnorm(n)
)

ordered_fit <- vine(
  survey,
  copula_controls = list(family_set = "indep")
)
ordered_draws <- rvine(5, ordered_fit)
str(ordered_draws)
#> 'data.frame':    5 obs. of  2 variables:
#>  $ rating: Ord.factor w/ 3 levels "low"<"middle"<..: 3 3 1 3 1
#>  $ score : num  -2.83583 -0.00419 -0.99283 -0.23898 0.62379
```

No `var_types` entry is needed because `ordered()` carries the
declaration in the column class.

## Zero-inflated variables

`zero_inflated()` marks a continuous variable with an atom at zero. It
behaves like a numeric vector in a data frame and survives ordinary
subsetting.

``` r
amount <- zero_inflated(c(rep(0, 25), rexp(75)))
zi_data <- data.frame(amount = amount, score = rnorm(100))
inherits(zi_data$amount, "zero_inflated")
#> [1] TRUE

zi_fit <- vine(
  zi_data,
  copula_controls = list(family_set = "indep")
)
summary(zi_fit)$margins[, c("name", "type")]
#> # A data.frame: 2 x 2 
#>    name type
#>  amount   zi
#>   score    c
```

The atom is always located at zero. A custom zero-inflated margin must
return the atom probability from `dmargin(0, margin)`. rvinecopulib then
computes the left limit as

    pmargin(0, margin) - dmargin(0, margin)

and uses the ordinary CDF away from zero. Alternatively, declare the
type with `var_types = "zi"`.

## Copula-scale data layouts

Let `d` be the dimension and `k` the number of discrete variables.
Copula methods accept two layouts:

  - **expanded:** `2d` columns, with all `F(x)` columns followed by all
    `F(x-)` columns;
  - **compact:** `d + k` columns, with the `F(x)` block followed only by
    left-limit columns for discrete variables, in variable order.

Consider a three-dimensional model with discrete variables 1 and 3.

``` r
n <- 80
raw <- data.frame(
  first = rpois(n, 2),
  middle = rnorm(n),
  last = rpois(n, 4)
)

values <- cbind(
  ppois(raw$first, 2),
  pnorm(raw$middle),
  ppois(raw$last, 4)
)
left_discrete <- cbind(
  ppois(raw$first - 1, 2),
  ppois(raw$last - 1, 4)
)
compact <- cbind(values, left_discrete)

left_all <- cbind(
  ppois(raw$first - 1, 2),
  pnorm(raw$middle),
  ppois(raw$last - 1, 4)
)
expanded <- cbind(values, left_all)

copula_fit <- vinecop(
  compact,
  var_types = c("d", "c", "d"),
  family_set = "indep"
)
stopifnot(isTRUE(all.equal(
  dvinecop(compact, copula_fit),
  dvinecop(expanded, copula_fit)
)))
```

For a mixed bivariate model, this means three columns: the two CDF
values and the left limit of the discrete variable. Fully discrete
bivariate models need four columns.

## What the density means

`dbicop()` and `dvinecop()` return the copula density contribution
defined as the joint density or mass divided by the marginal densities
or masses. This definition is consistent across continuous, discrete,
and mixed models.

For a bivariate model, write
\[u_{j} = F_{j}(x_{j}),\qquad u_{j}^{-} = F_{j}(x_{j}^{-}),\qquad\Delta_{j} = u_{j} - u_{j}^{-}.\]
If both variables are discrete, the copula contribution is the rectangle
probability divided by the two marginal masses:
\[c^{*}(u_{1},u_{2}) = \frac{C(u_{1},u_{2}) - C(u_{1}^{-},u_{2}) - C(u_{1},u_{2}^{-}) + C(u_{1}^{-},u_{2}^{-})}{\Delta_{1}\Delta_{2}}.\]
If only the first variable is discrete, the corresponding mixed
contribution is
\[c^{*}(u_{1},u_{2}) = \frac{\partial_{2}C(u_{1},u_{2}) - \partial_{2}C(u_{1}^{-},u_{2})}{\Delta_{1}}.\]
The continuous density `c(u1, u2)` is recovered when both jump widths
vanish. Vine likelihoods apply these continuous, mixed, or rectangle
contributions edge by edge to the appropriate conditional CDF values.

`dvine()` multiplies that contribution by the fitted marginal densities
or masses and therefore returns the full joint density or probability
mass on the original scale. In general,
\[f_{X}(x) = c^{*}(u,u^{-})\prod\limits_{j = 1}^{d}f_{j}(x_{j}),\] where
each `f_j` denotes a density or a probability mass as appropriate.

## Simulation

`rvinecop()` always returns one uniform-scale value per variable. It
does not invent a marginal distribution for a discrete copula variable.
Use `rvine()` when simulated values should be returned on the original
scale with their integer, ordered, or zero-inflated margins.

``` r
copula_draws <- rvinecop(4, copula_fit)
dim(copula_draws)
#> [1] 4 3

original_draws <- rvine(4, ordered_fit)
is.ordered(original_draws$rating)
#> [1] TRUE
```

## Rosenblatt transforms with atoms

For an atom, the Rosenblatt transform is an interval rather than a
single value. By default, `rosenblatt()` randomizes uniformly within
that interval so the transformed variable is uniform under a correct
model. Set `randomize_discrete = FALSE` to return the deterministic
upper endpoint.

More precisely, if `G_j(· | x_<j)` is the relevant conditional CDF and
`V_j` is an independent standard uniform variable, the randomized
component is
\[Z_{j} = G_{j}(x_{j}^{-} \mid x_{< j}) + V_{j}\{ G_{j}(x_{j} \mid x_{< j}) - G_{j}(x_{j}^{-} \mid x_{< j})\}.\]
The non-randomized version returns the upper endpoint `G_j(x_j | x_<j)`.

``` r
set.seed(402)
randomized <- rosenblatt(compact, copula_fit)
upper_endpoint <- rosenblatt(
  compact,
  copula_fit,
  randomize_discrete = FALSE
)
head(randomized)
#>            [,1]      [,2]      [,3]
#> [1,] 0.05289451 0.0958401 0.7055710
#> [2,] 0.33647905 0.7438199 0.7051919
#> [3,] 0.74612047 0.1121317 0.2203118
#> [4,] 0.13426465 0.5942987 0.1240232
#> [5,] 0.73703700 0.6269100 0.8771097
#> [6,] 0.51260896 0.1971273 0.7066150
head(upper_endpoint)
#>           [,1]      [,2]      [,3]
#> [1,] 0.1353353 0.0958401 0.7851304
#> [2,] 0.4060058 0.7438199 0.7851304
#> [3,] 0.8571235 0.1121317 0.2381033
#> [4,] 0.1353353 0.5942987 0.2381033
#> [5,] 0.8571235 0.6269100 0.8893260
#> [6,] 0.6766764 0.1971273 0.7851304
```

The randomization uses R’s random-number state and is reproducible with
`set.seed()`.

## Missing observations

Copula fitting discards incomplete observations pair by pair, retaining
the maximum information available for each edge. Margin implementations
decide how to handle missing values during marginal fitting; built-in
KDE margins support the package’s established missing-data workflow.
Custom `margin_family()` fitters should either handle missing values or
fail with a useful message.

## Related documentation

  - [Marginal
    models](https://vinecopulib.github.io/rvinecopulib/articles/marginal-models.md)
    documents how `vine()` fits and evaluates integer, ordered, and
    zero-inflated margins.
  - [Conditional simulation and Rosenblatt
    transforms](https://vinecopulib.github.io/rvinecopulib/articles/conditional-simulation.md)
    develops the multivariate transform and conditional-sampling
    workflow.
  - The [`vine()` distribution
    reference](https://vinecopulib.github.io/rvinecopulib/reference/vine_methods.md)
    and [vine-copula distribution
    reference](https://vinecopulib.github.io/rvinecopulib/reference/vinecop_methods.md)
    specify the accepted evaluation layouts.
