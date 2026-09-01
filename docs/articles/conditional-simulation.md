# Conditional simulation and Rosenblatt transforms

A vine is simulated along an order. Variables to condition on must form
the tail of an admissible sampling order. rvinecopulib can enforce this
during structure selection or transiently reorient a compatible fitted
model during simulation.

## Select a conditioning-aware structure

Pass variable indices or names through `conditioning_set` when fitting a
copula. The selected structure places those variables at the end of its
order.

``` r
n <- 180
z <- rnorm(n)
x <- cbind(
  response1 = z + rnorm(n),
  response2 = -0.5 * z + rnorm(n),
  driver1 = z + rnorm(n, sd = 0.5),
  driver2 = 0.4 * z + rnorm(n)
)
u <- pseudo_obs(x)

fit <- vinecop(
  u,
  family_set = "onepar",
  conditioning_set = c("driver1", "driver2")
)
get_structure(fit)
#> 4-dimensional R-vine structure ('rvine_structure')
#> 3 3 3 3
#> 1 4 4  
#> 4 1    
#> 2
```

Conditioning-aware selection supports fixed and automatically selected
truncation levels. It requires an MST tree algorithm (`"mst_prim"` or
`"mst_kruskal"`); random tree algorithms do not guarantee the required
order.

## Simulate conditionally on the copula scale

`u_cond` contains one value per conditioning variable. A vector or
one-row matrix is repeated for all simulations; an `n`-row matrix
specifies different conditions for every output row.

``` r
condition <- c(0.25, 0.8)
draws <- rvinecop(
  100,
  fit,
  u_cond = condition,
  conditioning_set = c("driver1", "driver2")
)
stopifnot(
  isTRUE(all.equal(draws[, "driver1"], rep(condition[1], 100))),
  isTRUE(all.equal(draws[, "driver2"], rep(condition[2], 100)))
)
head(draws)
#>       response1 response2 driver1 driver2
#> [1,] 0.63870946 0.4514150    0.25     0.8
#> [2,] 0.31615622 0.3556760    0.25     0.8
#> [3,] 0.07066392 0.3225717    0.25     0.8
#> [4,] 0.54276866 0.8714729    0.25     0.8
#> [5,] 0.87804170 0.9514778    0.25     0.8
#> [6,] 0.40480490 0.9892094    0.25     0.8
```

When `conditioning_set` is omitted, the columns of `u_cond` refer to the
last variables in the model’s current order.

``` r
conditions <- cbind(
  driver1 = seq(0.1, 0.9, length.out = 5),
  driver2 = rep(0.5, 5)
)
row_draws <- rvinecop(
  5,
  fit,
  u_cond = conditions,
  conditioning_set = c("driver1", "driver2")
)
stopifnot(isTRUE(all.equal(row_draws[, 3:4], conditions)))
```

## Condition on the original data scale

`vine()` exposes the same structure control through `copula_controls`.
`rvine()` accepts `x_cond` on the original scale and uses the fitted
margins to compute all necessary CDF values and left limits.

``` r
full_fit <- vine(
  as.data.frame(x),
  margins_controls = list(family_set = "kde1d"),
  copula_controls = list(
    family_set = "onepar",
    conditioning_set = c("driver1", "driver2")
  )
)

original_condition <- data.frame(driver1 = -0.5, driver2 = 1.25)
original_draws <- rvine(
  6,
  full_fit,
  x_cond = original_condition,
  conditioning_set = c("driver1", "driver2")
)
stopifnot(
  isTRUE(all.equal(original_draws[, "driver1"], rep(-0.5, 6))),
  isTRUE(all.equal(original_draws[, "driver2"], rep(1.25, 6)))
)
```

Column names are strongly recommended: they make the mapping explicit
and protect code from changes in variable order.

## Reorient an existing model transiently

Supplying an explicit `conditioning_set` to `rvinecop()` or `rvine()`
asks the backend to evaluate an equivalent orientation without modifying
the model. Only conditioning sets that can form an admissible tail of
the fitted structure are possible. If a specific set is central to the
analysis, request it during fitting; selection then guarantees a usable
structure.

## Discrete conditioning variables

On the copula scale, a discrete conditioning variable needs its value
`F(x)` and left limit `F(x-)`. As with model evaluation, `u_cond`
accepts:

  - a compact layout with one value column per conditioning variable
    followed by one left-limit column per discrete variable; or
  - an expanded layout with all value columns followed by all left-limit
    columns, where a continuous variable’s two columns must be equal.

`rvine()` handles this internally from `x_cond`, including ordered and
zero-inflated margins.

## Rosenblatt transforms

The Rosenblatt transform maps observations from a fitted distribution to
independent uniforms. For continuous variables in a sampling sequence
`s1, ..., sd`, its components are
\[Z_{j} = F_{s_{j} \mid s_{1},\ldots,s_{j - 1}}\left( X_{s_{j}} \mid X_{s_{1}},\ldots,X_{s_{j - 1}} \right),\qquad j = 1,\ldots,d,\]
where the first component is unconditional. If the fitted model is
correct, the `Z_j` are independent standard uniforms. The inverse
transform applies the corresponding conditional quantiles recursively
and therefore maps independent uniforms back to the model.

``` r
simulated <- rvinecop(100, fit)
independent <- rosenblatt(simulated, fit)
reconstructed <- inverse_rosenblatt(independent, fit)
max(abs(simulated - reconstructed))
#> [1] 7.21645e-16
```

The transform follows the model’s sampling order. An explicit
`conditioning_set` requests an admissible order ending in those
variables, which is useful when conditional quantiles must follow the
same orientation as conditional simulation.

``` r
transformed <- rosenblatt(
  simulated,
  fit,
  conditioning_set = c("driver1", "driver2")
)
inverse_rosenblatt(
  transformed,
  fit,
  conditioning_set = c("driver1", "driver2")
)[1:3, ]
#>      response1 response2    driver1   driver2
#> [1,] 0.1180315 0.6240082 0.02953261 0.0577602
#> [2,] 0.8016499 0.6850150 0.30646778 0.7316434
#> [3,] 0.1138517 0.4059435 0.29654106 0.8204806
```

For discrete variables, the transform randomizes within the jump
interval by default. See the [discrete-data
guide](https://vinecopulib.github.io/rvinecopulib/articles/discrete-data.md)
for the randomized formula, deterministic upper endpoints, and required
copula layouts.

## Reproducibility

Conditional simulation and randomized discrete transforms use R’s
random-number state. Call `set.seed()` immediately before the operation
when a reproducible draw is required. The original model object is never
mutated by a conditional simulation or transform.

## Related documentation

  - [Vine copula models and
    structures](https://vinecopulib.github.io/rvinecopulib/articles/vine-copula-models.md)
    explains the structure order and selection controls used here.
  - [Marginal
    models](https://vinecopulib.github.io/rvinecopulib/articles/marginal-models.md)
    documents the transformation between original observations and the
    copula scale.
  - The [`rvinecop()`
    reference](https://vinecopulib.github.io/rvinecopulib/reference/vinecop_methods.md),
    [`rvine()`
    reference](https://vinecopulib.github.io/rvinecopulib/reference/vine_methods.md),
    and [Rosenblatt
    reference](https://vinecopulib.github.io/rvinecopulib/reference/rosenblatt.md)
    give complete argument and return value details.
