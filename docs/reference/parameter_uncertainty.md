# Parameter uncertainty

Covariance matrices and Wald confidence intervals for the estimated
pair-copula parameters.

## Usage

``` r
# S3 method for class 'vinecop_dist'
vcov(object, newdata = NULL, step_wise = TRUE, cores = 1, ...)

# S3 method for class 'bicop_dist'
vcov(object, newdata = NULL, step_wise = TRUE, cores = 1, ...)

# S3 method for class 'vinecop_dist'
confint(object, parm, level = 0.95, ...)

# S3 method for class 'bicop_dist'
confint(object, parm, level = 0.95, ...)
```

## Arguments

  - object:
    
    an object of class `"bicop"` or `"vinecop"`, or of class
    `"bicop_dist"` or `"vinecop_dist"` together with `newdata`.

  - newdata:
    
    copula data to evaluate the derivatives at. Defaults to the data
    stored in `object`, which requires `keep_data = TRUE` at fitting
    time.

  - step\_wise:
    
    if `TRUE` (default), derivatives of the step-wise estimating
    equations; if `FALSE`, of the joint log-likelihood. See `scores()`.

  - cores:
    
    number of cores to use; if larger than one, computations are done in
    parallel on `cores` batches.

  - ...:
    
    passed on to `vcov()`.

  - parm:
    
    a specification of which parameters to give intervals for: a vector
    of names or of positions. Defaults to all of them.

  - level:
    
    the confidence level.

## Value

`vcov()` returns a covariance matrix with one row and column per
estimated pair-copula parameter, ordered by tree, then edge, then
parameter within the pair copula, with dimnames giving that position.
`confint()` returns a matrix with columns for the lower and upper
endpoints.

## Details

The estimator solves \\(\\bar s(\\hat\\theta) = 0\\) with \\(\\bar
s(\\theta) = n^{-1} \\sum\_i s\_i(\\theta)\\), so it is an M-estimator.
Writing $$A = -\\partial \\bar s(\\theta) / \\partial \\theta^\\top,
\\qquad B = n^{-1} \\sum\_i s\_i s\_i^\\top,$$ `vcov()` returns the
sandwich \\(A^{-1} B A^{-\\top} / n\\) of White (1982).

### Which objective is differentiated

`step_wise` selects the estimating equation, and it changes the nature
of \\(A\\).

With `step_wise = FALSE`, \\(\\bar s\\) is the gradient of the joint
log-likelihood and \\(A\\) is minus its averaged Hessian: a symmetric
matrix, and the classical setting for both formulas.

With `step_wise = TRUE` (the default), the pseudo-observations entering
each pair copula are held fixed. This is the estimating equation that
sequential, tree-by-tree estimation actually solves, but it is not the
gradient of any objective: \\(A\\) is a Jacobian rather than a Hessian,
and it is block triangular rather than symmetric, which is why the
sandwich is written with \\(A^{-\\top}\\) on the right.

### Restrictions

No derivatives are implemented for nonparametric pair copulas, so a
model containing one is refused. For models with discrete variables the
backend differentiates by central finite differences rather than in
closed form; the result is correct but carries the accuracy of a
finite-difference approximation.

### Why there is no method for vine distributions

These functions treat the copula data as observed, so they quantify
uncertainty in the pair-copula parameters *given* the margins. A model
fitted by `vine()` estimates its margins too, and the second stage
inherits the first stage's sampling variability: intervals that ignore
it are too short. In a small experiment (three-dimensional Gaussian
D-vine, equicorrelation 0.5, `n = 1000`, 1000 replications) nominal 95%
intervals attained 0.956 when the true margins were used and 0.929 when
they were estimated by ranks.

Correcting for it needs the influence function of the marginal
estimator, and `fit_margin()` is a user-supplied black box, so there is
nothing to differentiate. We therefore provide no method for `"vine"`
objects rather than a silently anticonservative one. To quantify
uncertainty for a complete vine distribution, resample:

    boot <- replicate(B, {
      idx <- sample(nrow(x), replace = TRUE)
      fit <- vine(x[idx, ], margins_controls = mc, copula_controls = cc)
      unlist(get_all_parameters(fit$copula))
    })
    apply(boot, 1, quantile, c(0.025, 0.975))

## References

White H (1982). "Maximum Likelihood Estimation of Misspecified Models."
*Econometrica*, **50**(1), 1–25.
[doi:10.2307/1912526](https://doi.org/10.2307/1912526)

Stoeber J, Schepsmeier U (2013). "Estimating Standard Errors in Regular
Vine Copula Models." *Computational Statistics*, **28**(6), 2679–2707.
[doi:10.1007/s00180-013-0423-8](https://doi.org/10.1007/s00180-013-0423-8)

## See also

`scores()`, `hessian()`, `vinecop()`, `vine()`

## Examples

``` r
# bivariate copula ------------------------------------------
u <- rbicop(500, "gaussian", 0, 0.5)
fit <- bicop(u, family_set = "parametric", keep_data = TRUE)
sqrt(diag(vcov(fit)))
#>   bb1.par1   bb1.par2 
#> 0.07111260 0.06084133 
confint(fit)
#>               2.5 %    97.5 %
#> bb1.par1 0.04609626 0.3248526
#> bb1.par2 1.16197193 1.4004656

# vine copula -----------------------------------------------
u <- pseudo_obs(matrix(rnorm(500 * 3), 500, 3))
fit <- vinecop(u, family_set = "parametric", keep_data = TRUE)
vcov(fit)
#>             T2.2,1;3
#> T2.2,1;3 0.002138162
confint(fit, level = 0.9)
#>                  5 %      95 %
#> T2.2,1;3 0.002526306 0.1546433

# derivatives of the joint log-likelihood instead of the step-wise ones
sqrt(diag(vcov(fit, step_wise = FALSE)))
#>   T2.2,1;3 
#> 0.04624027 
```
