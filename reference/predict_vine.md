# Predictions and fitted values for a vine copula model

Predictions of the density and distribution function for a vine copula
model.

## Usage

``` r
# S3 method for class 'vine'
predict(object, newdata, what = "pdf", n_mc = 10^4, cores = 1, ...)

# S3 method for class 'vine'
fitted(object, what = "pdf", n_mc = 10^4, cores = 1, ...)
```

## Arguments

- object:

  a `vine` object.

- newdata:

  points where the fit shall be evaluated.

- what:

  what to predict, either `"pdf"` or `"cdf"`.

- n_mc:

  number of samples used for quasi Monte Carlo integration when
  `what = "cdf"`.

- cores:

  number of cores to use; if larger than one, computations are done in
  parallel on `cores` batches .

- ...:

  unused.

## Value

[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
[`predict()`](https://rdrr.io/r/stats/predict.html) have return values
similar to
[`dvine()`](https://vinecopulib.github.io/rvinecopulib/reference/vine_methods.md)
and
[`pvine()`](https://vinecopulib.github.io/rvinecopulib/reference/vine_methods.md).

## Examples

``` r
x <- sapply(1:5, function(i) rnorm(50))
fit <- vine(x, copula_controls = list(family_set = "par"), keep_data = TRUE)
all.equal(predict(fit, x), fitted(fit), check.environment = FALSE)
#> [1] TRUE
```
