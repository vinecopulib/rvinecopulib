# Predictions and fitted values for a vine copula model

Predictions of the density and distribution function for a vine copula
model.

## Usage

``` r
# S3 method for class 'vinecop'
predict(object, newdata, what = "pdf", n_mc = 10^4, cores = 1, ...)

# S3 method for class 'vinecop'
fitted(object, what = "pdf", n_mc = 10^4, cores = 1, ...)
```

## Arguments

- object:

  a `vinecop` object.

- newdata:

  points where the fit shall be evaluated.

- what:

  what to predict, either `"pdf"` or `"cdf"`.

- n_mc:

  number of samples used for quasi Monte Carlo integration when
  `what = "cdf"`.

- cores:

  number of cores to use; if larger than one, computations are done in
  parallel on `cores` batches.

- ...:

  unused.

## Value

[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
[`predict()`](https://rdrr.io/r/stats/predict.html) have return values
similar to
[`dvinecop()`](https://vinecopulib.github.io/rvinecopulib/reference/vinecop_methods.md)
and
[`pvinecop()`](https://vinecopulib.github.io/rvinecopulib/reference/vinecop_methods.md).

## Details

[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) can only be
called if the model was fit with the `keep_data = TRUE` option.

### Discrete variables

When at least one variable is discrete, two types of "observations" are
required in `newdata`: the first \\n \\ x \\ d\\ block contains
realizations of \\F\_{X_j}(X_j)\\. The second \\n \\ x \\ d\\ block
contains realizations of \\F\_{X_j}(X_j^-)\\. The minus indicates a
left-sided limit of the cdf. For, e.g., an integer-valued variable, it
holds \\F\_{X_j}(X_j^-) = F\_{X_j}(X_j - 1)\\. For continuous variables
the left limit and the cdf itself coincide. Respective columns can be
omitted in the second block.

## Examples

``` r
u <- sapply(1:5, function(i) runif(50))
fit <- vinecop(u, family_set = "par", keep_data = TRUE)
all.equal(predict(fit, u), fitted(fit), check.environment = FALSE)
#> [1] TRUE
```
