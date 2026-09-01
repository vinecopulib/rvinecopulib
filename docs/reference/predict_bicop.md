# Predictions and fitted values for a bivariate copula model

Predictions of the density, distribution function, h-functions (with
their inverses) for a bivariate copula model.

## Usage

``` r
# S3 method for class 'bicop_dist'
predict(object, newdata, what = "pdf", ...)

# S3 method for class 'bicop'
fitted(object, what = "pdf", ...)
```

## Arguments

  - object:
    
    a `bicop` object.

  - newdata:
    
    points where the fit shall be evaluated.

  - what:
    
    what to predict, one of `"pdf"`, `"cdf"`, `"hfunc1"`, `"hfunc2"`,
    `"hinv1"`, `"hinv2"`.

  - ...:
    
    unused.

## Value

`fitted()` and `logLik()` have return values similar to `dbicop()`,
`pbicop()`, and `hbicop()`.

## Details

`fitted()` can only be called if the model was fit with the `keep_data =
TRUE` option.

### Discrete variables

When at least one variable is discrete, more than two columns are
required for `newdata`: the first \\(n \\times 2\\) block contains
realizations of \\(F\_{X\_1}(x\_1), F\_{X\_2}(x\_2)\\). The second \\(n
\\times 2\\) block contains realizations of \\(F\_{X\_1}(x\_1^-),
F\_{X\_2}(x\_2^-)\\). The minus indicates a left-sided limit of the cdf.
For, e.g., an integer-valued variable, it holds \\(F\_{X\_1}(x\_1^-) =
F\_{X\_1}(x\_1 - 1)\\). For continuous variables the left limit and the
cdf itself coincide. Respective columns can be omitted in the second
block.

## Examples

``` r
# Simulate and fit a bivariate copula model
u <- rbicop(500, "gauss", 0, 0.5)
fit <- bicop(u, family_set = "par", keep_data = TRUE)

# Predictions
all.equal(predict(fit, u, "hfunc1"), fitted(fit, "hfunc1"),
          check.environment = FALSE)
#> [1] TRUE
```
