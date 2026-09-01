# Modified vine copula Bayesian information criterion (mBICv)

Calculates the modified vine copula Bayesian information criterion.

## Usage

``` r
mBICV(object, psi0 = 0.9, newdata = NULL)
```

## Arguments

  - object:
    
    a fitted `vinecop` object.

  - psi0:
    
    baseline prior probability of a non-independence copula.

  - newdata:
    
    optional; a new data set.

## Details

The modified vine copula Bayesian information criterion (mBICv) is
defined as

$$BIC = -2 loglik + \\nu log(n) - 2 \\sum\_{t=1}^{d - 1} (q\_t
log(\\psi\_0^t) + (d - t - q\_t) log(1 - \\psi\_0^t)) $$

where \\(\\mathrm{loglik}\\) is the log-likelihood and \\(\\nu\\) is the
(effective) number of parameters of the model, \\(t\\) is the tree level
\\(\\psi\_0\\) is the prior probability of having a non-independence
copula and \\(q\_t\\) is the number of non-independence copulas in tree
\\(t\\). The mBICv is a consistent model selection criterion for
parametric sparse vine copula models.

## References

Nagler, T., Bumann, C., Czado, C. (2019). Model selection for sparse
high-dimensional vine copulas with application to portfolio risk.
*Journal of Multivariate Analysis, in press*
(<https://arxiv.org/pdf/1801.09739>)

## Examples

``` r
u <- sapply(1:5, function(i) runif(50))
fit <- vinecop(u, family_set = "par", keep_data = TRUE)
mBICV(fit, 0.9) # with a 0.9 prior probability of a non-independence copula
#> [1] 30.65039
mBICV(fit, 0.1) # with a 0.1 prior probability of a non-independence copula
#> [1] 32.48713
```
