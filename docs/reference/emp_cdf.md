# Corrected Empirical CDF

The empirical CDF with tail correction, ensuring that its output is
never 0 or 1.

## Usage

``` r
emp_cdf(x)
```

## Arguments

  - x:
    
    numeric vector of observations

## Value

A function with signature `function(x)` that returns \\(F\_n(x)\\).

## Details

The corrected empirical CDF is defined as $$ F\_n(x) = \\frac{1}{n + 1}
\\max\\biggl\\{1, \\sum\_{i = 1}^n 1(X\_i \\le x)\\biggr\\} $$

## Examples

``` r
# fit ECDF on simulated data
x <- rnorm(100)
cdf <- emp_cdf(x)

# output is bounded away from 0 and 1
cdf(-50)
#> [1] 0.00990099
cdf(50)
#> [1] 0.990099
```
