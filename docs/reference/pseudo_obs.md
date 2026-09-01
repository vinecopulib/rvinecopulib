# Pseudo-Observations

Compute the pseudo-observations for the given data matrix.

## Usage

``` r
pseudo_obs(x, ties_method = "average", lower_tail = TRUE)
```

## Arguments

  - x:
    
    vector or matrix random variates to be converted (column wise) to
    pseudo-observations.

  - ties\_method:
    
    similar to `ties.method` of `rank()` (only `"average"`, `"first"`
    and `"random"` currently available).

  - lower\_tail:
    
    `logical` which, if \`FALSE“, returns the pseudo-observations when
    applying the empirical marginal survival functions.

## Value

a vector of matrix of the same dimension as the input containing the
pseudo-observations.

## Details

Given `n` realizations \\(x\_i=(x\_{i1}, \\ldots,x\_{id})\\), \\(i \\in
\\left\\lbrace 1, \\ldots,n \\right\\rbrace \\) of a random vector `X`,
the pseudo-observations are defined via \\(u\_{ij}=r\_{ij}/(n+1)\\) for
\\(i \\in \\left\\lbrace 1, \\ldots,n \\right\\rbrace\\) and \\(j \\in
\\left\\lbrace 1, \\ldots,d \\right\\rbrace \\), where \\(r\_{ij}\\)
denotes the rank of \\(x\_{ij}\\) among all \\(x\_{kj}\\), \\(k \\in
\\left\\lbrace 1, \\ldots,n \\right\\rbrace \\).

The pseudo-observations can thus also be computed by component-wise
applying the empirical distribution functions to the data and scaling
the result by \\(n/(n+1)\\). This asymptotically negligible scaling
factor is used to force the variates to fall inside the open unit
hypercube, for example, to avoid problems with density evaluation at the
boundaries.

When `lower_tail = FALSE`, then `pseudo_obs()` simply returns `1 -
pseudo_obs()`.

## Examples

``` r
# pseudo-observations for a vector
pseudo_obs(rnorm(10))
#>  [1] 0.90909091 0.45454545 0.81818182 0.09090909 0.63636364 0.72727273
#>  [7] 0.27272727 0.18181818 0.54545455 0.36363636

# pseudo-observations for a matrix
pseudo_obs(cbind(rnorm(10), rnorm(10)))
#>             [,1]       [,2]
#>  [1,] 0.63636364 0.63636364
#>  [2,] 0.72727273 0.81818182
#>  [3,] 0.27272727 0.90909091
#>  [4,] 0.45454545 0.09090909
#>  [5,] 0.09090909 0.36363636
#>  [6,] 0.18181818 0.27272727
#>  [7,] 0.54545455 0.45454545
#>  [8,] 0.90909091 0.72727273
#>  [9,] 0.36363636 0.18181818
#> [10,] 0.81818182 0.54545455
```
