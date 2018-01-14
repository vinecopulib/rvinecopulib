rvinecopulib
==========

[![Build status Linux](https://travis-ci.org/vinecopulib/rvinecopulib.svg?branch=master)](https://travis-ci.org/vinecopulib/rvinecopulib)
[![Windows Build status](http://ci.appveyor.com/api/projects/status/github/vinecopulib/rvinecopulib?svg=true)](https://ci.appveyor.com/project/vinecopulib/rvinecopulib)
[![Coverage Status](https://img.shields.io/codecov/c/github/vinecopulib/rvinecopulib/master.svg)](https://codecov.io/github/vinecopulib/rvinecopulib?branch=master)
[![CRAN version](http://www.r-pkg.org/badges/version/rvinecopulib)](https://cran.r-project.org/package=rvinecopulib) 
[![CRAN downloads](http://cranlogs.r-pkg.org/badges/rvinecopulib)](https://cran.r-project.org/package=rvinecopulib)

Vine copulas are a flexible class of dependence models consisting of bivariate 
building blocks (see e.g., Aas et al., 2009). You can find a comprehensive 
list of publications and other materials on [vine-copula.org](http://www.statistics.ma.tum.de/en/research/vine-copula-models/).

This package is the [R](https://cran.r-project.org/) API to the C++ library 
[vinecopulib](https://github.com/vinecopulib/vinecopulib), a header-only 
C++ library for vine copula models based on [Boost](http://www.boost.org/) and 
[Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page).

It provides high-performance implementations of the core features of the popular
[VineCopula R library](https://github.com/tnagler/VineCopula), in particular
inference algorithms for both vine copula and bivariate copula models.
Advantages over VineCopula are  
* a sleaker and more modern API,
* shorter runtimes, especially in high dimensions,
* nonparametric and multi-parameter families.

As VineCopula, the package is primarily made for the statistical analysis of 
**vine copula models**. The package includes tools for parameter estimation, 
model selection, simulation, and visualization. Tools for estimation, selection 
and exploratory data analysis of **bivariate copula** models are also provided. 
Please see the [API documentation](https://vinecopulib.github.io/rvinecopulib/) 
for a detailed description of all functions.

Table of contents
-----------------

- [How to install](#how-to-install)
- [Package overview](#package-overview)
	- [Bivariate copula modeling: bicop_dist and bicop](#bivariate-copula-modeling-bicop_dist-and-bicop)
	- [Vine copula modeling: vinecop_dist and vinecop](#vine-copula-modeling-vinecop_dist-and-vinecop)
	- [Bivariate copula families](#bivariate-copula-families)
- [References](#references)

------------------------------------------------------------------------


How to install
--------------


You can install:

-   the stable release on CRAN:

    ``` r
    install.packages("rvinecopulib")
    ```

-   the latest development version:

    ``` r
    devtools::install_github("vinecopulib/rvinecopulib")
    ```

------------------------------------------------------------------------

Package overview
----------------

Below, we list most functions and features you should know about. As usual in 
copula models, data are assumed to be serially independent and lie in the unit
hypercube. 

### Bivariate copula modeling: bicop_dist and bicop

  * `bicop_dist`: Creates a bivariate copula by specifying the family, rotation 
    and parameters. Returns an object of class `bicop_dist`. The class has the
    following methods:
     
     * `print`: a brief overview of the bivariate copula. 
            
     * `plot`, `contour`: surface/perspective and contour plots of the copula
        density. Possibly coupled with standard normal margins (default for
        `contour`). 
        
  * `dbicop`, `pbicop`, `rbicop`, `hbicop`: Density, distribution function, random 
    generation and H-functions (with their inverses) for bivariate copula 
    distributions. Additionally to the evaluation points, you can provide 
    either `family`, `rotation` and `parameter`, or an object of class 
    `bicop_dist`.

  * `bicop`: Estimates parameters of a bivariate copula. Estimation can be done 
    by maximum likelihood (`par_method = "mle"`) or inversion of the empirical 
    Kendall's tau (`par_method = "itau"`, only available for one-parameter 
    families) for parametric families, and using local-likelihood 
    approximations of order zero/one/two for nonparametric models 
    (`nonpar_method="constant"`/`nonpar_method="linear"`/`nonpar_method="quadratic"`). 
    If `family_set` is a vector of families, then the family is selected using
    `selcrit="loglik"`, `selcrit="aic"` or `selcrit="bic"`. The function 
    returns an object of classes `bicop` and `bicop_dist`.
    The class `bicop` has the following following methods:
    
     * `print`: a more comprehensive overview of the bivariate copula model 
       with fit statistics. 
            
     * `predict`, `fitted`: predictions and fitted values for a bivariate 
       copula model.
       
     * `nobs`, `logLik`, `AIC`, `BIC`: usual fit statistics.

### Vine copula modeling: vinecop_dist and vinecop

  * `vinecop_dist`: Creates a vine copula by specifying a nested list of 
    `bicop_dist` objects and a quadratic structure matrix. 
    Returns an object of class `vinecop_dist`. The class has the
    following methods:
     
     * `print`, `summary`: a brief and more comprehensive overview of the vine 
       copula. 
            
     * `plot`: plots of the vine structure. 
        
  * `dvinecop`, `pvinecop`, `rvinecop`: Density, distribution function, random 
    generation for vine copula distributions. 

  * `vinecop`: automated fitting for vine copula models. The function inherits 
    the parameters of `bicop`. Optionally, a quadratic `matrix` can be used as 
    input to pre-specify the vine structure. `tree_crit` describes the 
    criterion for tree selection, one of `"tau"`, `"rho"`, `"hoeffd"` for 
    Kendall's tau, Spearman's rho, and Hoeffding's D, respectively.
    Additionally, `threshold` allows to threshold the `tree_crit` and 
    `trunc_lvl` to truncate the vine copula, with `threshold_sel` and 
    `trunc_lvl_sel` to automatically select both parameters. The function 
    returns an object of classes `vinecop` and `vinecop_dist`.
    The class has the `vinecop` has the following following methods:
    
     * `print`, `summary`: a brief and more comprehensive overview of the vine 
       copula with additional fit statistics information.
            
     * `predict`, `fitted`: predictions and fitted values for a vine 
       copula model.
       
     * `nobs`, `logLik`, `AIC`, `BIC`: usual fit statistics.

### Bivariate copula families

In this package several bivariate copula families are included for bivariate 
and multivariate analysis using vine copulas. It provides 
functionality of elliptical (Gaussian and Student-t) as well as Archimedean 
(Clayton, Gumbel, Frank, Joe, BB1, BB6, BB7 and BB8) copulas to cover a large
range of dependence patterns. For Archimedean copula families,
rotated versions are included to cover negative dependence as well. 
Additionally, nonparametric families are also supported.

| type          | name                  | name in R     |
|---------------|-----------------------|---------------|
| -             | Independence          | "indep"       |
| Elliptical    | Gaussian              | "gaussian"    |
| "             | Student t             | "student"     |
| Archimedean   | Clayton               | "clayton"     |
| "             | Gumbel                | "gumbel"      |
| "             | Frank                 | "frank"       |
| "             | Joe                   | "joe"         |
| "             | Clayton-Gumbel (BB1)  | "bb1"         |
| "             | Joe-Gumbel (BB6)      | "bb6"         |
| "             | Joe-Clayton (BB7)     | "bb7"         |
| "             | Joe-Frank (BB8)       | "bb8"         |
| Nonparametric | Transformation kernel | "tll"         |

Note that several convenience vectors of families are included:
* `"all"` contains all the families
* `"parametric"` contains the parametric families (all except `"tll"`)
* `"nonparametric"` contains the nonparametric families (`"indep"` and `"tll"`)
* `"one_par"` contains the parametric families with a single parameter
(`"gaussian"`, `"clayton"`, `"gumbel"`, `"frank"`, and `"joe"`)
* `"two_par"` contains the parametric families with two parameters
(`"student"`, `"bb1"`, `"bb6"`, `"bb7"`, and `"bb8"`)
* `"elliptical"` contains the elliptical families
* `"archimedean"` contains the archimedean families
* `"BB"` contains the BB families
* `"itau"` families for which estimation by Kendall's tau inversion is available
(`"indep"`,`"gaussian"`, `"student"`,`"clayton"`, `"gumbel"`, `"frank"`, `"joe"`)

The following table shows the parameter ranges of bivariate copula families with 
one or two parameters:

| Copula family                    | `par[1]`        | `par[2]`  |
|:---------------------------------|:-------------|:-------------|
| Gaussian                         | `(-1, 1)`    | -            |
| Student t                        | `(-1, 1)`    | `(2,Inf)`    |
| Clayton                          | `(0, Inf)`   | -            |
| Gumbel                           | `[1, Inf)`   | -            |
| Frank                            | `R \ {0}`    | -            |
| Joe                              | `(1, Inf)`   | -            |
| Clayton-Gumbel (BB1)             | `(0, Inf)`   | `[1, Inf)`   |
| Joe-Gumbel (BB6)                 | `[1 ,Inf)`   | `[1, Inf)`   |
| Joe-Clayton (BB7)                | `[1, Inf)`   | `(0, Inf)`   |
| Joe-Frank (BB8)                  | `[1, Inf)`   | `(0, 1]`     |

------------------------------------------------------------------------

References
----------

Aas, K., C. Czado, A. Frigessi, and H. Bakken (2009). Pair-copula constructions of multiple dependence. Insurance: Mathematics and Economics 44 (2), 182-198.
