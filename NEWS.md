# rvinecopulib 0.7.1.1.2

BUG FIX

* Fixes "deprecated-literal-operator" warning on clang20.

# rvinecopulib 0.7.1.1.1

BUG FIX

* fix handling of discrete variables in `vine()` models and related functions.

# rvinecopulib 0.7.1.1.0

Update following a upgrade of the C++ backend vinecopulib to 0.7.1, see
https://github.com/vinecopulib/vinecopulib/blob/main/NEWS.md.

The main changes on the R end are:

* improved documentation,

* support for zero-inflated variables,

* added new Tawn copula family,

* new argument `allow_rotations` to disable rotations of copula families,

* added variable names to vinecop summary (#276)

* fixed handling of logistic distribution (#275)

* fix NA handling in vine() control checks (#266)

* allow bicop_dist() with tll (#268)

# rvinecopulib 0.6.3.1.1

- add `-D_HAS_AUTO_PTR_ETC=0` flag to disable deprecated features used in boost.

# rvinecopulib 0.6.3.1.0

- fix `NA` handling in `to_pseudo_obs()` (#260) 

- add `emp_cdf()` for the tail corrected empirical cdf (#261) 

# rvinecopulib 0.6.2.1.3 (December 3, 2022)

- fix marginal PIT for discrete variables (see issue #257, thanks @rplzzz)

# rvinecopulib 0.6.2.1.2 (October 16, 2022)

- fix warning about C++17 attribute extension 'nodiscard' (#255)

# rvinecopulib 0.6.2.1.1 (August 30, 2022)

- replace bitwise operations on Boolean variables.

# rvinecopulib 0.6.2.1.0 (August 26, 2022)

Release following the updates of vinecopulib to 0.6.2, see 
https://github.com/vinecopulib/vinecopulib/releases. The most 
relevant changes are summarized below.

* improved documentation (discrete and missing dat, Rosenblatt transforms)

* better parallelization when there is a small number of edges (#555)


# rvinecopulib 0.6.1.1.3 (March 17, 2022)

- avoid bitwise operations on Boolean variables.


# rvinecopulib 0.6.1.1.2 (March 14, 2022)

- enforce use on recent version of RcppThread for proper handling of linker flags


# rvinecopulib 0.6.1.1.1 (October 6, 2021)

### BUG FIXES

- remove illegal pragmas from json header (#245) 

- allow tree restriction in summary functions (#244)


# rvinecopulib 0.6.1.1.0 (July 13, 2021)

Release following the updates of vinecopulib to 0.6.1, see 
https://github.com/vinecopulib/vinecopulib/releases. The most 
relevant changes are summarized below.

### BUG FIXES

- improved documentation (#241, #239)

- use `num_threads` in recursive calls to the inverse Rosenblatt

- force TLL to be nonnegative (#238)

- fix number of parameters for TLL


# rvinecopulib 0.5.5.1.1 (December 15, 2020)

Maintenance release following the changes to `all.equal()` in R 4.1.x.


# rvinecopulib 0.5.5.1.0 (November 24, 2020)

Release following the updates of vinecopulib to 0.5.5, see 
https://github.com/vinecopulib/vinecopulib/releases. The most 
relevant changes are summarized below.

### BUG FIXES

  * fix little bug in copula selection based on mBIC

  * stabilize BB7 copula pdf

  * fix threshold selection for (near-)independent data

  * fix vine copula selection for 1-dimensional models with discrete variables

  * fix user-visible variable types


# rvinecopulib 0.5.4.1.0 (September 30, 2020)

Release following the updates of vinecopulib to 0.5.4, see 
https://github.com/vinecopulib/vinecopulib/releases. The most 
relevant changes are summarized below.

### BUG FIXES

  * fix uninitialized number of parameters for TLL family

  * fix Kendall's tau of Frank copula for par <= 3

  * fix `dvinecop()` when discrete variables are present (#222)


# rvinecopulib 0.5.3.1.0 (August 11, 2020)

Release following the updates of vinecopulib to 0.5.3, see 
https://github.com/vinecopulib/vinecopulib/releases. The most 
relevant changes are summarized below.

### NEW FEATURES

  * allow 1-dimensional models (#216) 
  
  * make AIC default selection criterion (#213)

### BUG FIXES

  * catch na in ktau_to_par (#214)
  
  * make Bicop/Vinecop objects indepent of copied-from-objects

  * enforce parameters bounds in tau_to_parameters for Archimedean families


# rvinecopulib 0.5.2.1.0 (May 7, 2020)

Release following the updates of vinecopulib to 0.5.2, see 
https://github.com/vinecopulib/vinecopulib/releases. The most 
relevant changes are summarized below.

NEW FEATURES

  * single-integer constructors for `dvine_structure()`/`cvine_structure()`. (#203)

  * add `var_names = "hide"` option in `plot.vinecop_dist()`. (#203)

  * add function `plot.rvine_matrix()`. (#203)

BUG FIXES

  * fix bug for (negative) tau to parameter conversion for Frank family. (#207)

  * fix rare error `rvine_structure_sim()/rvine_matrix_sim()` 

  * safeguard `"tll"` family against comonotonic data. 

  * stabilize archimedean h-functions near independence.


# rvinecopulib 0.5.1.1.0 (November 25, 2019)

Release following the updates of vinecopulib to 0.5.1, see 
https://github.com/vinecopulib/vinecopulib/releases. 

BUG FIX

  * fix out of range bug for weighted TLL influence when sample size is small.  


# rvinecopulib 0.5.0.1.0 (November 25, 2019)

Release following the updates of vinecopulib to 0.5.0, see 
https://github.com/vinecopulib/vinecopulib/releases. The most 
relevant changes are summarized below.

NEW FEATURES

  * modelling discrete variables with bivariate or vine copulas. (#195)

  * selection of partially specified R-vine structures. (#195)

  * convenience classes `dvine_structure()`/`cvine_structure()` for D- and 
    C-vine structures. (#195)

  * new criterion for tree selection: `"joe"` corresponds to -log(1-r^2), 
    where r is the pairwise partial correlation. (#195)

  * random sampling of R-vine structures. (#197)

  * add `weights` argument to `vine()`. (#188)

  * parallelized fitting of margins in `vine()`. (#198)
  
API BREAK

  * The new `var_types` argument for discrete models has been placed
    early in `bicop()/vinecop()` due to its importance. This might break old
    code calling these functions with unnamed arguments.

BUG FIXES AND OTHER IMPROVEMENTS

  * better support for 0-truncated structures. (#195)

  * ensure consistency of TLL likelihood during and after fit. (#195)

  * fixed order of ranks in `pseudo_obs(.., ties.method = "first")`. (#195)

  * safer computation of multivariate normal cdf. (#195)

  * improved memory efficiency. (#195)


# rvinecopulib 0.3.1.1.0 (July 4, 2019)

Release following the updates of vinecopulib to 0.3.2, see 
https://github.com/vinecopulib/vinecopulib/releases. The most 
relevant changes are summarized below.

NEW FEATURES

   * improved extensibility for packages including on the C++-headers of vinecopulib (#178)

   * new EDA function `pairs_copula_data()` (#181).

BUG FIXES

   * ensure that input and output type of `pseudo_obs()` match (#182).
   
   * fix printing of `"tll"` family in `summary.vinecop() (#183).


# rvinecopulib 0.3.1.1.0 (April 19, 2019)

Release following the updates of vinecopulib to 0.3.1, see 
https://github.com/vinecopulib/vinecopulib/releases. The most 
relevant changes are summarized below.

NEW FEATURES

   * import vinecopulib 0.3.1 (#171)

   * improve documentation (#168, #170)

   * warning message for wrong family in itau method (#169)

   * refactoring for enhanced extensibility of the class `Vinecop` (vinecopulib #407)

   * simplify algorithms by reversing definition of natural order (vinecopulib #387)

   * improve selection of truncation level (vinecopulib #373)

   * add truncate methods for `TriangularArray`, `RVineStructure` and `Vinecop` (vinecopulib #372)

BUG FIXES

   * don't strip debug symbols unconditionally on linux (#174)

# rvinecopulib 0.3.0.1.1 (August 22, 2018)

BUG FIXES

   * fix non-portable use of `log()` in C++ code (#147).
   
   * remove parallelized unit test to avoid segfault on Solaris (#147).

# rvinecopulib 0.3.0.1.0 (August 9, 2018)

Release following the updates of vinecopulib to 0.3.0, see 
https://github.com/vinecopulib/vinecopulib/releases. The most 
relevant changes are summarized below.

NEW FEATURES

   * new generic `truncate_model()` to truncated vine structures and models 
     (#144).

   * new functions `rosenblatt()` and `inverse_rosenblatt()` for computing the 
     (inverse) Rosenblatt transformation (#142).
  
   * faster algorithms for nonparametric copulas based on bilinear interpolation.

   * refactor vine structures and related algorithms with triangular arrays 
     to improve efficiency of truncated models (#136).
     
   * new classes `rvine_structure` and `rvine_matrix` for storing the vine 
     structure including `as_`- and `is.`-generics (#136).

   * allow for generating quasi-random  numbers (#126).

   * improved parallelization: faster of fitting vine copula models and 
     parallelized versions of many algorithms including pdf, cdf and simulation 
     (#339, #363).
     
   * allow weights for observations (#118).

   * faster compilation using only a single wrapper file (#124).

   * improved print and summary generics (#131).

BUG FIXES

   * fix cdf of StudentBicop

   * improved numerical stability.

   * fix gcc-8 warning.
   
   * fix missing variable names for class `vine`.


# rvinecopulib 0.2.8.1.0 (May 8, 2018)

Release following the updates of vinecopulib to 0.2.8, see 
https://github.com/vinecopulib/vinecopulib/releases. The most 
relevant changes are summarized below.

NEW FEATURES
   
   * new `vine_dist` and `vine` classes for data with non-uniform margins (#97).
   
   * new function `truncate_model()` for `vinecop_dist` and `vine_dist` 
   objects (#95, #97).
   
   * new convenience functions `get_pair_copula()`, `get_parameters()`, 
   `get_ktau()`, `get_family()` for `bicop_dist`, `vinecop_dist` and 
   `vine_dist` objects (#95, #107, #109).
   
   * new convenience functions `get_matrix()`, `get_all_pair_copulas()`,
   `get_all_parameters()`, `get_all_ktaus()`, `get_all_families()` for 
   `vinecop_dist` and `vine_dist` objects (#95, #107, #109).
   
   * new (`dim`) and improved (`print`, `summary` and `logLik`) generic methods 
   for `vinecop_dist` and `vine_dist` objects (#104, #109, #110).
   
   * new function `pseudo_obs` to compute pseudo-observations (#108).

   * improved documentation (#98, #100).

   * improved sanity checks and error messages (#99, #102).
   
BUG FIXES
   
   * make mcor correction less aggressive (#103).
   
   * fix truncation of pdf values (#103).
   
   * use increased search interval for parameter estimation when initial fit is 
     unreasonable (#103).
     
   * ensure that boundaries are respected for Joe's `hinv` methods (#103).

   * improve numerical stability by more restrictive parameter bounds for Joe 
     and BB7 copulas (#103). 

# rvinecopulib 0.2.7.1.0 (March 1, 2018)

Release following the updates of vinecopulib to 0.2.7, see 
https://github.com/vinecopulib/vinecopulib/releases. The most 
relevant changes are summarized below.

NEW FEATURES

   * new criterion for tree selection `"mcor"`.

BUG FIXES
   
   * fix bandwidth scaling for family `"tll"`.


# rvinecopulib 0.2.6.1.1 (February 24, 2018)

Patch of rvinecopulib 0.2.6.1.0.

BUG FIXES
   
   * corrected documentation items.


# rvinecopulib 0.2.6.1.0 (February 23, 2018)

Release following the updates of vinecopulib to 0.2.6, see 
https://github.com/vinecopulib/vinecopulib/releases. The most 
relevant changes are summarized below.

NEW FEATURES

   * add checks for data in (0, 1).

   * improved bandwidth selection for method `"tll"` by scaling with 
     maximum-correlation coefficient.
     
   * add mBICV criterion to select the truncation level and threshold along with
     new function `mBICV()`.

BUG FIXES
   
   * fix calculation of Hoeffding's D.


# rvinecopulib 0.2.5.1.0 (January 14, 2017)

Release following the updates of vinecopulib to 0.2.4 and 0.2.5, see 
https://github.com/vinecopulib/vinecopulib/releases. The most 
relevant changes are summarized below.

NEW FEATURES

   * faster simulation and pdf functions for truncated vines.

   * speed up vine copula algorithms by pre-computing information related to 
     the vine structure.
     
   * the selected threshold parameter can be returned from an `vinecop` 
     object.
     
BUG FIXES

   * make bb8 lower bound ensure feasible computations in `par_to_tau()`.

   * default initialize `Rcout` (#277).

   * fix storage order of pair copulas when structure is fixed.
   
   * fixed selection algorithm for threshold and truncation level.


# rvinecopulib 0.2.3.1.0 (November 18, 2017)

Release following the update of vinecopulib to 0.2.3, see 
https://github.com/vinecopulib/vinecopulib/releases/tag/v0.2.3. The most 
relevant changes are summarized below.

NEW FEATURES

   * faster implementation of Archimedean pdfs.

BUG FIXES

   * add safeguards for `bicop()`/`vinecop()` called with 
     insufficient data.

   * fix segfault issue in completing a truncated vine fit.

   * make `par_method = "itau"` respect the parameter bounds.


# rvinecopulib 0.2.2.1.0 (November 9, 2017)

Release following the updates of vinecopulib to 0.2.1 and 0.2.2, see 
https://github.com/vinecopulib/vinecopulib/releases/tag/v0.2.1 and
https://github.com/vinecopulib/vinecopulib/releases/tag/v0.2.2. The most 
relevant changes are summarized below.

NEW FEATURES

   * faster vine copula estimation and selection by parallelizing further 
     sub-routines.

   * enhanced cross-platform compatibility.
        
   * increased precision of maximum-likelihood estimators.
   
   * allow `"loglik"` as selection criterion.
   
BUG FIXES
   
   * fixed `itau` estimation method for Frank copulas (only allowed for positive
     parameters).

  * make interpolation grid symmetric around (0.5, 0.5) again (for `"tll"` 
    estimator).
  

# rvinecopulib 0.2.0.1.0 (October 30, 2017)

Release following the update of vinecopulib to 0.2.0, see 
https://github.com/vinecopulib/vinecopulib/releases/tag/v0.2.0. The most 
relevant changes are summarized below.

PACKAGING/DEPENDECY

   * the C++ core of the library (besides wrappers) is now header only, so 
     other R packages can access its functionality easily via LinkingTo.

   * removed dependency on `NLopt`.

NEW FEATURES
 
   * NA handling.
   
   * parallelized selection/estimation of (pair-) copulas, see the
     `cores` argument in `bicop()` and `vinecop()`.
   
   * efficient storage and fitting of truncated vines.
   
   * Brent line search for (profile-) maximum-likelihood estimation of 
     one-parameter families.
     
   * more restrictive parameter bounds for Archimedean families, ensuring 
     their numerical stability.

BUG FIXES

   * error thrown whenever `vinecop()` or `bicop()` are called with
     data sets containing a single row.
     
   * made order of `rvinecop(..., U)` consistent for d = 2 and d > 2.
   
   * fixed bug in interpolation of kernel estimators near upper right corner.
   
   * interpolation grid is now symmetric around (0.5, 0.5).
   
   * stabilized quadratic tll estimator near zero.
   
   * stabilized Archimedean pdfs.


# rvinecopulib 0.1.0.1.1 (September 1, 2017)

BUG FIXES

   * Improve portability when using mathematical functions

# rvinecopulib 0.1.0.1.0 (August 29, 2017)

Initial release.
