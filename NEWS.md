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
   
   * make mcor correction less agressive (#103).
   
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
