# vinecopulib 0.2.5.1.0 (January 14, 2017)

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


# vinecopulib 0.2.3.1.0 (November 18, 2017)

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


# vinecopulib 0.2.2.1.0 (November 9, 2017)

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
  

# vinecopulib 0.2.0.1.0 (October 30, 2017)

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
