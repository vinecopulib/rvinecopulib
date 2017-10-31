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
