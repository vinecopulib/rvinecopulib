# Changelog

## rvinecopulib 1.0.0.1.0

The first stable release. It bundles vinecopulib 1.0.0 and collects a
large backend and frontend update: analytic derivatives, conditional
simulation, conditioning-aware transforms and structure selection,
observation-specific parameters, faster evaluation and fitting, and a
modernized C++17 build. See the [vinecopulib 1.0.0
NEWS](https://github.com/vinecopulib/vinecopulib/blob/009a06da1f54dc7690420b5d4c167ba30f32dbca/NEWS.md)
for the complete backend changes.

#### BREAKING API CHANGES

  - Require R \>= 4.3.0, C++17, and Boost headers from BH \>= 1.75.0-0.

  - R-vine structures now follow the backend convention with the
    conditioned variable on the diagonal. Consequently, the matrix,
    order, structure array, and edge orientation representing a model
    can differ from earlier releases; densities and log-likelihoods are
    unchanged.

  - Marginal fitting now uses explicit margin-family and fitted-margin
    S3 protocols. KDE options (`xmin`, `xmax`, `mult`, `bw`, and `deg`)
    are configured with `kde1d_family()` instead of being entries in
    `margins_controls`; variable types are declared with the top-level
    `var_types` argument. Custom `margin_family()` fitters now receive
    `x`, `weights`, and `type` on every call.

#### BEHAVIOR CHANGES

  - TLL fits can change because the backend no longer clamps
    interpolation grid endpoints at the boundary of the unit square,
    uses recovered latent samples for discrete inputs, and normalizes
    interpolation grids symmetrically.

  - Conditioning-aware selection, Rosenblatt transforms, and conditional
    simulation now support fixed, automatically selected, and zero
    truncation levels. Conditioning-aware selection continues to require
    an MST tree algorithm.

  - Kendall’s tau for the BB6, BB7, BB8, and Tawn families incorporates
    several numerical fixes. Maximum-likelihood estimates can also shift
    in the low digits after the backend optimizer changed from BOBYQA to
    Brent/BFGS.

  - Compact `d + k` and expanded `2d` layouts for discrete variables are
    handled consistently across bivariate and vine-copula evaluation,
    Rosenblatt transforms, and conditional simulation.

#### NEW FEATURES

  - Add fitted-margin generics for distribution evaluation and model
    information, together with the `margin_dist()` constructor. Add the
    complementary `fit_margin()` family protocol. Both fitted margins
    and family specifications expose their metadata through
    `margin_info()`.

  - Add `kde1d_family()`, `univariateML_family()`, and `stats_margin()`
    adapters. Fixed `stats_margin()` objects retain the parameter counts
    of their distributions. Core vine fitting and evaluation use only
    the two protocols and contain no backend-specific dispatch.

  - Add `zero_inflated()` and a top-level `var_types` argument to
    `vine()`; continuous, integer-valued discrete, and zero-inflated
    left-limit CDFs are handled centrally by the margin protocol.

  - Allow `vine()` to select marginal families through
    `margins_controls$family_set`. Parametric families from
    `univariateML` and user-defined `margin_family()` candidates are
    supported alongside the existing `kde1d` margins.

  - Fit margins in parallel with forked processes when `cores > 1` on
    non-Windows systems. Margin fitting remains serial on Windows and
    can be controlled separately with `margins_controls$cores`.
    Stochastic custom fitters can use one margin-fitting core when
    results must be invariant to the core count.

  - Add conditional simulation to `rvinecop()` and `rvine()`.
    Conditioning values can be common or observation-specific, and
    `conditioning_set` accepts indices or names. `vinecop()` and
    `vine()` also accept `conditioning_set` for conditioning-aware
    structure selection. Simulation remains synchronized with R’s
    `set.seed()` for pseudo- and quasi-random generation.

  - Add `conditioning_set` to `rosenblatt()` and `inverse_rosenblatt()`.
    The model is transiently evaluated in a compatible order and is not
    modified.

  - Add `scores()` and `hessian()` for bivariate copulas and extend both
    functions for vine copulas. They support observation-specific
    parameter matrices for continuous parametric models.

  - Add first- and second-order density, log-density, and h-function
    derivatives to `dbicop()` and `hbicop()` through the `deriv`
    argument.

  - Add `keep_all` to `dvinecop()` to return per-edge densities and
    h-functions, and allow observation-specific parameters in
    `dvinecop()`.

  - Support observation-specific parameters passed directly to
    `dbicop()`, `pbicop()`, `hbicop()`, and `rbicop()`. Parameter rows
    are not recycled; `bicop_dist()` objects continue to store one fixed
    parameter set.

  - Allow an R function as `tree_crit` in `vinecop()` and `vine()`. The
    callback is serialized on the calling thread, while pair-copula
    fitting may still use multiple cores.

  - Add `tail_dep()` and `blomqvist_beta()` for bivariate copula models
    and include these dependence summaries in printed model output.

#### PERFORMANCE

  - Incorporate broad backend speedups for bivariate and vine
    evaluation, analytic derivative cascades, TLL fitting/interpolation,
    structure selection, pseudo-observations, integration, and shared
    Eigen/thread primitives.

#### BUG FIXES

  - Validate and report failed marginal candidates consistently,
    preserve their warnings, reject fitted supports that exclude
    observations, and detect failed forked workers and malformed margin
    core counts early.

  - Make inverse Rosenblatt transforms thread-safe and custom tree
    criteria safe under multithreaded fitting.

  - Do not retain transformed copula data when partial `copula_controls`
    omit `keep_data`; retaining data now requires an explicit opt-in.

  - Preserve rows with missing unordered-factor values when expanding
    factors in `vine()`.

  - Fix TLL CDF integration and inversion, use the recovered latent
    sample when fitting discrete TLL models, and evaluate narrow
    discrete intervals at their midpoint.

  - Ensure weighted Wilson random spanning trees terminate when all
    candidate edge strengths are zero.

  - Correctly evaluate and refit zero-truncated models and models whose
    omitted pair copulas represent implicit independence.

  - Fix starting parameters for discrete models and edge-case per-row
    parameter shapes.

  - Preserve and correctly trim variable names for discrete copula data.

  - Return standard `logLik` objects from fitted bivariate copula, vine
    copula, and vine distribution models.

#### BUILD SYSTEM AND DEPENDENCIES

  - Vendor the complete vinecopulib 1.0.0 header tree, while keeping the
    package wrapper header as the R-specific integration point so
    downstream packages can include all public backend headers.

  - Compile the package as C++17 and require BH \>= 1.75.0-0. The
    backend update also reduces its Boost surface to Graph, Math, and
    Random.

  - Require univariateML \>= 1.5.0 for optional parametric margin
    support.

  - Modernize `src/update_vinecopulib.sh`: it accepts a branch, tag, or
    exact commit, imports through a temporary checkout, preserves
    package-owned wrappers, copies all public headers, and records the
    imported commit in the vendored tree.

## rvinecopulib 0.7.3.1.0

CRAN release: 2025-06-13

#### NEW FEATURES

  - Allow for random spanning trees as alternatives to the MST-based
    structure selection using `tree_algorithm` in `vine` and `vinecop`.
    Options are `"mst_prim"`, `"mst_kruskal"`, `"random_weighted"` or
    `"random_unweighted"`
    ([\#307](https://github.com/vinecopulib/rvinecopulib/issues/307)).

#### BUG FIXES

  - Decouple edge insertion from criterion computation fix randomness
    issues in structure selection when using multiple threads
    ([\#640](https://github.com/vinecopulib/vinecopulib/pull/640))

## rvinecopulib 0.7.2.1.0

CRAN release: 2025-03-24

BUG FIX

  - fix TLL speed issues related to FFT
    ([\#305](https://github.com/vinecopulib/rvinecopulib/issues/305)).

## rvinecopulib 0.7.1.1.2

CRAN release: 2025-03-03

BUG FIX

  - Fixes “deprecated-literal-operator” warning on clang20.

## rvinecopulib 0.7.1.1.1

CRAN release: 2025-02-10

BUG FIX

  - fix handling of discrete variables in `vine()` models and related
    functions.

## rvinecopulib 0.7.1.1.0

CRAN release: 2025-01-21

Update following a upgrade of the C++ backend vinecopulib to 0.7.1, see
<https://github.com/vinecopulib/vinecopulib/blob/main/NEWS.md>.

The main changes on the R end are:

  - improved documentation,

  - support for zero-inflated variables,

  - added new Tawn copula family,

  - new argument `allow_rotations` to disable rotations of copula
    families,

  - added variable names to vinecop summary
    ([\#276](https://github.com/vinecopulib/rvinecopulib/issues/276))

  - fixed handling of logistic distribution
    ([\#275](https://github.com/vinecopulib/rvinecopulib/issues/275))

  - fix NA handling in vine() control checks
    ([\#266](https://github.com/vinecopulib/rvinecopulib/issues/266))

  - allow bicop\_dist() with tll
    ([\#268](https://github.com/vinecopulib/rvinecopulib/issues/268))

## rvinecopulib 0.6.3.1.1

CRAN release: 2023-02-23

  - add `-D_HAS_AUTO_PTR_ETC=0` flag to disable deprecated features used
    in boost.

## rvinecopulib 0.6.3.1.0

CRAN release: 2023-02-20

  - fix `NA` handling in `to_pseudo_obs()`
    ([\#260](https://github.com/vinecopulib/rvinecopulib/issues/260))

  - add `emp_cdf()` for the tail corrected empirical cdf
    ([\#261](https://github.com/vinecopulib/rvinecopulib/issues/261))

## rvinecopulib 0.6.2.1.3 (December 3, 2022)

CRAN release: 2022-12-06

  - fix marginal PIT for discrete variables (see issue
    [\#257](https://github.com/vinecopulib/rvinecopulib/issues/257),
    thanks [@rplzzz](https://github.com/rplzzz))

## rvinecopulib 0.6.2.1.2 (October 16, 2022)

CRAN release: 2022-10-17

  - fix warning about C++17 attribute extension ‘nodiscard’
    ([\#255](https://github.com/vinecopulib/rvinecopulib/issues/255))

## rvinecopulib 0.6.2.1.1 (August 30, 2022)

CRAN release: 2022-09-04

  - replace bitwise operations on Boolean variables.

## rvinecopulib 0.6.2.1.0 (August 26, 2022)

CRAN release: 2022-08-26

Release following the updates of vinecopulib to 0.6.2, see
<https://github.com/vinecopulib/vinecopulib/releases>. The most relevant
changes are summarized below.

  - improved documentation (discrete and missing dat, Rosenblatt
    transforms)

  - better parallelization when there is a small number of edges
    ([\#555](https://github.com/vinecopulib/rvinecopulib/issues/555))

## rvinecopulib 0.6.1.1.3 (March 17, 2022)

CRAN release: 2022-03-18

  - avoid bitwise operations on Boolean variables.

## rvinecopulib 0.6.1.1.2 (March 14, 2022)

CRAN release: 2022-03-14

  - enforce use on recent version of RcppThread for proper handling of
    linker flags

## rvinecopulib 0.6.1.1.1 (October 6, 2021)

CRAN release: 2021-10-07

#### BUG FIXES

  - remove illegal pragmas from json header
    ([\#245](https://github.com/vinecopulib/rvinecopulib/issues/245))

  - allow tree restriction in summary functions
    ([\#244](https://github.com/vinecopulib/rvinecopulib/issues/244))

## rvinecopulib 0.6.1.1.0 (July 13, 2021)

Release following the updates of vinecopulib to 0.6.1, see
<https://github.com/vinecopulib/vinecopulib/releases>. The most relevant
changes are summarized below.

#### BUG FIXES

  - improved documentation
    ([\#241](https://github.com/vinecopulib/rvinecopulib/issues/241),
    [\#239](https://github.com/vinecopulib/rvinecopulib/issues/239))

  - use `num_threads` in recursive calls to the inverse Rosenblatt

  - force TLL to be nonnegative
    ([\#238](https://github.com/vinecopulib/rvinecopulib/issues/238))

  - fix number of parameters for TLL

## rvinecopulib 0.5.5.1.1 (December 15, 2020)

CRAN release: 2021-01-06

Maintenance release following the changes to `all.equal()` in R 4.1.x.

## rvinecopulib 0.5.5.1.0 (November 24, 2020)

CRAN release: 2020-11-23

Release following the updates of vinecopulib to 0.5.5, see
<https://github.com/vinecopulib/vinecopulib/releases>. The most relevant
changes are summarized below.

#### BUG FIXES

  - fix little bug in copula selection based on mBIC

  - stabilize BB7 copula pdf

  - fix threshold selection for (near-)independent data

  - fix vine copula selection for 1-dimensional models with discrete
    variables

  - fix user-visible variable types

## rvinecopulib 0.5.4.1.0 (September 30, 2020)

CRAN release: 2020-10-03

Release following the updates of vinecopulib to 0.5.4, see
<https://github.com/vinecopulib/vinecopulib/releases>. The most relevant
changes are summarized below.

#### BUG FIXES

  - fix uninitialized number of parameters for TLL family

  - fix Kendall’s tau of Frank copula for par \<= 3

  - fix `dvinecop()` when discrete variables are present
    ([\#222](https://github.com/vinecopulib/rvinecopulib/issues/222))

## rvinecopulib 0.5.3.1.0 (August 11, 2020)

CRAN release: 2020-08-12

Release following the updates of vinecopulib to 0.5.3, see
<https://github.com/vinecopulib/vinecopulib/releases>. The most relevant
changes are summarized below.

#### NEW FEATURES

  - allow 1-dimensional models
    ([\#216](https://github.com/vinecopulib/rvinecopulib/issues/216))

  - make AIC default selection criterion
    ([\#213](https://github.com/vinecopulib/rvinecopulib/issues/213))

#### BUG FIXES

  - catch na in ktau\_to\_par
    ([\#214](https://github.com/vinecopulib/rvinecopulib/issues/214))

  - make Bicop/Vinecop objects indepent of copied-from-objects

  - enforce parameters bounds in tau\_to\_parameters for Archimedean
    families

## rvinecopulib 0.5.2.1.0 (May 7, 2020)

CRAN release: 2020-05-07

Release following the updates of vinecopulib to 0.5.2, see
<https://github.com/vinecopulib/vinecopulib/releases>. The most relevant
changes are summarized below.

NEW FEATURES

  - single-integer constructors for
    `dvine_structure()`/`cvine_structure()`.
    ([\#203](https://github.com/vinecopulib/rvinecopulib/issues/203))

  - add `var_names = "hide"` option in `plot.vinecop_dist()`.
    ([\#203](https://github.com/vinecopulib/rvinecopulib/issues/203))

  - add function `plot.rvine_matrix()`.
    ([\#203](https://github.com/vinecopulib/rvinecopulib/issues/203))

BUG FIXES

  - fix bug for (negative) tau to parameter conversion for Frank family.
    ([\#207](https://github.com/vinecopulib/rvinecopulib/issues/207))

  - fix rare error `rvine_structure_sim()/rvine_matrix_sim()`

  - safeguard `"tll"` family against comonotonic data.

  - stabilize archimedean h-functions near independence.

## rvinecopulib 0.5.1.1.0 (November 25, 2019)

CRAN release: 2019-11-26

Release following the updates of vinecopulib to 0.5.1, see
<https://github.com/vinecopulib/vinecopulib/releases>.

BUG FIX

  - fix out of range bug for weighted TLL influence when sample size is
    small.

## rvinecopulib 0.5.0.1.0 (November 25, 2019)

CRAN release: 2019-11-25

Release following the updates of vinecopulib to 0.5.0, see
<https://github.com/vinecopulib/vinecopulib/releases>. The most relevant
changes are summarized below.

NEW FEATURES

  - modelling discrete variables with bivariate or vine copulas.
    ([\#195](https://github.com/vinecopulib/rvinecopulib/issues/195))

  - selection of partially specified R-vine structures.
    ([\#195](https://github.com/vinecopulib/rvinecopulib/issues/195))

  - convenience classes `dvine_structure()`/`cvine_structure()` for D-
    and C-vine structures.
    ([\#195](https://github.com/vinecopulib/rvinecopulib/issues/195))

  - new criterion for tree selection: `"joe"` corresponds to
    -log(1-r^2), where r is the pairwise partial correlation.
    ([\#195](https://github.com/vinecopulib/rvinecopulib/issues/195))

  - random sampling of R-vine structures.
    ([\#197](https://github.com/vinecopulib/rvinecopulib/issues/197))

  - add `weights` argument to `vine()`.
    ([\#188](https://github.com/vinecopulib/rvinecopulib/issues/188))

  - parallelized fitting of margins in `vine()`.
    ([\#198](https://github.com/vinecopulib/rvinecopulib/issues/198))

API BREAK

  - The new `var_types` argument for discrete models has been placed
    early in `bicop()/vinecop()` due to its importance. This might break
    old code calling these functions with unnamed arguments.

BUG FIXES AND OTHER IMPROVEMENTS

  - better support for 0-truncated structures.
    ([\#195](https://github.com/vinecopulib/rvinecopulib/issues/195))

  - ensure consistency of TLL likelihood during and after fit.
    ([\#195](https://github.com/vinecopulib/rvinecopulib/issues/195))

  - fixed order of ranks in `pseudo_obs(.., ties.method = "first")`.
    ([\#195](https://github.com/vinecopulib/rvinecopulib/issues/195))

  - safer computation of multivariate normal cdf.
    ([\#195](https://github.com/vinecopulib/rvinecopulib/issues/195))

  - improved memory efficiency.
    ([\#195](https://github.com/vinecopulib/rvinecopulib/issues/195))

## rvinecopulib 0.3.1.1.0 (July 4, 2019)

CRAN release: 2019-04-19

Release following the updates of vinecopulib to 0.3.2, see
<https://github.com/vinecopulib/vinecopulib/releases>. The most relevant
changes are summarized below.

NEW FEATURES

  - improved extensibility for packages including on the C++-headers of
    vinecopulib
    ([\#178](https://github.com/vinecopulib/rvinecopulib/issues/178))

  - new EDA function `pairs_copula_data()`
    ([\#181](https://github.com/vinecopulib/rvinecopulib/issues/181)).

BUG FIXES

  - ensure that input and output type of `pseudo_obs()` match
    ([\#182](https://github.com/vinecopulib/rvinecopulib/issues/182)).

  - fix printing of `"tll"` family in \`summary.vinecop()
    ([\#183](https://github.com/vinecopulib/rvinecopulib/issues/183)).

## rvinecopulib 0.3.1.1.0 (April 19, 2019)

CRAN release: 2019-04-19

Release following the updates of vinecopulib to 0.3.1, see
<https://github.com/vinecopulib/vinecopulib/releases>. The most relevant
changes are summarized below.

NEW FEATURES

  - import vinecopulib 0.3.1
    ([\#171](https://github.com/vinecopulib/rvinecopulib/issues/171))

  - improve documentation
    ([\#168](https://github.com/vinecopulib/rvinecopulib/issues/168),
    [\#170](https://github.com/vinecopulib/rvinecopulib/issues/170))

  - warning message for wrong family in itau method
    ([\#169](https://github.com/vinecopulib/rvinecopulib/issues/169))

  - refactoring for enhanced extensibility of the class `Vinecop`
    (vinecopulib
    [\#407](https://github.com/vinecopulib/rvinecopulib/issues/407))

  - simplify algorithms by reversing definition of natural order
    (vinecopulib
    [\#387](https://github.com/vinecopulib/rvinecopulib/issues/387))

  - improve selection of truncation level (vinecopulib
    [\#373](https://github.com/vinecopulib/rvinecopulib/issues/373))

  - add truncate methods for `TriangularArray`, `RVineStructure` and
    `Vinecop` (vinecopulib
    [\#372](https://github.com/vinecopulib/rvinecopulib/issues/372))

BUG FIXES

  - don’t strip debug symbols unconditionally on linux
    ([\#174](https://github.com/vinecopulib/rvinecopulib/issues/174))

## rvinecopulib 0.3.0.1.1 (August 22, 2018)

CRAN release: 2018-08-22

BUG FIXES

  - fix non-portable use of `log()` in C++ code
    ([\#147](https://github.com/vinecopulib/rvinecopulib/issues/147)).

  - remove parallelized unit test to avoid segfault on Solaris
    ([\#147](https://github.com/vinecopulib/rvinecopulib/issues/147)).

## rvinecopulib 0.3.0.1.0 (August 9, 2018)

CRAN release: 2018-08-09

Release following the updates of vinecopulib to 0.3.0, see
<https://github.com/vinecopulib/vinecopulib/releases>. The most relevant
changes are summarized below.

NEW FEATURES

  - new generic `truncate_model()` to truncated vine structures and
    models
    ([\#144](https://github.com/vinecopulib/rvinecopulib/issues/144)).

  - new functions `rosenblatt()` and `inverse_rosenblatt()` for
    computing the (inverse) Rosenblatt transformation
    ([\#142](https://github.com/vinecopulib/rvinecopulib/issues/142)).

  - faster algorithms for nonparametric copulas based on bilinear
    interpolation.

  - refactor vine structures and related algorithms with triangular
    arrays to improve efficiency of truncated models
    ([\#136](https://github.com/vinecopulib/rvinecopulib/issues/136)).

  - new classes `rvine_structure` and `rvine_matrix` for storing the
    vine structure including `as_`- and `is.`-generics
    ([\#136](https://github.com/vinecopulib/rvinecopulib/issues/136)).

  - allow for generating quasi-random numbers
    ([\#126](https://github.com/vinecopulib/rvinecopulib/issues/126)).

  - improved parallelization: faster of fitting vine copula models and
    parallelized versions of many algorithms including pdf, cdf and
    simulation
    ([\#339](https://github.com/vinecopulib/rvinecopulib/issues/339),
    [\#363](https://github.com/vinecopulib/rvinecopulib/issues/363)).

  - allow weights for observations
    ([\#118](https://github.com/vinecopulib/rvinecopulib/issues/118)).

  - faster compilation using only a single wrapper file
    ([\#124](https://github.com/vinecopulib/rvinecopulib/issues/124)).

  - improved print and summary generics
    ([\#131](https://github.com/vinecopulib/rvinecopulib/issues/131)).

BUG FIXES

  - fix cdf of StudentBicop

  - improved numerical stability.

  - fix gcc-8 warning.

  - fix missing variable names for class `vine`.

## rvinecopulib 0.2.8.1.0 (May 8, 2018)

CRAN release: 2018-05-09

Release following the updates of vinecopulib to 0.2.8, see
<https://github.com/vinecopulib/vinecopulib/releases>. The most relevant
changes are summarized below.

NEW FEATURES

  - new `vine_dist` and `vine` classes for data with non-uniform margins
    ([\#97](https://github.com/vinecopulib/rvinecopulib/issues/97)).

  - new function `truncate_model()` for `vinecop_dist` and `vine_dist`
    objects
    ([\#95](https://github.com/vinecopulib/rvinecopulib/issues/95),
    [\#97](https://github.com/vinecopulib/rvinecopulib/issues/97)).

  - new convenience functions `get_pair_copula()`, `get_parameters()`,
    `get_ktau()`, `get_family()` for `bicop_dist`, `vinecop_dist` and
    `vine_dist` objects
    ([\#95](https://github.com/vinecopulib/rvinecopulib/issues/95),
    [\#107](https://github.com/vinecopulib/rvinecopulib/issues/107),
    [\#109](https://github.com/vinecopulib/rvinecopulib/issues/109)).

  - new convenience functions `get_matrix()`, `get_all_pair_copulas()`,
    `get_all_parameters()`, `get_all_ktaus()`, `get_all_families()` for
    `vinecop_dist` and `vine_dist` objects
    ([\#95](https://github.com/vinecopulib/rvinecopulib/issues/95),
    [\#107](https://github.com/vinecopulib/rvinecopulib/issues/107),
    [\#109](https://github.com/vinecopulib/rvinecopulib/issues/109)).

  - new (`dim`) and improved (`print`, `summary` and `logLik`) generic
    methods for `vinecop_dist` and `vine_dist` objects
    ([\#104](https://github.com/vinecopulib/rvinecopulib/issues/104),
    [\#109](https://github.com/vinecopulib/rvinecopulib/issues/109),
    [\#110](https://github.com/vinecopulib/rvinecopulib/issues/110)).

  - new function `pseudo_obs` to compute pseudo-observations
    ([\#108](https://github.com/vinecopulib/rvinecopulib/issues/108)).

  - improved documentation
    ([\#98](https://github.com/vinecopulib/rvinecopulib/issues/98),
    [\#100](https://github.com/vinecopulib/rvinecopulib/issues/100)).

  - improved sanity checks and error messages
    ([\#99](https://github.com/vinecopulib/rvinecopulib/issues/99),
    [\#102](https://github.com/vinecopulib/rvinecopulib/issues/102)).

BUG FIXES

  - make mcor correction less aggressive
    ([\#103](https://github.com/vinecopulib/rvinecopulib/issues/103)).

  - fix truncation of pdf values
    ([\#103](https://github.com/vinecopulib/rvinecopulib/issues/103)).

  - use increased search interval for parameter estimation when initial
    fit is unreasonable
    ([\#103](https://github.com/vinecopulib/rvinecopulib/issues/103)).

  - ensure that boundaries are respected for Joe’s `hinv` methods
    ([\#103](https://github.com/vinecopulib/rvinecopulib/issues/103)).

  - improve numerical stability by more restrictive parameter bounds for
    Joe and BB7 copulas
    ([\#103](https://github.com/vinecopulib/rvinecopulib/issues/103)).

## rvinecopulib 0.2.7.1.0 (March 1, 2018)

CRAN release: 2018-03-02

Release following the updates of vinecopulib to 0.2.7, see
<https://github.com/vinecopulib/vinecopulib/releases>. The most relevant
changes are summarized below.

NEW FEATURES

  - new criterion for tree selection `"mcor"`.

BUG FIXES

  - fix bandwidth scaling for family `"tll"`.

## rvinecopulib 0.2.6.1.1 (February 24, 2018)

CRAN release: 2018-02-25

Patch of rvinecopulib 0.2.6.1.0.

BUG FIXES

  - corrected documentation items.

## rvinecopulib 0.2.6.1.0 (February 23, 2018)

Release following the updates of vinecopulib to 0.2.6, see
<https://github.com/vinecopulib/vinecopulib/releases>. The most relevant
changes are summarized below.

NEW FEATURES

  - add checks for data in (0, 1).

  - improved bandwidth selection for method `"tll"` by scaling with
    maximum-correlation coefficient.

  - add mBICV criterion to select the truncation level and threshold
    along with new function `mBICV()`.

BUG FIXES

  - fix calculation of Hoeffding’s D.

## rvinecopulib 0.2.5.1.0 (January 14, 2017)

CRAN release: 2018-01-14

Release following the updates of vinecopulib to 0.2.4 and 0.2.5, see
<https://github.com/vinecopulib/vinecopulib/releases>. The most relevant
changes are summarized below.

NEW FEATURES

  - faster simulation and pdf functions for truncated vines.

  - speed up vine copula algorithms by pre-computing information related
    to the vine structure.

  - the selected threshold parameter can be returned from an `vinecop`
    object.

BUG FIXES

  - make bb8 lower bound ensure feasible computations in `par_to_tau()`.

  - default initialize `Rcout`
    ([\#277](https://github.com/vinecopulib/rvinecopulib/issues/277)).

  - fix storage order of pair copulas when structure is fixed.

  - fixed selection algorithm for threshold and truncation level.

## rvinecopulib 0.2.3.1.0 (November 18, 2017)

CRAN release: 2017-11-18

Release following the update of vinecopulib to 0.2.3, see
<https://github.com/vinecopulib/vinecopulib/releases/tag/v0.2.3>. The
most relevant changes are summarized below.

NEW FEATURES

  - faster implementation of Archimedean pdfs.

BUG FIXES

  - add safeguards for `bicop()`/`vinecop()` called with insufficient
    data.

  - fix segfault issue in completing a truncated vine fit.

  - make `par_method = "itau"` respect the parameter bounds.

## rvinecopulib 0.2.2.1.0 (November 9, 2017)

Release following the updates of vinecopulib to 0.2.1 and 0.2.2, see
<https://github.com/vinecopulib/vinecopulib/releases/tag/v0.2.1> and
<https://github.com/vinecopulib/vinecopulib/releases/tag/v0.2.2>. The
most relevant changes are summarized below.

NEW FEATURES

  - faster vine copula estimation and selection by parallelizing further
    sub-routines.

  - enhanced cross-platform compatibility.

  - increased precision of maximum-likelihood estimators.

  - allow `"loglik"` as selection criterion.

BUG FIXES

  - fixed `itau` estimation method for Frank copulas (only allowed for
    positive parameters).

  - make interpolation grid symmetric around (0.5, 0.5) again (for
    `"tll"` estimator).

## rvinecopulib 0.2.0.1.0 (October 30, 2017)

CRAN release: 2017-10-30

Release following the update of vinecopulib to 0.2.0, see
<https://github.com/vinecopulib/vinecopulib/releases/tag/v0.2.0>. The
most relevant changes are summarized below.

PACKAGING/DEPENDECY

  - the C++ core of the library (besides wrappers) is now header only,
    so other R packages can access its functionality easily via
    LinkingTo.

  - removed dependency on `NLopt`.

NEW FEATURES

  - NA handling.

  - parallelized selection/estimation of (pair-) copulas, see the
    `cores` argument in `bicop()` and `vinecop()`.

  - efficient storage and fitting of truncated vines.

  - Brent line search for (profile-) maximum-likelihood estimation of
    one-parameter families.

  - more restrictive parameter bounds for Archimedean families, ensuring
    their numerical stability.

BUG FIXES

  - error thrown whenever `vinecop()` or `bicop()` are called with data
    sets containing a single row.

  - made order of `rvinecop(..., U)` consistent for d = 2 and d \> 2.

  - fixed bug in interpolation of kernel estimators near upper right
    corner.

  - interpolation grid is now symmetric around (0.5, 0.5).

  - stabilized quadratic tll estimator near zero.

  - stabilized Archimedean pdfs.

## rvinecopulib 0.1.0.1.1 (September 1, 2017)

CRAN release: 2017-09-01

BUG FIXES

  - Improve portability when using mathematical functions

## rvinecopulib 0.1.0.1.0 (August 29, 2017)

CRAN release: 2017-08-30

Initial release.
