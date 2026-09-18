# Changelog

## rvinecopulib 1.0.0.1.0

The first stable release, based on vinecopulib 1.0.0. Highlights include
a new marginal-modeling interface, conditional simulation and
transforms, analytic derivatives, observation-specific parameters, and
faster evaluation and fitting. See the [vinecopulib 1.0.0
NEWS](https://github.com/vinecopulib/vinecopulib/blob/009a06da1f54dc7690420b5d4c167ba30f32dbca/NEWS.md)
for the backend changes.

#### BREAKING CHANGES

- Require R \>= 4.3.0, C++17, Boost headers from BH \>= 1.75.0-0, and
  wdm \>= 0.3.0. Optional parametric margin fitting requires
  univariateML \>= 1.5.0.

- R-vine structures now follow the backend convention with the
  conditioned variable on the diagonal. Consequently, the matrix, order,
  structure array, and edge orientation representing a model can differ
  from earlier releases; densities and log-likelihoods are unchanged.

- Marginal fitting now uses explicit margin-family and fitted-margin S3
  protocols. Configure KDE options with
  [`kde1d_family()`](https://vinecopulib.github.io/rvinecopulib/reference/kde1d_family.md)
  and variable types with the top-level `var_types` argument. Custom
  [`margin_family()`](https://vinecopulib.github.io/rvinecopulib/reference/margin_family.md)
  fitters now receive `x`, `weights`, and `type` on every call.

#### BEHAVIOR CHANGES

- TLL fits, particularly for discrete data, can change after fixes to
  CDF integration, inversion, and boundary handling.

- Kendall’s tau for the BB6, BB7, BB8, and Tawn families incorporates
  numerical fixes. Maximum-likelihood estimates can also shift slightly
  because of optimizer improvements.

- Compact `d + k` and expanded `2d` layouts for discrete variables are
  handled consistently across evaluation, Rosenblatt transforms, and
  conditional simulation.

#### NEW FEATURES

- Add extensible protocols for fitted margins and margin families,
  including
  [`margin_dist()`](https://vinecopulib.github.io/rvinecopulib/reference/margin_dist.md)
  and
  [`margin_family()`](https://vinecopulib.github.io/rvinecopulib/reference/margin_family.md),
  distribution and quantile generics, and model metadata through
  [`margin_info()`](https://vinecopulib.github.io/rvinecopulib/reference/margin_protocol.md).

- [`dvinecop()`](https://vinecopulib.github.io/rvinecopulib/reference/vinecop_methods.md)
  and
  [`dvine()`](https://vinecopulib.github.io/rvinecopulib/reference/vine_methods.md)
  gain a `log` argument returning the log-density, and
  `dvinecop(keep_all = TRUE)` reports `logpdf` alongside `pdf`. A vine
  density is a product of one factor per edge, so it underflows to `0`
  in high dimensions or under strong dependence while its logarithm is
  still an ordinary double.

- [`dvine()`](https://vinecopulib.github.io/rvinecopulib/reference/vine_methods.md)
  accumulates the margins and the copula density in log space and
  exponentiates once at the end, so it returns a representable density
  where it previously returned `0`: the copula factor used to be
  exponentiated before being multiplied by the marginal densities,
  discarding a joint density that concentrated margins bring back into
  range.

- Add
  [`kde1d_family()`](https://vinecopulib.github.io/rvinecopulib/reference/kde1d_family.md),
  [`univariateML_family()`](https://vinecopulib.github.io/rvinecopulib/reference/univariateML_family.md),
  and
  [`stats_margin()`](https://vinecopulib.github.io/rvinecopulib/reference/stats_margin.md)
  adapters.
  [`vine()`](https://vinecopulib.github.io/rvinecopulib/reference/vine.md)
  can select among these and user-defined candidates through
  `margins_controls$family_set` while reporting and skipping failed
  candidates.

- Add
  [`zero_inflated()`](https://vinecopulib.github.io/rvinecopulib/reference/zero_inflated.md)
  and `var_types` to distinguish continuous, discrete, and zero-inflated
  variables in
  [`vine()`](https://vinecopulib.github.io/rvinecopulib/reference/vine.md).

- Fit margins in parallel on non-Windows systems, controlled separately
  through `margins_controls$cores`.

- Add conditioning-aware structure selection, conditional simulation,
  and Rosenblatt transforms through `conditioning_set`. Conditioning
  values can be common or observation-specific, and fixed, automatically
  selected, and zero truncation levels are supported.

- Add
  [`scores()`](https://vinecopulib.github.io/rvinecopulib/reference/vinecop_methods.md)
  and
  [`hessian()`](https://vinecopulib.github.io/rvinecopulib/reference/vinecop_methods.md)
  for bivariate and vine copulas, and first- and second-order
  derivatives to
  [`dbicop()`](https://vinecopulib.github.io/rvinecopulib/reference/bicop_methods.md)
  and
  [`hbicop()`](https://vinecopulib.github.io/rvinecopulib/reference/bicop_methods.md)
  through `deriv`.

- Add sandwich covariance estimates and Wald confidence intervals for
  fitted bivariate and vine copulas through
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) and
  [`confint()`](https://rdrr.io/r/stats/confint.html).

- Support observation-specific parameters in bivariate copula functions
  and
  [`dvinecop()`](https://vinecopulib.github.io/rvinecopulib/reference/vinecop_methods.md).
  The latter can also return per-edge densities and h-functions through
  `keep_all`.

- Allow custom R functions and symmetrized Chatterjee’s xi (`"cxi"`) as
  tree-selection criteria in
  [`vinecop()`](https://vinecopulib.github.io/rvinecopulib/reference/vinecop.md)
  and
  [`vine()`](https://vinecopulib.github.io/rvinecopulib/reference/vine.md).

- Add
  [`tail_dep()`](https://vinecopulib.github.io/rvinecopulib/reference/bicop_dependence.md)
  and
  [`blomqvist_beta()`](https://vinecopulib.github.io/rvinecopulib/reference/bicop_dependence.md)
  for bivariate copula models and include these dependence summaries in
  printed model output.

#### PERFORMANCE

- Speed up bivariate and vine evaluation, fitting, structure selection,
  pseudo-observations, integration, derivatives, and TLL interpolation.

#### BUG FIXES

- Ensure weighted Wilson random spanning trees terminate when all
  candidate edge strengths are zero.

- `mBICV(object, newdata = )` no longer returns `Inf` when the density
  underflows. It computed the log-likelihood as
  `sum(log(dvinecop(newdata)))`, which collapses as soon as one
  observation underflows to `0`; it now uses the backend’s log-space
  log-likelihood, as the stored `object$loglik` always did.

- Correctly evaluate and refit zero-truncated models and models whose
  omitted pair copulas represent implicit independence.

- Correct the mBICV sparsity prior and require its prior probability
  `psi0` to lie strictly between zero and one.

- Preserve rows with missing unordered-factor values and variable names
  for discrete copula data.

- Do not retain transformed copula data unless `keep_data = TRUE` is
  explicitly requested.

- Return standard `logLik` objects from fitted bivariate copula, vine
  copula, and vine distribution models.

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

- fix handling of discrete variables in
  [`vine()`](https://vinecopulib.github.io/rvinecopulib/reference/vine.md)
  models and related functions.

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

- allow bicop_dist() with tll
  ([\#268](https://github.com/vinecopulib/rvinecopulib/issues/268))

## rvinecopulib 0.6.3.1.1

CRAN release: 2023-02-23

- add `-D_HAS_AUTO_PTR_ETC=0` flag to disable deprecated features used
  in boost.

## rvinecopulib 0.6.3.1.0

CRAN release: 2023-02-20

- fix `NA` handling in `to_pseudo_obs()`
  ([\#260](https://github.com/vinecopulib/rvinecopulib/issues/260))

- add
  [`emp_cdf()`](https://vinecopulib.github.io/rvinecopulib/reference/emp_cdf.md)
  for the tail corrected empirical cdf
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

Maintenance release following the changes to
[`all.equal()`](https://rdrr.io/r/base/all.equal.html) in R 4.1.x.

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

- fix
  [`dvinecop()`](https://vinecopulib.github.io/rvinecopulib/reference/vinecop_methods.md)
  when discrete variables are present
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

- catch na in ktau_to_par
  ([\#214](https://github.com/vinecopulib/rvinecopulib/issues/214))

- make Bicop/Vinecop objects indepent of copied-from-objects

- enforce parameters bounds in tau_to_parameters for Archimedean
  families

## rvinecopulib 0.5.2.1.0 (May 7, 2020)

CRAN release: 2020-05-07

Release following the updates of vinecopulib to 0.5.2, see
<https://github.com/vinecopulib/vinecopulib/releases>. The most relevant
changes are summarized below.

NEW FEATURES

- single-integer constructors for
  [`dvine_structure()`](https://vinecopulib.github.io/rvinecopulib/reference/rvine_structure.md)/[`cvine_structure()`](https://vinecopulib.github.io/rvinecopulib/reference/rvine_structure.md).
  ([\#203](https://github.com/vinecopulib/rvinecopulib/issues/203))

- add `var_names = "hide"` option in
  [`plot.vinecop_dist()`](https://vinecopulib.github.io/rvinecopulib/reference/plot.vinecop_dist.md).
  ([\#203](https://github.com/vinecopulib/rvinecopulib/issues/203))

- add function
  [`plot.rvine_matrix()`](https://vinecopulib.github.io/rvinecopulib/reference/plot.rvine_structure.md).
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

- convenience classes
  [`dvine_structure()`](https://vinecopulib.github.io/rvinecopulib/reference/rvine_structure.md)/[`cvine_structure()`](https://vinecopulib.github.io/rvinecopulib/reference/rvine_structure.md)
  for D- and C-vine structures.
  ([\#195](https://github.com/vinecopulib/rvinecopulib/issues/195))

- new criterion for tree selection: `"joe"` corresponds to -log(1-r^2),
  where r is the pairwise partial correlation.
  ([\#195](https://github.com/vinecopulib/rvinecopulib/issues/195))

- random sampling of R-vine structures.
  ([\#197](https://github.com/vinecopulib/rvinecopulib/issues/197))

- add `weights` argument to
  [`vine()`](https://vinecopulib.github.io/rvinecopulib/reference/vine.md).
  ([\#188](https://github.com/vinecopulib/rvinecopulib/issues/188))

- parallelized fitting of margins in
  [`vine()`](https://vinecopulib.github.io/rvinecopulib/reference/vine.md).
  ([\#198](https://github.com/vinecopulib/rvinecopulib/issues/198))

API BREAK

- The new `var_types` argument for discrete models has been placed early
  in `bicop()/vinecop()` due to its importance. This might break old
  code calling these functions with unnamed arguments.

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

- new EDA function
  [`pairs_copula_data()`](https://vinecopulib.github.io/rvinecopulib/reference/pairs_copula_data.md)
  ([\#181](https://github.com/vinecopulib/rvinecopulib/issues/181)).

BUG FIXES

- ensure that input and output type of
  [`pseudo_obs()`](https://vinecopulib.github.io/rvinecopulib/reference/pseudo_obs.md)
  match
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

- fix non-portable use of [`log()`](https://rdrr.io/r/base/Log.html) in
  C++ code
  ([\#147](https://github.com/vinecopulib/rvinecopulib/issues/147)).

- remove parallelized unit test to avoid segfault on Solaris
  ([\#147](https://github.com/vinecopulib/rvinecopulib/issues/147)).

## rvinecopulib 0.3.0.1.0 (August 9, 2018)

CRAN release: 2018-08-09

Release following the updates of vinecopulib to 0.3.0, see
<https://github.com/vinecopulib/vinecopulib/releases>. The most relevant
changes are summarized below.

NEW FEATURES

- new generic
  [`truncate_model()`](https://vinecopulib.github.io/rvinecopulib/reference/truncate_model.md)
  to truncated vine structures and models
  ([\#144](https://github.com/vinecopulib/rvinecopulib/issues/144)).

- new functions
  [`rosenblatt()`](https://vinecopulib.github.io/rvinecopulib/reference/rosenblatt.md)
  and
  [`inverse_rosenblatt()`](https://vinecopulib.github.io/rvinecopulib/reference/rosenblatt.md)
  for computing the (inverse) Rosenblatt transformation
  ([\#142](https://github.com/vinecopulib/rvinecopulib/issues/142)).

- faster algorithms for nonparametric copulas based on bilinear
  interpolation.

- refactor vine structures and related algorithms with triangular arrays
  to improve efficiency of truncated models
  ([\#136](https://github.com/vinecopulib/rvinecopulib/issues/136)).

- new classes `rvine_structure` and `rvine_matrix` for storing the vine
  structure including `as_`- and `is.`-generics
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

- new function
  [`truncate_model()`](https://vinecopulib.github.io/rvinecopulib/reference/truncate_model.md)
  for `vinecop_dist` and `vine_dist` objects
  ([\#95](https://github.com/vinecopulib/rvinecopulib/issues/95),
  [\#97](https://github.com/vinecopulib/rvinecopulib/issues/97)).

- new convenience functions
  [`get_pair_copula()`](https://vinecopulib.github.io/rvinecopulib/reference/getters.md),
  [`get_parameters()`](https://vinecopulib.github.io/rvinecopulib/reference/getters.md),
  [`get_ktau()`](https://vinecopulib.github.io/rvinecopulib/reference/getters.md),
  [`get_family()`](https://vinecopulib.github.io/rvinecopulib/reference/getters.md)
  for `bicop_dist`, `vinecop_dist` and `vine_dist` objects
  ([\#95](https://github.com/vinecopulib/rvinecopulib/issues/95),
  [\#107](https://github.com/vinecopulib/rvinecopulib/issues/107),
  [\#109](https://github.com/vinecopulib/rvinecopulib/issues/109)).

- new convenience functions
  [`get_matrix()`](https://vinecopulib.github.io/rvinecopulib/reference/getters.md),
  [`get_all_pair_copulas()`](https://vinecopulib.github.io/rvinecopulib/reference/getters.md),
  [`get_all_parameters()`](https://vinecopulib.github.io/rvinecopulib/reference/getters.md),
  [`get_all_ktaus()`](https://vinecopulib.github.io/rvinecopulib/reference/getters.md),
  [`get_all_families()`](https://vinecopulib.github.io/rvinecopulib/reference/getters.md)
  for `vinecop_dist` and `vine_dist` objects
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

- add mBICV criterion to select the truncation level and threshold along
  with new function
  [`mBICV()`](https://vinecopulib.github.io/rvinecopulib/reference/mBICV.md).

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

- add safeguards for
  [`bicop()`](https://vinecopulib.github.io/rvinecopulib/reference/bicop.md)/[`vinecop()`](https://vinecopulib.github.io/rvinecopulib/reference/vinecop.md)
  called with insufficient data.

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

- make interpolation grid symmetric around (0.5, 0.5) again (for `"tll"`
  estimator).

## rvinecopulib 0.2.0.1.0 (October 30, 2017)

CRAN release: 2017-10-30

Release following the update of vinecopulib to 0.2.0, see
<https://github.com/vinecopulib/vinecopulib/releases/tag/v0.2.0>. The
most relevant changes are summarized below.

PACKAGING/DEPENDECY

- the C++ core of the library (besides wrappers) is now header only, so
  other R packages can access its functionality easily via LinkingTo.

- removed dependency on `NLopt`.

NEW FEATURES

- NA handling.

- parallelized selection/estimation of (pair-) copulas, see the `cores`
  argument in
  [`bicop()`](https://vinecopulib.github.io/rvinecopulib/reference/bicop.md)
  and
  [`vinecop()`](https://vinecopulib.github.io/rvinecopulib/reference/vinecop.md).

- efficient storage and fitting of truncated vines.

- Brent line search for (profile-) maximum-likelihood estimation of
  one-parameter families.

- more restrictive parameter bounds for Archimedean families, ensuring
  their numerical stability.

BUG FIXES

- error thrown whenever
  [`vinecop()`](https://vinecopulib.github.io/rvinecopulib/reference/vinecop.md)
  or
  [`bicop()`](https://vinecopulib.github.io/rvinecopulib/reference/bicop.md)
  are called with data sets containing a single row.

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
