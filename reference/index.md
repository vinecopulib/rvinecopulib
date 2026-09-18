# Package index

## Full distributions and margins

- [`vine()`](https://vinecopulib.github.io/rvinecopulib/reference/vine.md)
  [`vine_dist()`](https://vinecopulib.github.io/rvinecopulib/reference/vine.md)
  : Vine copula models
- [`dvine()`](https://vinecopulib.github.io/rvinecopulib/reference/vine_methods.md)
  [`pvine()`](https://vinecopulib.github.io/rvinecopulib/reference/vine_methods.md)
  [`rvine()`](https://vinecopulib.github.io/rvinecopulib/reference/vine_methods.md)
  : Vine based distributions
- [`predict(`*`<vine>`*`)`](https://vinecopulib.github.io/rvinecopulib/reference/predict_vine.md)
  [`fitted(`*`<vine>`*`)`](https://vinecopulib.github.io/rvinecopulib/reference/predict_vine.md)
  : Predictions and fitted values for a vine copula model
- [`as_margin()`](https://vinecopulib.github.io/rvinecopulib/reference/as_margin.md)
  : Normalize an object to the fitted-margin protocol
- [`margin_dist()`](https://vinecopulib.github.io/rvinecopulib/reference/margin_dist.md)
  : Create a fitted or fixed custom margin
- [`stats_margin()`](https://vinecopulib.github.io/rvinecopulib/reference/stats_margin.md)
  : Create a fixed margin from a stats distribution
- [`margin_family()`](https://vinecopulib.github.io/rvinecopulib/reference/margin_family.md)
  : Define a custom margin family
- [`kde1d_family()`](https://vinecopulib.github.io/rvinecopulib/reference/kde1d_family.md)
  : Define a kde1d margin family
- [`univariateML_family()`](https://vinecopulib.github.io/rvinecopulib/reference/univariateML_family.md)
  : Define a univariateML margin family
- [`dmargin()`](https://vinecopulib.github.io/rvinecopulib/reference/margin_protocol.md)
  [`pmargin()`](https://vinecopulib.github.io/rvinecopulib/reference/margin_protocol.md)
  [`qmargin()`](https://vinecopulib.github.io/rvinecopulib/reference/margin_protocol.md)
  [`margin_info()`](https://vinecopulib.github.io/rvinecopulib/reference/margin_protocol.md)
  : Fitted marginal distribution protocol
- [`fit_margin()`](https://vinecopulib.github.io/rvinecopulib/reference/margin_family_protocol.md)
  : Margin-family fitting protocol
- [`zero_inflated()`](https://vinecopulib.github.io/rvinecopulib/reference/zero_inflated.md)
  : Declare zero-inflated data

## Vine copula models

- [`vinecop()`](https://vinecopulib.github.io/rvinecopulib/reference/vinecop.md)
  : Fitting vine copula models
- [`vinecop_dist()`](https://vinecopulib.github.io/rvinecopulib/reference/vinecop_dist.md)
  : Vine copula models
- [`dvinecop()`](https://vinecopulib.github.io/rvinecopulib/reference/vinecop_methods.md)
  [`scores()`](https://vinecopulib.github.io/rvinecopulib/reference/vinecop_methods.md)
  [`hessian()`](https://vinecopulib.github.io/rvinecopulib/reference/vinecop_methods.md)
  [`pvinecop()`](https://vinecopulib.github.io/rvinecopulib/reference/vinecop_methods.md)
  [`rvinecop()`](https://vinecopulib.github.io/rvinecopulib/reference/vinecop_methods.md)
  : Vine copula distributions
- [`predict(`*`<vinecop>`*`)`](https://vinecopulib.github.io/rvinecopulib/reference/predict_vinecop.md)
  [`fitted(`*`<vinecop>`*`)`](https://vinecopulib.github.io/rvinecopulib/reference/predict_vinecop.md)
  : Predictions and fitted values for a vine copula model
- [`rosenblatt()`](https://vinecopulib.github.io/rvinecopulib/reference/rosenblatt.md)
  [`inverse_rosenblatt()`](https://vinecopulib.github.io/rvinecopulib/reference/rosenblatt.md)
  : Rosenblatt and inverse Rosenblatt transforms
- [`truncate_model()`](https://vinecopulib.github.io/rvinecopulib/reference/truncate_model.md)
  : Truncate a vine copula model
- [`mBICV()`](https://vinecopulib.github.io/rvinecopulib/reference/mBICV.md)
  : Modified vine copula Bayesian information criterion (mBICv)
- [`vcov(`*`<vinecop_dist>`*`)`](https://vinecopulib.github.io/rvinecopulib/reference/parameter_uncertainty.md)
  [`vcov(`*`<bicop_dist>`*`)`](https://vinecopulib.github.io/rvinecopulib/reference/parameter_uncertainty.md)
  [`confint(`*`<vinecop_dist>`*`)`](https://vinecopulib.github.io/rvinecopulib/reference/parameter_uncertainty.md)
  [`confint(`*`<bicop_dist>`*`)`](https://vinecopulib.github.io/rvinecopulib/reference/parameter_uncertainty.md)
  : Parameter uncertainty

## Bivariate copula models

- [`bicop()`](https://vinecopulib.github.io/rvinecopulib/reference/bicop.md)
  : Fit and select bivariate copula models

- [`bicop_dist()`](https://vinecopulib.github.io/rvinecopulib/reference/bicop_dist.md)
  : Bivariate copula models

- [`dbicop()`](https://vinecopulib.github.io/rvinecopulib/reference/bicop_methods.md)
  [`pbicop()`](https://vinecopulib.github.io/rvinecopulib/reference/bicop_methods.md)
  [`rbicop()`](https://vinecopulib.github.io/rvinecopulib/reference/bicop_methods.md)
  [`scores(`*`<bicop_dist>`*`)`](https://vinecopulib.github.io/rvinecopulib/reference/bicop_methods.md)
  [`hessian(`*`<bicop_dist>`*`)`](https://vinecopulib.github.io/rvinecopulib/reference/bicop_methods.md)
  [`hbicop()`](https://vinecopulib.github.io/rvinecopulib/reference/bicop_methods.md)
  : Bivariate copula distributions

- [`tail_dep()`](https://vinecopulib.github.io/rvinecopulib/reference/bicop_dependence.md)
  [`blomqvist_beta()`](https://vinecopulib.github.io/rvinecopulib/reference/bicop_dependence.md)
  : Dependence measures of a bivariate copula

- [`predict(`*`<bicop_dist>`*`)`](https://vinecopulib.github.io/rvinecopulib/reference/predict_bicop.md)
  [`fitted(`*`<bicop>`*`)`](https://vinecopulib.github.io/rvinecopulib/reference/predict_bicop.md)
  : Predictions and fitted values for a bivariate copula model

- [`plot(`*`<bicop_dist>`*`)`](https://vinecopulib.github.io/rvinecopulib/reference/plot.bicop_dist.md)
  [`plot(`*`<bicop>`*`)`](https://vinecopulib.github.io/rvinecopulib/reference/plot.bicop_dist.md)
  [`contour(`*`<bicop_dist>`*`)`](https://vinecopulib.github.io/rvinecopulib/reference/plot.bicop_dist.md)
  [`contour(`*`<bicop>`*`)`](https://vinecopulib.github.io/rvinecopulib/reference/plot.bicop_dist.md)
  :

  Plotting tools for `bicop_dist` and `bicop` objects

- [`par_to_ktau()`](https://vinecopulib.github.io/rvinecopulib/reference/par_to_ktau.md)
  [`ktau_to_par()`](https://vinecopulib.github.io/rvinecopulib/reference/par_to_ktau.md)
  : Conversion between Kendall's tau and parameters

- [`as.bicop()`](https://vinecopulib.github.io/rvinecopulib/reference/as.bicop.md)
  : Convert list to bicop object

## R-vine structures

- [`rvine_structure()`](https://vinecopulib.github.io/rvinecopulib/reference/rvine_structure.md)
  [`cvine_structure()`](https://vinecopulib.github.io/rvinecopulib/reference/rvine_structure.md)
  [`dvine_structure()`](https://vinecopulib.github.io/rvinecopulib/reference/rvine_structure.md)
  [`rvine_matrix()`](https://vinecopulib.github.io/rvinecopulib/reference/rvine_structure.md)
  : R-vine structure

- [`as_rvine_structure()`](https://vinecopulib.github.io/rvinecopulib/reference/as_rvine_structure.md)
  [`as_rvine_matrix()`](https://vinecopulib.github.io/rvinecopulib/reference/as_rvine_structure.md)
  : Coerce various kind of objects to R-vine structures and matrices

- [`rvine_structure_sim()`](https://vinecopulib.github.io/rvinecopulib/reference/rvine_structure_sim.md)
  [`rvine_matrix_sim()`](https://vinecopulib.github.io/rvinecopulib/reference/rvine_structure_sim.md)
  : Simulate R-vine structures

- [`plot(`*`<rvine_structure>`*`)`](https://vinecopulib.github.io/rvinecopulib/reference/plot.rvine_structure.md)
  [`plot(`*`<rvine_matrix>`*`)`](https://vinecopulib.github.io/rvinecopulib/reference/plot.rvine_structure.md)
  : Plotting R-vine structures

- [`get_structure()`](https://vinecopulib.github.io/rvinecopulib/reference/getters.md)
  [`get_pair_copula()`](https://vinecopulib.github.io/rvinecopulib/reference/getters.md)
  [`get_parameters()`](https://vinecopulib.github.io/rvinecopulib/reference/getters.md)
  [`get_ktau()`](https://vinecopulib.github.io/rvinecopulib/reference/getters.md)
  [`get_family()`](https://vinecopulib.github.io/rvinecopulib/reference/getters.md)
  [`get_all_pair_copulas()`](https://vinecopulib.github.io/rvinecopulib/reference/getters.md)
  [`get_all_parameters()`](https://vinecopulib.github.io/rvinecopulib/reference/getters.md)
  [`get_all_ktaus()`](https://vinecopulib.github.io/rvinecopulib/reference/getters.md)
  [`get_all_families()`](https://vinecopulib.github.io/rvinecopulib/reference/getters.md)
  :

  Extracts components of `bicop_dist` and `vinecop_dist` objects

## Diagnostics and utilities

- [`pseudo_obs()`](https://vinecopulib.github.io/rvinecopulib/reference/pseudo_obs.md)
  : Pseudo-Observations

- [`emp_cdf()`](https://vinecopulib.github.io/rvinecopulib/reference/emp_cdf.md)
  : Corrected Empirical CDF

- [`pairs_copula_data()`](https://vinecopulib.github.io/rvinecopulib/reference/pairs_copula_data.md)
  : Exploratory pairs plot for copula data

- [`plot(`*`<vinecop_dist>`*`)`](https://vinecopulib.github.io/rvinecopulib/reference/plot.vinecop_dist.md)
  [`plot(`*`<vinecop>`*`)`](https://vinecopulib.github.io/rvinecopulib/reference/plot.vinecop_dist.md)
  [`contour(`*`<vinecop_dist>`*`)`](https://vinecopulib.github.io/rvinecopulib/reference/plot.vinecop_dist.md)
  [`contour(`*`<vinecop>`*`)`](https://vinecopulib.github.io/rvinecopulib/reference/plot.vinecop_dist.md)
  :

  Plotting `vinecop_dist` and `vinecop` objects.

## Package

- [`rvinecopulib-package`](https://vinecopulib.github.io/rvinecopulib/reference/rvinecopulib.md)
  [`rvinecopulib`](https://vinecopulib.github.io/rvinecopulib/reference/rvinecopulib.md)
  : High Performance Algorithms for Vine Copula Modeling
