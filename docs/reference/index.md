# Package index

## Full distributions and margins

<!-- end list -->

  - `vine()` `vine_dist()` : Vine copula models
  - `dvine()` `pvine()` `rvine()` : Vine based distributions
  - `predict(<vine>)` `fitted(<vine>)` : Predictions and fitted values
    for a vine copula model
  - `as_margin()` : Normalize an object to the fitted-margin protocol
  - `margin_dist()` : Create a fitted or fixed custom margin
  - `stats_margin()` : Create a fixed margin from a stats distribution
  - `margin_family()` : Define a custom margin family
  - `kde1d_family()` : Define a kde1d margin family
  - `univariateML_family()` : Define a univariateML margin family
  - `dmargin()` `pmargin()` `qmargin()` `margin_info()` : Fitted
    marginal distribution protocol
  - `fit_margin()` : Margin-family fitting protocol
  - `zero_inflated()` : Declare zero-inflated data

## Vine copula models

<!-- end list -->

  - `vinecop()` : Fitting vine copula models
  - `vinecop_dist()` : Vine copula models
  - `dvinecop()` `scores()` `hessian()` `pvinecop()` `rvinecop()` : Vine
    copula distributions
  - `predict(<vinecop>)` `fitted(<vinecop>)` : Predictions and fitted
    values for a vine copula model
  - `rosenblatt()` `inverse_rosenblatt()` : Rosenblatt and inverse
    Rosenblatt transforms
  - `truncate_model()` : Truncate a vine copula model
  - `mBICV()` : Modified vine copula Bayesian information criterion
    (mBICv)
  - `vcov(<vinecop_dist>)` `vcov(<bicop_dist>)`
    `confint(<vinecop_dist>)` `confint(<bicop_dist>)` : Parameter
    uncertainty

## Bivariate copula models

<!-- end list -->

  - `bicop()` : Fit and select bivariate copula models

  - `bicop_dist()` : Bivariate copula models

  - `dbicop()` `pbicop()` `rbicop()` `scores(<bicop_dist>)`
    `hessian(<bicop_dist>)` `hbicop()` : Bivariate copula distributions

  - `tail_dep()` `blomqvist_beta()` : Dependence measures of a bivariate
    copula

  - `predict(<bicop_dist>)` `fitted(<bicop>)` : Predictions and fitted
    values for a bivariate copula model

  - `plot(<bicop_dist>)` `plot(<bicop>)` `contour(<bicop_dist>)`
    `contour(<bicop>)` :
    
    Plotting tools for `bicop_dist` and `bicop` objects

  - `par_to_ktau()` `ktau_to_par()` : Conversion between Kendall's tau
    and parameters

  - `as.bicop()` : Convert list to bicop object

## R-vine structures

<!-- end list -->

  - `rvine_structure()` `cvine_structure()` `dvine_structure()`
    `rvine_matrix()` : R-vine structure

  - `as_rvine_structure()` `as_rvine_matrix()` : Coerce various kind of
    objects to R-vine structures and matrices

  - `rvine_structure_sim()` `rvine_matrix_sim()` : Simulate R-vine
    structures

  - `plot(<rvine_structure>)` `plot(<rvine_matrix>)` : Plotting R-vine
    structures

  - `get_structure()` `get_pair_copula()` `get_parameters()`
    `get_ktau()` `get_family()` `get_all_pair_copulas()`
    `get_all_parameters()` `get_all_ktaus()` `get_all_families()` :
    
    Extracts components of `bicop_dist` and `vinecop_dist` objects

## Diagnostics and utilities

<!-- end list -->

  - `pseudo_obs()` : Pseudo-Observations

  - `emp_cdf()` : Corrected Empirical CDF

  - `pairs_copula_data()` : Exploratory pairs plot for copula data

  - `plot(<vinecop_dist>)` `plot(<vinecop>)` `contour(<vinecop_dist>)`
    `contour(<vinecop>)` :
    
    Plotting `vinecop_dist` and `vinecop` objects.

## Package

<!-- end list -->

  - `rvinecopulib-package` `rvinecopulib` : High Performance Algorithms
    for Vine Copula Modeling
