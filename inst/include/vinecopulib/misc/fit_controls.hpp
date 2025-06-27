// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <Eigen/Dense>
#include <string>
#include <vector>
#include <vinecopulib/misc/tools_optional.hpp>

namespace vinecopulib {

//! Configuration options for initializing a FitControlsVinecop object.
//!
//! This struct provides a flexible way to initialize a `FitControlsVinecop`
//! object. Each field is optional, and default values are applied if the field
//! is not set.
struct FitControlsConfig
{
  //! The set of bicop families to consider. Defaults to all families.
  optional::optional<std::vector<BicopFamily>> family_set;

  //! Method for parametric estimation (e.g., "mle"). Default: "mle".
  optional::optional<std::string> parametric_method;

  //! Method for nonparametric estimation (e.g., "constant"). Default:
  //! "constant".
  optional::optional<std::string> nonparametric_method;

  //! Method for nonparametric estimation (e.g., "constant"). Default:
  //! "constant".
  optional::optional<double> nonparametric_mult;

  //! Criterion for model selection (e.g., "aic"). Default: "aic".
  optional::optional<std::string> selection_criterion;

  //! Observation weights. Default: an empty Eigen::VectorXd.
  optional::optional<Eigen::VectorXd> weights;

  //! Threshold for the psi0 parameter. Default: 0.9.
  optional::optional<double> psi0;

  //! Whether to preselect families based on preliminary criteria. Default:
  //! true.
  optional::optional<bool> preselect_families;

  //! Whether to allow rotations for the families. Default: true.
  optional::optional<bool> allow_rotations;

  //! Number of threads to use during fitting. Default: 1.
  optional::optional<size_t> num_threads;

  //! Truncation level for truncated vines. Default: no truncation.
  optional::optional<size_t> trunc_lvl;

  //! The criterion for selecting the maximum spanning tree. Default: "tau".
  optional::optional<std::string> tree_criterion;

  //! Threshold for thresholded vines. Default: 0.
  optional::optional<double> threshold;

  //! Whether to select the truncation level automatically. Default: false.
  optional::optional<bool> select_threshold;

  //! Whether to select the threshold automatically. Default: false.
  optional::optional<bool> select_trunc_lvl;

  //! Whether to select the families automatically. Default: true.
  optional::optional<bool> select_families;

  //! Whether to show a trace of the building progress. Default: false.
  optional::optional<bool> show_trace;

  //! The algorithm for building the maximum spanning tree. Default: "mst_prim".
  optional::optional<std::string> tree_algorithm;

  //! A vector of random seeds for the random number generator
  optional::optional<std::vector<int>> seeds;
};

}