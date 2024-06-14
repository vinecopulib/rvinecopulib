// Copyright © 2016-2023 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/bicop/abstract.hpp>

namespace vinecopulib {

namespace tools_interpolation {
class InterpolationGrid;
}

class BoxCovering {
public:
  BoxCovering(const Eigen::MatrixXd& u, size_t K = 40) : u_(u), K_(K) {
    boxes_.resize(K);
    for (size_t k = 0; k < K; k++) {
      boxes_[k].resize(K);
      for (size_t j = 0; j < K; j++) {
        boxes_[k][j] = std::make_unique<Box>(
          std::vector<double> {static_cast<double>(k) / K, static_cast<double>(j) / K},
          std::vector<double> {static_cast<double>(k + 1) / K, static_cast<double>(j + 1) / K}
        );
      }
    }

    size_t k, j;
    n_ = u.rows();
    for (size_t i = 0; i < n_; i++) {
      k = std::floor(u(i, 0) * K);
      j = std::floor(u(i, 1) * K);
      boxes_[k][j]->indices_.insert(i);
    }
  }

  std::vector<size_t> get_box_indices(const Eigen::VectorXd& lower,
                                      const Eigen::VectorXd& upper)
  {
    std::vector<size_t> indices;
    indices.reserve(n_);
    auto l0 = std::floor(lower(0) * K_);
    auto l1 = std::floor(lower(1) * K_);
    auto u0 = std::ceil(upper(0) * K_);
    auto u1 = std::ceil(upper(1) * K_);

    for (size_t k = l0; k < u0; k++) {
      for (size_t j = l1; j < u1; j++) {
        for (auto& i : boxes_[k][j]->indices_) {
          if (k == l0 | k == u0 - 1) {
            if (u_(i, 0) < lower(0) | u_(i, 0) > upper(0))
              continue;
          }
          if (j == l1 | j == u1 - 1) {
            if (u_(i, 1) < lower(1) | u_(i, 1) > upper(1))
              continue;
          }
          indices.push_back(i);
        }
      }
    }

    return indices;
  }

  void swap_sample(size_t i, const Eigen::VectorXd& new_sample)
  {
    auto k = std::floor(u_(i, 0) * K_);
    auto j = std::floor(u_(i, 1) * K_);
    boxes_[k][j]->indices_.erase(i);

    u_.row(i) = new_sample;
    k = std::floor(new_sample(0) * K_);
    j = std::floor(new_sample(1) * K_);
    boxes_[k][j]->indices_.insert(i);
  }

private:
  struct Box {
  public:
    Box(std::vector<double> lower, std::vector<double> upper) :
    lower_(lower), upper_(upper)
    {}

    std::vector<double> lower_;
    std::vector<double> upper_;
    std::set<size_t> indices_;
  };

  size_t n_;
  size_t K_;
  Eigen::MatrixXd u_;
  std::vector<std::vector<std::unique_ptr<Box>>> boxes_;

};

Eigen::MatrixXd find_latent_sample(const Eigen::MatrixXd& u, double b, size_t niter = 3,
                                   const std::vector<int>& seeds = std::vector<int>());


//! @brief An abstract class for kernel copulas.
//!
//! Evaluation functions of kernel estimators are implemented efficiently
//! using spline interpolation, see Nagler (2016).
//!
//! This class is used in the implementation underlying the Bicop class.
//! Users should not use AbstractBicop or derived classes directly, but
//! always work with the Bicop interface.
//!
//! @literature
//! Nagler, Thomas. *kdecopula: An R Package for the Kernel Estimation of
//! Copula Densities*. arXiv:1603.04229 [stat.CO], 2016
class KernelBicop : public AbstractBicop
{
public:
  KernelBicop();

protected:
  Eigen::VectorXd pdf_raw(const Eigen::MatrixXd& u) override;

  Eigen::VectorXd pdf(const Eigen::MatrixXd& u) override;

  Eigen::VectorXd cdf(const Eigen::MatrixXd& u) override;

  Eigen::VectorXd hfunc1_raw(const Eigen::MatrixXd& u) override;

  Eigen::VectorXd hfunc2_raw(const Eigen::MatrixXd& u) override;

  Eigen::VectorXd hfunc1(const Eigen::MatrixXd& u) override;

  Eigen::VectorXd hfunc2(const Eigen::MatrixXd& u) override;

  Eigen::VectorXd hinv1_raw(const Eigen::MatrixXd& u) override;

  Eigen::VectorXd hinv2_raw(const Eigen::MatrixXd& u) override;

  double get_npars() const override;

  void set_npars(const double& npars) override;

  Eigen::MatrixXd get_parameters() const override;

  Eigen::MatrixXd get_parameters_lower_bounds() const override;

  Eigen::MatrixXd get_parameters_upper_bounds() const override;

  void set_parameters(const Eigen::MatrixXd& parameters) override;

  double parameters_to_tau(const Eigen::MatrixXd& parameters) override;

  void flip() override;

  Eigen::MatrixXd tau_to_parameters(const double& tau) override;

  Eigen::VectorXd make_normal_grid(size_t m = 30);


  std::shared_ptr<tools_interpolation::InterpolationGrid> interp_grid_;
  double npars_;
};


}

#include <vinecopulib/bicop/implementation/kernel.ipp>
