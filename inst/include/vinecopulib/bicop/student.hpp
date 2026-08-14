// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/bicop/elliptical.hpp>

namespace vinecopulib {
//! @brief The Student t copula.
//!
//! This class is used in the implementation underlying the Bicop class.
//! Users should not use AbstractBicop or derived classes directly, but
//! always work with the Bicop interface.
//!
//! @literature
//! Joe, Harry. Dependence modeling with copulas. CRC Press, 2014.
class StudentBicop : public EllipticalBicop
{
public:
  // constructor
  StudentBicop();

private:
  // evaluation leaves (`parameters` is m x 2, m in {1, n}); these loop per row
  // because the t-distribution helpers take scalar parameters
  Eigen::VectorXd pdf_raw(const Eigen::MatrixXd& u,
                          const Eigen::MatrixXd& parameters) override;

  Eigen::VectorXd cdf(const Eigen::MatrixXd& u,
                      const Eigen::MatrixXd& parameters) override;

  Eigen::VectorXd hfunc1_raw(const Eigen::MatrixXd& u,
                             const Eigen::MatrixXd& parameters) override;

  Eigen::VectorXd hinv1_raw(const Eigen::MatrixXd& u,
                            const Eigen::MatrixXd& parameters) override;

  // analytic derivative leaves (canonical selectors; see tools_deriv),
  // ported from the VineCopula C sources (tcopuladeriv.c,
  // tcopuladeriv_new.c, logderiv.c)
  Eigen::VectorXd pdf_deriv_raw(const Eigen::MatrixXd& u,
                                const Eigen::MatrixXd& parameters,
                                const std::string& deriv) override;

  Eigen::VectorXd pdf_deriv2_raw(const Eigen::MatrixXd& u,
                                 const Eigen::MatrixXd& parameters,
                                 const std::string& deriv) override;

  Eigen::VectorXd hfunc1_deriv_raw(const Eigen::MatrixXd& u,
                                   const Eigen::MatrixXd& parameters,
                                   const std::string& deriv) override;

  Eigen::VectorXd hfunc1_deriv2_raw(const Eigen::MatrixXd& u,
                                    const Eigen::MatrixXd& parameters,
                                    const std::string& deriv) override;

  Eigen::VectorXd logpdf_deriv_raw(const Eigen::MatrixXd& u,
                                   const Eigen::MatrixXd& parameters,
                                   const std::string& deriv) override;

  Eigen::VectorXd logpdf_deriv2_raw(const Eigen::MatrixXd& u,
                                    const Eigen::MatrixXd& parameters,
                                    const std::string& deriv) override;

  // applies a scalar kernel to each row of `u` with per-row parameters
  static Eigen::VectorXd apply_kernel(
    const Eigen::MatrixXd& u,
    const Eigen::MatrixXd& parameters,
    double (*kernel)(double, double, double, double));

  // scalar derivative kernels backing the leaves above; arguments follow
  // vinecopulib conventions (`hfunc1` conditions on `u1`, so the
  // VineCopula h-kernels are called with their `u := u2`, `v := u1`)
  static double pdf_scalar(double u1, double u2, double rho, double nu);
  static double diff_lpdf_rho(double u1, double u2, double rho, double nu);
  static double diff_lpdf_nu(double u1, double u2, double rho, double nu);
  static double diff2_lpdf_rho(double u1, double u2, double rho, double nu);
  static double diff2_lpdf_nu(double u1, double u2, double rho, double nu);
  static double diff2_lpdf_rho_nu(double u1, double u2, double rho, double nu);
  static double diff_pdf_rho(double u1, double u2, double rho, double nu);
  static double diff_pdf_nu(double u1, double u2, double rho, double nu);
  static double diff_pdf_u1(double u1, double u2, double rho, double nu);
  static double diff2_pdf_rho(double u1, double u2, double rho, double nu);
  static double diff2_pdf_nu(double u1, double u2, double rho, double nu);
  static double diff2_pdf_rho_nu(double u1, double u2, double rho, double nu);
  static double diff2_pdf_rho_u1(double u1, double u2, double rho, double nu);
  static double diff2_pdf_nu_u1(double u1, double u2, double rho, double nu);
  static double diff2_pdf_u1(double u1, double u2, double rho, double nu);
  static double diff2_pdf_u1_u2(double u1, double u2, double rho, double nu);
  static double diff_hfunc1_rho(double u1, double u2, double rho, double nu);
  static double diff_hfunc1_nu(double u1, double u2, double rho, double nu);
  static double diff_hfunc1_u1(double u1, double u2, double rho, double nu);
  static double diff2_hfunc1_rho(double u1, double u2, double rho, double nu);
  static double diff2_hfunc1_nu(double u1, double u2, double rho, double nu);
  static double diff2_hfunc1_rho_nu(double u1,
                                    double u2,
                                    double rho,
                                    double nu);
  static double diff2_hfunc1_rho_u1(double u1,
                                    double u2,
                                    double rho,
                                    double nu);
  static double diff2_hfunc1_nu_u1(double u1, double u2, double rho, double nu);
  static double diff2_hfunc1_u1(double u1, double u2, double rho, double nu);

  // single source of truth for the math (shared by the state-based and
  // parameter-aware leaves)
  static Eigen::VectorXd pdf_impl(const Eigen::MatrixXd& u,
                                  double rho,
                                  double nu);
  static Eigen::VectorXd cdf_impl(const Eigen::MatrixXd& u,
                                  double rho,
                                  double nu);
  static Eigen::VectorXd hfunc1_impl(const Eigen::MatrixXd& u,
                                     double rho,
                                     double nu);
  static Eigen::VectorXd hinv1_impl(const Eigen::MatrixXd& u,
                                    double rho,
                                    double nu);

  Eigen::MatrixXd tau_to_parameters(const double& tau) override;

  Eigen::MatrixXd parameters_to_taildep(
    const Eigen::MatrixXd& parameters) override;

  Eigen::VectorXd get_start_parameters(const double tau) override;
};
}

#include <vinecopulib/bicop/implementation/student.ipp>
