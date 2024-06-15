// Copyright © 2016-2023 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_interpolation.hpp>
#include <vinecopulib/misc/tools_stats.hpp>
#include <vinecopulib/misc/quadtree.hpp>
#include <wdm/eigen.hpp>

namespace vinecopulib {

 inline
 Eigen::MatrixXd find_latent_sample(
  const Eigen::MatrixXd& u, double b, size_t niter,
  const std::vector<int>& seeds)
 {
   using namespace tools_stats;

   // for better cache behavior
   using MatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

   size_t n = u.rows();
   auto seed_state = SeedState(seeds);

   auto w = simulate_uniform(n, 2, false, seed_state.next());
   MatrixXd uu = w.array() * u.leftCols(2).array() +
     (1 - w.array()) * u.rightCols(2).array();

   MatrixXd x = qnorm(uu);
   auto depth = std::ceil(std::log2(10 / b));

   MatrixXd lb = safe_qnorm(u.rightCols(2)).array() - b;
   MatrixXd ub = safe_qnorm(u.leftCols(2)).array() + b;

   MatrixXd norm_sim(n, 2);
   Point new_sample, old_sample;

   QuadTree quadtree(BoundingBox(-4.5, -4.5, 4.5, 4.5), depth, seed_state.next());
   for (size_t i = 0; i < n; i++) {
     quadtree.insert(Point{x(i, 0), x(i, 1), i});
   }

   for (size_t it = 0; it < niter; it++) {

     norm_sim = simulate_normal(n, 2, seed_state.next()).array() * b;

     for (size_t i = 0; i < n; i++) {
       try {
        new_sample = quadtree.sample(
          BoundingBox(lb(i, 0), lb(i, 1), ub(i, 0), ub(i, 1))
        );
        quadtree.remove(Point{x(i, 0), x(i, 1), i});

        new_sample.x += norm_sim(i, 0);
        new_sample.y += norm_sim(i, 1);
        new_sample.index = i;
        x(i, 0) = new_sample.x;
        x(i, 1) = new_sample.y;
        quadtree.insert(new_sample);
       } catch (std::runtime_error& e) {
         // if we can't find a new sample, we just keep the old one
       }
     }
   }

   return to_pseudo_obs(x);
 }


inline KernelBicop::KernelBicop()
{
  // construct default grid (equally spaced on Gaussian scale)
  size_t m = 30;
  auto grid_points = this->make_normal_grid(m);

  // move boundary points to 0/1, so we don't have to extrapolate
  grid_points(0) = 0.0;
  grid_points(m - 1) = 1.0;

  interp_grid_ = std::make_shared<tools_interpolation::InterpolationGrid>(
    grid_points, Eigen::MatrixXd::Constant(m, m, 1.0) // independence
  );
  npars_ = 0.0;
}

inline Eigen::VectorXd
KernelBicop::pdf_raw(const Eigen::MatrixXd& u)
{
  auto pdf = interp_grid_->interpolate(u);
  tools_eigen::trim(pdf, 1e-20, DBL_MAX);
  return pdf;
}

inline Eigen::VectorXd
KernelBicop::pdf(const Eigen::MatrixXd& u)
{
  if (u.cols() == 4) {
    // evaluate jittered density at mid rank for stability
    return pdf_raw((u.leftCols(2) + u.rightCols(2)).array() / 2.0);
  }
  return pdf_raw(u);
}

inline Eigen::VectorXd
KernelBicop::cdf(const Eigen::MatrixXd& u)
{
  return interp_grid_->integrate_2d(u);
}

inline Eigen::VectorXd
KernelBicop::hfunc1_raw(const Eigen::MatrixXd& u)
{
  return interp_grid_->integrate_1d(u, 1);
}

inline Eigen::VectorXd
KernelBicop::hfunc2_raw(const Eigen::MatrixXd& u)
{
  return interp_grid_->integrate_1d(u, 2);
}

inline Eigen::VectorXd
KernelBicop::hfunc1(const Eigen::MatrixXd& u)
{
  if (u.cols() == 4) {
    auto u_avg = u;
    u_avg.col(0) = (u.col(0) + u.col(2)).array() / 2.0;
    return hfunc1_raw(u_avg.leftCols(2));
  }
  return hfunc1_raw(u);
}

inline Eigen::VectorXd
KernelBicop::hfunc2(const Eigen::MatrixXd& u)
{
  if (u.cols() == 4) {
    auto u_avg = u;
    u_avg.col(1) = (u.col(1) + u.col(3)).array() / 2.0;
    return hfunc2_raw(u_avg.leftCols(2));
  }
  return hfunc2_raw(u);
}

inline Eigen::VectorXd
KernelBicop::hinv1_raw(const Eigen::MatrixXd& u)
{
  return hinv1_num(u);
}

inline Eigen::VectorXd
KernelBicop::hinv2_raw(const Eigen::MatrixXd& u)
{
  return hinv2_num(u);
}

inline double
KernelBicop::parameters_to_tau(const Eigen::MatrixXd& parameters)
{
  auto oldpars = this->get_parameters();
  auto old_types = var_types_;
  this->set_parameters(parameters);
  var_types_ = { "c", "c" };

  std::vector<int> seeds = {
    204967043, 733593603, 184618802, 399707801, 290266245
  };
  auto u = tools_stats::ghalton(1000, 2, seeds);
  u.col(1) = hinv1_raw(u);

  this->set_parameters(oldpars);
  var_types_ = old_types;
  return wdm::wdm(u, "tau")(0, 1);
}

inline double
KernelBicop::get_npars() const
{
  return npars_;
}

inline void
KernelBicop::set_npars(const double& npars)
{
  if (npars < 0) {
    throw std::runtime_error("npars must be positive.");
  }
  npars_ = npars;
}

inline Eigen::MatrixXd
KernelBicop::get_parameters() const
{
  return interp_grid_->get_values();
}

inline Eigen::MatrixXd
KernelBicop::get_parameters_lower_bounds() const
{
  return Eigen::MatrixXd::Constant(30, 30, 0.0);
}

inline Eigen::MatrixXd
KernelBicop::get_parameters_upper_bounds() const
{
  return Eigen::MatrixXd::Constant(30, 30, 1e4);
}

inline void
KernelBicop::set_parameters(const Eigen::MatrixXd& parameters)
{
  if (parameters.minCoeff() < 0) {
    std::stringstream message;
    message << "density should be larger than 0. ";
    throw std::runtime_error(message.str().c_str());
  }
  // don't normalize again!
  interp_grid_->set_values(parameters, 0);
}

inline void
KernelBicop::flip()
{
  interp_grid_->flip();
}

inline Eigen::MatrixXd
KernelBicop::tau_to_parameters(const double& tau)
{
  return no_tau_to_parameters(tau);
}

// construct default grid (equally spaced on Gaussian scale)
inline Eigen::VectorXd
KernelBicop::make_normal_grid(size_t m)
{
  Eigen::VectorXd grid_points(m);
  for (size_t i = 0; i < m; ++i)
    grid_points(i) =
      -3.25 + static_cast<double>(i) * (6.5 / static_cast<double>(m - 1));
  grid_points = tools_stats::pnorm(grid_points);

  return grid_points;
}
}
