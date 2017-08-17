// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/bicop/tll0.hpp>
#include <vinecopulib/bicop/family.hpp>
#include <vinecopulib/misc/tools_stats.hpp>

namespace vinecopulib
{
    Tll0Bicop::Tll0Bicop()
    {
        family_ = BicopFamily::tll0;
    }

    Eigen::VectorXd Tll0Bicop::gaussian_kernel_2d(
        const Eigen::Matrix<double, Eigen::Dynamic, 2>& x
    )
    {
        return tools_stats::dnorm(x).rowwise().prod();
    }

    void Tll0Bicop::fit(const Eigen::Matrix<double, Eigen::Dynamic, 2>& data, 
                        std::string, double mult)
    {
        // construct default grid (equally spaced on Gaussian scale)
        size_t m = 30;
        Eigen::VectorXd grid_points(m);
        for (size_t i = 0; i < m; ++i)
            grid_points(i) = - 3.25 + i * (6.25 / (double) m);
        grid_points = tools_stats::pnorm(grid_points);

        // expand the interpolation grid; a matrix with two columns where each row
        // contains one combination of the grid points
        auto grid_2d = tools_eigen::expand_grid(grid_points);

        // transform evaluation grid and data by inverse Gaussian cdf
        Eigen::Matrix<double, Eigen::Dynamic, 2> z = tools_stats::qnorm(grid_2d);
        Eigen::Matrix<double, Eigen::Dynamic, 2> z_data = tools_stats::qnorm(data);

        // apply normal density to z (used later for normalization)
        Eigen::Matrix<double, Eigen::Dynamic, 2> phi = tools_stats::dnorm(z);
        Eigen::Matrix<double, Eigen::Dynamic, 2> phi_data = tools_stats::dnorm(z_data);

        // find bandwitools_stats::dth matrix
        size_t n = data.rows();
        Eigen::Matrix<double, Eigen::Dynamic, 2> centered =
            z_data.rowwise() - z_data.colwise().mean();
        Eigen::Matrix2d cov = (centered.adjoint() * centered) / double(n - 1);

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> takes_root(cov);
        Eigen::Matrix2d cov_root = takes_root.operatorSqrt();
        Eigen::Matrix2d B = 1.25 * std::pow(n, - 1.0 / 6.0) * cov_root.transpose();
        B *= mult;

        // apply bandwitools_stats::dth matrix
        z = (B.inverse() * z.transpose()).transpose();
        z_data = (B.inverse() * z_data.transpose()).transpose();

        // compute estimator on each evaluation point
        Eigen::VectorXd kernels(n);
        double det_B = B.determinant();
        size_t i = 0;
        size_t j = 0;
        Eigen::MatrixXd values(m, m);
        for (size_t k = 0; k < m * m; ++k) {
            kernels = gaussian_kernel_2d((z_data - z.row(k).replicate(n, 1)));
            values(j, i) = kernels.mean() / (det_B * phi.row(k).prod());
            ++i;
            if (i % m == 0) {
                i = 0;
                ++j;
            }
        }

        // for interpolation, we shift the limiting gridpoints to 0 and 1
        grid_points(0) = 0.0;
        grid_points(m - 1) = 1.0;
        interp_grid_ = tools_interpolation::InterpolationGrid(grid_points, values);

        // compute effective number of parameters
        double K0 = gaussian_kernel_2d(Eigen::MatrixXd::Constant(1, 2, 0.0))(0);
        Eigen::VectorXd scale = phi_data.rowwise().prod();
        npars_ =  K0 / det_B / (scale.array() * this->pdf(data).array()).mean();
    }
}
