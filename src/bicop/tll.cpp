// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/bicop/tll.hpp>
#include <vinecopulib/bicop/family.hpp>
#include <vinecopulib/misc/tools_stats.hpp>

namespace vinecopulib
{
    TllBicop::TllBicop()
    {
        family_ = BicopFamily::tll;
    }

    Eigen::VectorXd TllBicop::gaussian_kernel_2d(
        const Eigen::Matrix<double, Eigen::Dynamic, 2>& x
    )
    {
        return tools_stats::dnorm(x).rowwise().prod();
    }

    Eigen::Matrix2d TllBicop::bandwidth(
            const Eigen::Matrix<double, Eigen::Dynamic, 2>& x,
            std::string method
    )
    {
        size_t n = x.rows();
        Eigen::Matrix<double, Eigen::Dynamic, 2> centered =
                x.rowwise() - x.colwise().mean();
        Eigen::Matrix2d cov = (centered.adjoint() * centered) / double(n - 1);

        double mult;
        if (method == "constant") {
            mult = std::pow(n, - 1.0 / 3.0);
        } else {
            double degree;
            if (method == "linear") {
                degree = 1.0;
            } else {
                degree = 2.0;
            }
            mult = 9.0 * std::pow(n, - 1.0 / (2.0 * degree + 1.0));
        }

        return mult * cov;
    }

    Eigen::Matrix2d chol22(const Eigen::Matrix2d& B) {

        Eigen::Matrix2d rB;

        rB(0,0) = std::sqrt(B(0,0));
        rB(0,1) = 0.0;
        rB(1,0) = B(1,0) / rB(0,0);
        rB(1,1) = std::sqrt( B(1, 1) - rB(1, 0) * rB(1, 0));

        return rB;
    }

    Eigen::VectorXd TllBicop::ftll(
            const Eigen::Matrix<double, Eigen::Dynamic, 2>& x,
            const Eigen::Matrix<double, Eigen::Dynamic, 2>& x_data,
            const Eigen::Matrix2d& B,
            std::string method)
    {
        Eigen::Matrix2d irB = chol22(B).inverse();
        double det_irB = irB.determinant();

        Eigen::MatrixXd z = (irB * x.transpose()).transpose();
        Eigen::MatrixXd z_data = (irB * x_data.transpose()).transpose();

        size_t m = x.rows();
        size_t n = x_data.rows();
        double f0;
        Eigen::Vector2d b;
        Eigen::Matrix2d S(B);
        Eigen::VectorXd kernels(n);
        Eigen::MatrixXd zz(n, 2), zz2(n, 2);
        Eigen::VectorXd res = Eigen::VectorXd::Ones(m);

        for (size_t k = 0; k < m; ++k) {
            zz = z_data - z.row(k).replicate(n, 1);
            kernels = gaussian_kernel_2d(zz);
            f0 = kernels.mean();
            if (method != "constant") {
                zz = (irB * zz.transpose()).transpose();
                b = zz.cwiseProduct(kernels.replicate(1, 2)).colwise().mean() / f0;
                if (method == "quadratic") {
                    zz2 = zz.cwiseProduct(kernels.replicate(1, 2)) /
                          (f0 * (double) n);
                    b = B * b;
                    S = (B * (zz.transpose() * zz2) * B -
                         b * b.transpose()).inverse();
                    res(k) *= std::sqrt(S.determinant()) / det_irB;
                }
                res(k) *= std::exp(- 0.5 * double(b.transpose() * S * b));
            }
            res(k) *= f0 * det_irB;
        }
        return res;
    }

    void TllBicop::fit(const Eigen::Matrix<double, Eigen::Dynamic, 2>& data, 
                        std::string method, double mult)
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

        // find bandwidth matrix
        Eigen::Matrix2d B = bandwidth(z_data, method);
        B *= mult;

        // compute the density estimator
        Eigen::VectorXd f = ftll(z, z_data, B, method);

        // apply normal density to z (used later for normalization)
        f = f.cwiseQuotient(tools_stats::dnorm(z).rowwise().prod());

        // store values in mxm grid
        Eigen::MatrixXd values(m, m);
        values = Eigen::Map<Eigen::MatrixXd>(f.data(), m, m).transpose();
        
        // for interpolation, we shift the limiting gridpoints to 0 and 1
        grid_points(0) = 0.0;
        grid_points(m - 1) = 1.0;
        interp_grid_ = tools_interpolation::InterpolationGrid(grid_points, values);

        // compute effective number of parameters
        double K0 = gaussian_kernel_2d(Eigen::MatrixXd::Constant(1, 2, 0.0))(0);
        Eigen::VectorXd scale = tools_stats::dnorm(z_data).rowwise().prod();
        npars_ =  K0 / B.determinant() / (scale.array() * this->pdf(data).array()).mean();
    }
}
