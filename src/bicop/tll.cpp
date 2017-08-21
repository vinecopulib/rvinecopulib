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
            mult = std::pow(n, -1.0 / 3.0);
        } else {
            double degree;
            if (method == "linear") {
                degree = 1.0;
            } else {
                degree = 2.0;
            }
            mult = std::pow(n, - 1.0 / (2.0 * degree + 1.0));
        }

        return mult * cov;
    }

    Eigen::VectorXd TllBicop::ftll(
            const Eigen::Matrix<double, Eigen::Dynamic, 2>& x,
            const Eigen::Matrix<double, Eigen::Dynamic, 2>& x_data,
            Eigen::Matrix2d B,
            std::string method)
    {
        // compute inverse/root bandwidth matrices
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> takes_root(B);
        Eigen::Matrix2d rB = takes_root.operatorSqrt();
        Eigen::Matrix2d irB = rB.inverse();
        double det_irB = irB.determinant();
        Eigen::Matrix2d iB = B.inverse();
        
        // apply bandwidth matrix
        Eigen::MatrixXd z = (irB * x.transpose()).transpose();
        Eigen::MatrixXd z_data = (irB * x_data.transpose()).transpose();

        size_t n = x.rows();
        size_t m = x_data.rows();
        double f0;
        Eigen::Vector2d f1, b;
        Eigen::MatrixXd f2, S;
        Eigen::VectorXd kernels(m);
        Eigen::MatrixXd zz(m, 2), zz2(m, 2);
        Eigen::VectorXd res = Eigen::VectorXd::Ones(n);
        S = B;
        for (size_t k = 0; k < n; ++k) {
            zz = z_data - z.row(k).replicate(m, 1);
            kernels = gaussian_kernel_2d(zz);
            f0 = kernels.mean() * det_irB;
            res(k) *= f0;
            if (method != "constant") {
                zz = (rB * zz.transpose()).transpose();
                f1(0) = zz.col(0).cwiseProduct(kernels).mean() * det_irB;
                f1(1) = zz.col(1).cwiseProduct(kernels).mean() * det_irB;
                f1 = iB * f1;
                if (method == "linear") {
                    b(0) = f1(0) / f0;
                    b(1) = f1(1) / f0;
                    f1 = iB * f1;
                } 
                // else {
                //     zz2.col(0) = zz.col(0).cwiseProduct(kernels);
                //     zz2.col(1) = zz.col(1).cwiseProduct(kernels);
                //     f2 = zz.transpose() * zz2 * det_irB - iB * f0;
                //     b = B * f1 / f0;
                //     S = ((B * f2 * B) / f0 + B - b * b.transpose()).inverse();
                //     res(k) *= std::sqrt(S.determinant()) / det_irB;
                // }
                res(k) *= std::exp(-0.5 * (b.transpose() * S * b)(0));
            }
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

        // apply normal density to z (used later for normalization)
        Eigen::Matrix<double, Eigen::Dynamic, 2> phi = tools_stats::dnorm(z);
        Eigen::Matrix<double, Eigen::Dynamic, 2> phi_data = tools_stats::dnorm(z_data);

        // find bandwidth matrix
        Eigen::Matrix2d B = bandwidth(z_data, method);
        B *= mult;

        // compute the density estimator on z-scale
        Eigen::VectorXd f = ftll(z, z_data, B, method);
        // bring to copula sacle (divide through std normal marginal densities)
        f = f.array() / phi.rowwise().prod().array();
        // Rcpp::Rcout << phi.array().rowwise().prod() << std::endl;
        Eigen::MatrixXd values(m, m);
        values = Eigen::Map<Eigen::MatrixXd>(f.data(), m, m).transpose();

        // for interpolation, we shift the limiting gridpoints to 0 and 1
        grid_points(0) = 0.0;
        grid_points(m - 1) = 1.0;
        interp_grid_ = tools_interpolation::InterpolationGrid(grid_points, values);

        // compute effective number of parameters
        double K0 = gaussian_kernel_2d(Eigen::MatrixXd::Constant(1, 2, 0.0))(0);
        Eigen::VectorXd scale = phi_data.rowwise().prod();
        npars_ =  K0 / B.determinant() / (scale.array() * this->pdf(data).array()).mean();
    }
}
