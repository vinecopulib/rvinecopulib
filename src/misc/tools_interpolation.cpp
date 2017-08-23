// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_interpolation.hpp>
#include <vinecopulib/misc/tools_stats.hpp>
#include <stdexcept>
#include <iostream>

namespace vinecopulib {

namespace tools_interpolation {
    //! Constructor
    //!
    //! @param grid_points an ascending sequence of grid_points; used in both
    //! dimensions.
    //! @param values a dxd matrix of copula density values evaluated at
    //! (grid_points_i, grid_points_j).
    //! @param norm_times how many times the normalization routine should run.
    InterpolationGrid::InterpolationGrid(const Eigen::VectorXd& grid_points, 
                                         const Eigen::MatrixXd& values,
                                         int norm_times)
    {
        if (values.cols() != values.rows()) {
            throw std::runtime_error("values must be a quadratic matrix");
        }
        if (grid_points.size() != values.rows()) {
            throw std::runtime_error(
                    "number of grid_points must equal dimension of values");
        }

        grid_points_ = grid_points;
        values_ = values;
        normalize_margins(norm_times);
    }

    Eigen::MatrixXd InterpolationGrid::get_values() const
    {
        return values_;
    }

    void InterpolationGrid::set_values(const Eigen::MatrixXd& values,
                                       int norm_times)
    {
        if (values.size() != values_.size()) {
            if (values.rows() != values_.rows()) {
                std::stringstream message;
                message <<
                        "values have has wrong number of rows; " <<
                        "expected: " << values_.rows() << ", " <<
                        "actual: " << values.rows() << std::endl;
                throw std::runtime_error(message.str().c_str());
            }
            if (values.cols() != values_.cols()) {
                std::stringstream message;
                message <<
                        "values have wrong number of columns; " <<
                        "expected: " << values_.cols() << ", " <<
                        "actual: " << values.cols() << std::endl;
                throw std::runtime_error(message.str().c_str());
            }
        }

        values_ = values;
        normalize_margins(norm_times);
    }

    void InterpolationGrid::flip()
    {
        values_.transposeInPlace();
    }
    
    //! renormalizes the estimate to uniform margins
    //! 
    //! @param times how many times the normalization routine should run.
    void InterpolationGrid::normalize_margins(int times)
    {
        size_t m = grid_points_.size();        
        for(int k = 0; k < times; ++k) {
            for (size_t i = 0; i < m; ++i) {
                values_.row(i) /= int_on_grid(1.0, values_.row(i), grid_points_);
            }
            for (size_t j = 0; j < m; ++j) {
                values_.col(j) /= int_on_grid(1.0, values_.col(j), grid_points_);
            }
        }
    }

    //! Interpolation in two dimensions
    //!
    //! @param x mx2 matrix of evaluation points.
    Eigen::VectorXd InterpolationGrid::interpolate(const Eigen::MatrixXd& x)
    {
        ptrdiff_t N = x.rows();
        ptrdiff_t m = grid_points_.size();
        Eigen::VectorXd y(4), out(N), a(4), tmpgrid(4), tmpvals(4);
        ptrdiff_t i = 0;
        ptrdiff_t j = 0;
        ptrdiff_t i0, i3;

        for (ptrdiff_t n = 0; n < N; ++n) {
            // find cell
            bool found_i = false;
            bool found_j = false;
            for (ptrdiff_t k = 1; k < (m-1); ++k) {
                if ((x(n, 0) >= grid_points_(k))) { 
                    i = k;
                } else {
                    found_i = true;
                }
                if ((x(n, 1) >= grid_points_(k))) {
                    j = k;
                } else {
                    found_j = true;
                }
                if (found_i & found_j) {
                    break;
                }
            }

            // construct grid for first direction
            i0 = std::max(i-1, (ptrdiff_t) 0);
            i3 = std::min(i+2, m-1);
            tmpgrid(0) = grid_points_(i0);
            tmpgrid(1) = grid_points_(i);
            tmpgrid(2) = grid_points_(i+1);
            tmpgrid(3) = grid_points_(i3);

            // interpolate in one direction (four times)
            for (ptrdiff_t s = 0; s < 4; ++s) {
                i0 = std::max(i-1, (ptrdiff_t) 0);
                i3 = std::min(i+2, m-1);
                ptrdiff_t jj = std::min(m-1, j-1+s);
                jj = std::max((ptrdiff_t) 0, jj);

                tmpvals(0) = values_(i0,  jj);
                tmpvals(1) = values_(i,   jj);
                tmpvals(2) = values_(i+1, jj);
                tmpvals(3) = values_(i3,  jj);

                y(s) = interp_on_grid(x(n, 0), tmpvals, tmpgrid);
                y(s) = fmax(y(s), 0.0);
            }

            // use these four points to interpolate in the remaining direction#
            i0 = std::max(j-1, (ptrdiff_t) 0);
            i3 = std::min(j+2, m-1);
            tmpgrid(0) = grid_points_(i0);
            tmpgrid(1) = grid_points_(j);
            tmpgrid(2) = grid_points_(j+1);
            tmpgrid(3) = grid_points_(i3);
            
            out(n) = interp_on_grid(x(n, 1), y, tmpgrid);
            out(n) = fmax(out(n), 1e-15);
        }

        return out;
    }

    //! Integrate the grid along one axis
    //!
    //! @param u mx2 matrix of evaluation points
    //! @param cond_var either 1 or 2; the axis considered fixed.
    //!
    Eigen::VectorXd InterpolationGrid::intergrate_1d(const Eigen::MatrixXd& u, 
                                                     size_t cond_var)
    {
        ptrdiff_t n = u.rows();
        ptrdiff_t m = grid_points_.size();
        Eigen::VectorXd tmpvals(m), out(n), tmpa(4), tmpb(4);
        Eigen::MatrixXd tmpgrid(m, 2);
        double upr = 0.0;
        double tmpint, int1;
        tmpint = 0.0;
        
        for (ptrdiff_t i = 0; i < n; ++i) {
            if (cond_var == 1) {
                upr = u(i, 1);
                tmpgrid.col(0) = Eigen::VectorXd::Constant(m, u(i, 0));
                tmpgrid.col(1) = grid_points_;
            } else if (cond_var == 2) {
                upr = u(i, 0);
                tmpgrid.col(0) = grid_points_;
                tmpgrid.col(1) = Eigen::VectorXd::Constant(m, u(i, 1));
            }
            tmpvals = interpolate(tmpgrid);
            tmpint = int_on_grid(upr, tmpvals, grid_points_);
            int1 = int_on_grid(1.0, tmpvals, grid_points_);
            out(i) = tmpint/int1;
            out(i) = fmax(out(i), 1e-10);
            out(i) = fmin(out(i), 1-1e-10);
        }

        return out;
    }

    //! Integrate the grid along the two axis
    //!
    //! @param u mx2 matrix of evaluation points
    //!
    Eigen::VectorXd InterpolationGrid::intergrate_2d(const Eigen::MatrixXd& u)
    {

        double upr, tmpint, tmpint1;
        ptrdiff_t n = u.rows();
        ptrdiff_t m = grid_points_.size();
        Eigen::VectorXd tmpvals(m), tmpvals2(m), out(n);
        Eigen::MatrixXd tmpgrid(m, 2);
        tmpgrid.col(1) = grid_points_;

        for (ptrdiff_t i = 0; i < n; ++i) {
            upr = u(i, 1);
            for (ptrdiff_t k = 0; k < m-1; ++k) {
                tmpgrid.col(0) = Eigen::VectorXd::Constant(m,  grid_points_(k));
                tmpvals = interpolate(tmpgrid);
                tmpint = int_on_grid(upr, tmpvals, grid_points_);
                tmpvals2(i) = tmpint;
            }
            upr = u(i, 0);
            tmpint = int_on_grid(upr, tmpvals2, grid_points_);
            tmpint1 = int_on_grid(1.0, tmpvals2, grid_points_);
            out(i) = tmpint/tmpint1;
            out(i) = fmax(out(i), 1e-10);
            out(i) = fmin(out(i), 1-1e-10);
        }

        return out;
    }

// ---------------- Utility functions for spline interpolation ----------------

    //! Evaluate a cubic polynomial
    //!
    //! @param x evaluation point.
    //! @param a polynomial coefficients
    double InterpolationGrid::cubic_poly(const double& x, const Eigen::VectorXd& a)
    {
        double x2 = x*x;
        double x3 = x2*x;
        return a(0) + a(1)*x + a(2)*x2 + a(3)*x3;
    }

    //! Indefinite integral of a cubic polynomial
    //!
    //! @param x evaluation point.
    //! @param a polynomial coefficients.
    double InterpolationGrid::cubic_indef_integral(const double& x,
                                                   const Eigen::VectorXd& a)
    {
        double x2 = x*x;
        double x3 = x2*x;
        double x4 = x3*x;
        return a(0)*x + a(1)/2.0*x2 + a(2)/3.0*x3 + a(3)/4.0*x4;
    }

    //! Definite integral of a cubic polynomial
    //!
    //! @param lower lower limit of the integral.
    //! @param upper upper limit of the integral.
    //! @param a polynomial coefficients.
    double InterpolationGrid::cubic_integral(const double& lower, 
                                             const double& upper, 
                                             const Eigen::VectorXd& a)
    {
        return cubic_indef_integral(upper, a) - cubic_indef_integral(lower, a);
    }

    //! Calculate coefficients for cubic intrpolation spline
    //!
    //! @param vals length 4 vector of function values.
    //! @param grid length 4 vector of grid points.
    Eigen::VectorXd InterpolationGrid::find_coefs(const Eigen::VectorXd& vals, 
                                                  const Eigen::VectorXd& grid)
    {
        Eigen::VectorXd a(4);

        double dt0 = grid(1) - grid(0);
        double dt1 = grid(2) - grid(1);
        double dt2 = grid(3) - grid(2);

        /* check for repeated points (important for boundaries) */
        if (dt1 < 1e-4) dt1 = 1.0;
        if (dt0 < 1e-4) dt0 = dt1;
        if (dt2 < 1e-4) dt2 = dt1;

        // compute tangents when parameterized in (t1,t2)
        double dx1 = (vals(1) - vals(0)) / dt0;
        dx1 -= (vals(2) - vals(0)) / (dt0 + dt1);
        dx1 += (vals(2) - vals(1)) / dt1;
        double dx2 = (vals(2) - vals(1)) / dt1;
        dx2 -= (vals(3) - vals(1)) / (dt1 + dt2);
        dx2 += (vals(3) - vals(2)) / dt2;

        // rescale tangents for parametrization in (0,1)
        dx1 *= dt1;
        dx2 *= dt1;

        // compute coefficents
        a(0) = vals(1);
        a(1) = dx1;
        a(2) = -3*vals(1) + 3*vals(2) - 2*dx1 - dx2;
        a(3) = 2*vals(1) - 2*vals(2) + dx1 + dx2;

        return a;
    }

    //! Interpolate on 4 points
    //!
    //! @param x evaluation point.
    //! @param vals length 4 vector of function values.
    //! @param grid length 4 vector of grid points.
    double InterpolationGrid::interp_on_grid(const double& x, 
                                             const Eigen::VectorXd& vals, 
                                             const Eigen::VectorXd& grid)
    {
        Eigen::VectorXd a = find_coefs(vals, grid);
        double xev = fmax((x - grid(1)), 0) / (grid(2) - grid(1));
        return cubic_poly(xev, a);
    }


// ---------------- Utility functions for integration ----------------


    //! Integrate a spline interpolant
    //!
    //! @param upr upper limit of integration (lower is 0).
    //! @param vals vector of values to be interpolated and integrated.
    //! @param grid vector of grid points on which vals has been computed.
    //!
    //! @return Integral of interpolation spline defined by (vals, grid).
    double InterpolationGrid::int_on_grid(const double& upr,
                                          const Eigen::VectorXd& vals, 
                                          const Eigen::VectorXd& grid)
    {
        ptrdiff_t m = grid.size();
        Eigen::VectorXd tmpvals(4), tmpgrid(4), tmpa(4), a(4);
        double uprnew, newint;

        double tmpint = 0.0;

        if (upr > grid(0)) {
            // go up the grid and integrate
            for (ptrdiff_t k = 0; k < m-1; ++k) {
                // stop loop if fully integrated
                if (upr < grid(k)) break;

                // select length 4 subvectors and calculate spline coefficients
                tmpvals(0) = vals(std::max(k-1, (ptrdiff_t) 0));
                tmpvals(1) = vals(k);
                tmpvals(2) = vals(k+1);
                tmpvals(3) = vals(std::min(k+2, m-1));

                tmpgrid(0) = grid(std::max(k-1, (ptrdiff_t) 0));
                tmpgrid(1) = grid(k);
                tmpgrid(2) = grid(k+1);
                tmpgrid(3) = grid(std::min(k+2, m-1));

                tmpa = find_coefs(tmpvals, tmpgrid);

                // don't integrate over full cell if upr is in interior
                uprnew = (upr - grid(k)) / (grid(k+1) - grid(k));
                newint = cubic_integral(0.0, fmin(1.0, uprnew), tmpa);
                tmpint += newint * (grid(k+1) - grid(k));
            }
        }

        return tmpint;
    }
}

}
