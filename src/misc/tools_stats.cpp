// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "misc/tools_stats.hpp"
#include "misc/tools_stl.hpp"
#include "misc/tools_c.h"

//! @file misc/tools_stats.cpp

//! Utilities for statistical analysis
namespace tools_stats {
    
    //! simulates from the multivariate uniform distribution
    //!
    //! @param n number of observations.
    //! @param d dimension.
    //!
    //! @return An \f$ n \times d \f$ matrix of independent 
    //! \f$ \mathrm{U}[0, 1] \f$ random variables.
    Eigen::MatrixXd simulate_uniform(size_t n, size_t d)
    {
        if ((n < 1) | (d < 1)) {
            throw std::runtime_error("both n and d must be at least 1.");        
        }
        Eigen::MatrixXd U(n, d);
        std::random_device rd;
        std::default_random_engine generator(rd());
        std::uniform_real_distribution<double> distribution(0.0, 1.0);
        return U.unaryExpr([&](double) { return distribution(generator); });
    }
    
    //! applies the empirical probability integral transform to a data matrix.
    //!
    //! Gives pseudo-observations from the copula by applying the empirical
    //! distribution function (scaled by n + 1) to each margin/column.
    //!
    //! @param x a matrix of real numbers.
    //! @param ties_method indicates how to treat ties; same as in R, see
    //! https://stat.ethz.ch/R-manual/R-devel/library/base/html/rank.html.
    //! @return Psuedo-observations of the copula, i.e. F_X(X) (column-wise)
    Eigen::MatrixXd to_pseudo_obs(Eigen::MatrixXd x, std::string ties_method)
    {
        for (int j = 0; j < x.cols(); ++j)
            x.col(j) = to_pseudo_obs_1d((Eigen::VectorXd) x.col(j), ties_method);

        return x;
    }
    
    //! applies the empirical probability integral transform to a data vector.
    //!
    //! Gives pseudo-observations from the copula by applying the empirical
    //! distribution function (scaled by n + 1) to each margin/column.
    //!
    //! @param x a vector of real numbers.
    //! @param ties_method indicates how to treat ties; same as in R, see
    //! https://stat.ethz.ch/R-manual/R-devel/library/base/html/rank.html.
    //! @return Psuedo-observations of the copula, i.e. F_X(X) (column-wise)
    Eigen::VectorXd to_pseudo_obs_1d(Eigen::VectorXd x, std::string ties_method)
    {
        size_t n = x.size();
        std::vector<double> xvec(x.data(), x.data() + n);
        auto order = tools_stl::get_order(xvec);
        if (ties_method == "first") {
            for (auto i : order)
                x[order[i]] = (double) (i + 1);
        } else if (ties_method == "average") {
            for (size_t i = 0, reps; i < n; i += reps) {
                // find replications
                reps = 1;
                while ((i + reps < n) && (x[order[i]] == x[order[i + reps]]))
                    ++reps;
                // assign average rank of the tied values
                for (size_t k = 0; k < reps; ++k)
                    x[order[i + k]] = i + 1 + (reps - 1) / 2.0;
            }
        } else if (ties_method == "random") {
            // set up random number generator
            std::random_device rd;
            std::default_random_engine gen(rd());
            auto sim = [&] (int m) {
                std::uniform_int_distribution<> distr(0, m - 1);
                return distr(gen);
            };
            for (size_t i = 0, reps; i < n; i += reps) {
                // find replications
                reps = 1;
                while ((i + reps < n) && (x[order[i]] == x[order[i + reps]]))
                    ++reps;
                // assign random rank between ties
                std::vector<size_t> rvals(reps);
                std::iota(rvals.begin(), rvals.end(), 0);  // 0, 1, 2, ...
                std::random_shuffle(rvals.begin(), rvals.end(), sim);
                for (size_t k = 0; k < reps; ++k)
                    x[order[i + k]] = (double) (i + 1 + rvals[k]);
            }
        } else {
            std::stringstream msg;
            msg << "unknown ties method (" << ties_method << ")";
            throw std::runtime_error(msg.str().c_str());
        }

        return x / (x.size() + 1.0);
    }
    
    //! @name Pairwise dependence measures
    //! @param x an \f$ n \times 2 \f$ matrix of observations.
    //! @{
      
    //! calculates the pairwise Kendall's \f$ \tau \f$.
    double pairwise_ktau(Eigen::Matrix<double, Eigen::Dynamic, 2>& x)
    {
        double tau;
        int n = (int) x.rows();
        int two = 2;
        ktau_matrix(x.data(), &two, &n, &tau);
        return tau;
    }

    //! calculates the pairwise correlation.
    double pairwise_cor(const Eigen::Matrix<double, Eigen::Dynamic, 2>& x)
    {
        double rho;
        auto z = x.rowwise() - x.colwise().mean();
        Eigen::MatrixXd sigma = z.adjoint() * z;
        rho = sigma(1,0) / sqrt(sigma(0,0) * sigma(1,1));

        return rho;
    }

    //! calculates the pair-wise Hoeffding's D.
    double pairwise_hoeffd(Eigen::Matrix<double, Eigen::Dynamic, 2> x)
    {
        size_t n = x.rows();

        // Compute the ranks
        auto R = to_pseudo_obs(x);
        R = (n+1.0)*R;

        // Compute Q, with Qi the number of points with both columns less than
        // their ith value
        Eigen::VectorXd Q(n);
        Eigen::Matrix<double, Eigen::Dynamic, 2> tmp = Eigen::MatrixXd::Ones(n, 2);
        for(size_t i=0; i<n; i++) {
            tmp.col(0) = Eigen::VectorXd::Constant(n,x(i,0));
            tmp.col(1) = Eigen::VectorXd::Constant(n,x(i,1));
            tmp = (x-tmp).unaryExpr([](double v){
                double res = 0.0;
                if(v < 0.0) {
                    res = 1.0;
                }
                return res;
            });
            Q(i) = tmp.rowwise().prod().sum();
        }

        Eigen::Matrix<double, Eigen::Dynamic, 2> ones = Eigen::MatrixXd::Ones(n, 2);
        double A = (R-ones).cwiseProduct(R-2*ones).rowwise().prod().sum();
        double B = (R-2*ones).rowwise().prod().cwiseProduct(Q).sum();
        double C = Q.cwiseProduct(Q-ones.col(0)).sum();

        double D = (A - 2*(n-2)*B + (n-2)*(n-3)*C);
        D /= (n*(n-1)*(n-2)*(n-3)*(n-4));

        return 30.0*D;
    }
    
    //! @}
}
