// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_stats.hpp>
#include <vinecopulib/misc/tools_stl.hpp>
#include <vinecopulib/misc/tools_c.h>

//! @file misc/tools_stats.cpp

namespace vinecopulib {

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
    double pairwise_tau(Eigen::Matrix<double, Eigen::Dynamic, 2>& x)
    {
        double tau;
        int n = (int) x.rows();
        int two = 2;
        ktau_matrix_c(x.data(), &two, &n, &tau);
        return tau;
    }

    //! calculates the pairwise Pearson correlation.
    double pairwise_cor(const Eigen::Matrix<double, Eigen::Dynamic, 2>& x)
    {
        double rho;
        auto z = x.rowwise() - x.colwise().mean();
        Eigen::MatrixXd sigma = z.adjoint() * z;
        rho = sigma(1,0) / sqrt(sigma(0,0) * sigma(1,1));

        return rho;
    }

    //! calculates the pairwise Spearman's \f$ \rho \f$.
    double pairwise_rho(Eigen::Matrix<double, Eigen::Dynamic, 2> x)
    {
        x = to_pseudo_obs(x);
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
     
    //! calculates a matrix of pairwise dependence measures.
    //! @param x an \f$ n \times d \f$ matrix of observations.
    //! @param measure either `"cor"` for Pearson correlation, `"tau"` for
    //!     Kendall's \f$ \tau \f$, `"rho"` for Spearman's  \f$ \rho \f$,
    //!     or `"hoeffd"` for Hoeffding's \f$ D \f$.
    //! @return a quadratic matrix of pairwise dependence measures.
    Eigen::MatrixXd dependence_matrix(const Eigen::MatrixXd& x, 
                                      const std::string& measure)
    {
        int n = x.rows();
        int d = x.cols();
        Eigen::MatrixXd mat(d, d);
        mat.diagonal() = Eigen::VectorXd::Constant(d, 1.0);
        Eigen::Matrix<double, Eigen::Dynamic, 2> pair_data(n, 2);
        for (int i = 1; i < d; ++i) {
            for (int j = 0; j < i; ++j) {
                pair_data.col(0) = x.col(i);
                pair_data.col(1) = x.col(j);
                if (measure == "tau") {
                    mat(i, j) = pairwise_tau(pair_data);
                } else if (measure == "cor") {
                    mat(i, j) = pairwise_cor(pair_data);
                } else if (measure == "rho") {
                    mat(i, j) = pairwise_rho(pair_data);
                } else if (measure == "hoeffd") {
                    mat(i, j) = pairwise_hoeffd(pair_data);
                } else {
                    throw std::runtime_error("measure not implemented");
                }
                mat(j, i) = mat(i, j);
            }
        }
        
        return mat;
    }

    //! Maximal dimension allowed for generalized Halton quasi Monte Carlo.
    #define ghalton_max_dim 360

    //! Primes for ghalton()
    static Eigen::Matrix<int, ghalton_max_dim, 1> primes = [] {
        Eigen::Matrix<int, ghalton_max_dim, 1> tmp;
        tmp << 2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71, 73,79,83,89,97,101,
                103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,
                199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,
                313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,
                433,439,443,449,457,461,463,467,479,487,491,499,503,509,521,523,541,547,557,
                563,569,571,577,587,593,599,601,607,613,617,619,631,641,643,647,653,659,661,
                673,677,683,691,701,709,719,727,733,739,743,751,757,761,769,773,787,797,809,
                811,821,823,827,829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,
                941,947,953,967,971,977,983,991,997,1009,1013,1019,1021,1031,1033,1039,1049,
                1051,1061,1063,1069,1087,1091,1093,1097,1103,1109,1117,1123,1129,1151,1153,
                1163,1171,1181,1187,1193,1201,1213,1217,1223,1229,1231,1237,1249,1259,1277,
                1279,1283,1289,1291,1297,1301,1303,1307,1319,1321,1327,1361,1367,1373,1381,
                1399,1409,1423,1427,1429,1433,1439,1447,1451,1453,1459,1471,1481,1483,1487,
                1489,1493,1499,1511,1523,1531,1543,1549,1553,1559,1567,1571,1579,1583,1597,
                1601,1607,1609,1613,1619,1621,1627,1637,1657,1663,1667,1669,1693,1697,1699,
                1709,1721,1723,1733,1741,1747,1753,1759,1777,1783,1787,1789,1801,1811,1823,
                1831,1847,1861,1867,1871,1873,1877,1879,1889,1901,1907,1913,1931,1933,1949,
                1951,1973,1979,1987,1993,1997,1999,2003,2011,2017,2027,2029,2039,2053,2063,
                2069,2081,2083,2087,2089,2099,2111,2113,2129,2131,2137,2141,2143,2153,2161,
                2179,2203,2207,2213,2221,2237,2239,2243,2251,2267,2269,2273,2281,2287,2293,
                2297,2309,2311,2333,2339,2341,2347,2351,2357,2371,2377,2381,2383,2389,2393,
                2399,2411,2417,2423;
        return tmp;
    }();

    // Scrambling factors for ghalton()
    static Eigen::Matrix<int, ghalton_max_dim, 1> permTN2 = [] {
        Eigen::Matrix<int, ghalton_max_dim, 1> tmp;
        tmp << 1,1,3,3,4,9,7,5,9,18,18,8,13,31,9,19,36,33,21,44,43,61,60,56,26,71,32,77,26,95,
                92,47,29,61,57,69,115,63,92,31,104,126,50,80,55,152,114,80,83,97,95,150,148,55,
                80,192,71,76,82,109,105,173,58,143,56,177,203,239,196,143,278,227,87,274,264,84,
                226,163,231,177,95,116,165,131,156,105,188,142,105,125,269,292,215,182,294,152,
                148,144,382,194,346,323,220,174,133,324,215,246,159,337,254,423,484,239,440,362,
                464,376,398,174,149,418,306,282,434,196,458,313,512,450,161,315,441,549,555,431,
                295,557,172,343,472,604,297,524,251,514,385,531,663,674,255,519,324,391,394,533,
                253,717,651,399,596,676,425,261,404,691,604,274,627,777,269,217,599,447,581,640,
                666,595,669,686,305,460,599,335,258,649,771,619,666,669,707,737,854,925,818,424,
                493,463,535,782,476,451,520,886,340,793,390,381,274,500,581,345,363,1024,514,
                773,932,556,954,793,294,863,393,827,527,1007,622,549,613,799,408,856,601,1072,
                938,322,1142,873,629,1071,1063,1205,596,973,984,875,918,1133,1223,933,1110,1228,
                1017,701,480,678,1172,689,1138,1022,682,613,635,984,526,1311,459,1348,477,716,
                1075,682,1245,401,774,1026,499,1314,743,693,1282,1003,1181,1079,765,815,1350,
                1144,1449,718,805,1203,1173,737,562,579,701,1104,1105,1379,827,1256,759,540,
                1284,1188,776,853,1140,445,1265,802,932,632,1504,856,1229,1619,774,1229,1300,
                1563,1551,1265,905,1333,493,913,1397,1250,612,1251,1765,1303,595,981,671,1403,
                820,1404,1661,973,1340,1015,1649,855,1834,1621,1704,893,1033,721,1737,1507,1851,
                1006,994,923,872,1860;
        return tmp;
    }();

    //! simulates from the multivariate Generalized Halton Sequence
    //!
    //! For more information on Generalized Halton Sequence, see
    //! Faure, H., Lemieux, C. (2009). Generalized Halton Sequences in 2008:
    //! A Comparative Study. ACM-TOMACS 19(4), Article 15.
    //!
    //! @param n number of observations.
    //! @param d dimension.
    //!
    //! @return An \f$ n \times d \f$ matrix of quasi-random
    //! \f$ \mathrm{U}[0, 1] \f$ variables.
    Eigen::MatrixXd ghalton(size_t n, size_t d) {

        Eigen::MatrixXd res(d, n);

        // Coefficients of the shift
        Eigen::MatrixXi shcoeff(d, 32);
        Eigen::VectorXi base = primes.block(0,0,d,1);
        Eigen::MatrixXd u =  Eigen::VectorXd::Zero(d, 1);
        auto U = simulate_uniform(d,32);
        for(int k = 31; k >= 0; k--) {
            shcoeff.col(k) = (base.cast <double>()).cwiseProduct(U.block(0,k,d,1)).cast <int> ();
            u = (u + shcoeff.col(k).cast <double> ()).cwiseQuotient(base.cast <double>());
        }
        res.block(0,0,d,1) = u;

        Eigen::VectorXi perm = permTN2.block(0,0,d,1);
        Eigen::MatrixXi coeff(d,32);
        Eigen::VectorXi tmp(d);
        int k;
        auto mod = [] (const int& u, const int& v) {return u % v;};
        for(size_t i = 1; i < n; i++) {

            // Find i in the prime base
            tmp = Eigen::VectorXi::Constant(d, (int) i);
            coeff = Eigen::MatrixXi::Zero(d, 32);
            k = 0;
            while((tmp.maxCoeff()  > 0) && (k < 32)) {
                coeff.col(k) = tmp.binaryExpr(base, mod);
                tmp = tmp.cwiseQuotient(base);
                k++;
            }

            u =  Eigen::VectorXd::Zero(d);
            k = 31;
            while(k >= 0) {
                tmp = perm.cwiseProduct(coeff.col(k)) + shcoeff.col(k);
                u = u + tmp.binaryExpr(base, mod).cast <double>();
                u = u.cwiseQuotient(base.cast <double>());
                k--;
            }
            res.block(0,i,d,1) = u;
        }

        return res.transpose();
    }
}

}
