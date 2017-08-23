// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/bicop/tools_select.hpp>
#include <vinecopulib/misc/tools_stats.hpp>
#include <vinecopulib/misc/tools_stl.hpp>
#include <cmath>

std::vector<double> get_c1c2(const Eigen::Matrix<double, Eigen::Dynamic, 2>& data,
                             double tau)
{
    using namespace vinecopulib::tools_stats;
    size_t n = data.rows();
    Eigen::MatrixXd x = Eigen::MatrixXd::Zero(n, 2);
    Eigen::MatrixXd z1 = x;
    Eigen::MatrixXd z2 = x;
    x = qnorm(data);
    
    int count1 = 0, count2 = 0;
    for (size_t j = 0; j < n; ++j) {
        if (tau > 0) {
            if ((x(j, 0) > 0) && (x(j, 1) > 0)) {
                z1.row(count1) = x.row(j);
                ++count1;
            }
            if ((x(j, 0) < 0) && (x(j, 1) < 0)) {
                z2.row(count2) = x.row(j);
                ++count2;
            }
        } else {
            if ((x(j, 0) < 0) && (x(j, 1) > 0)) {
                z1.row(count1) = x.row(j);
                ++count1;
            }
            if ((x(j, 0) > 0) && (x(j, 1) < 0)) {
                z2.row(count2) = x.row(j);
                ++count2;
            }
        }
    }
    
    // if one of the quadrants is empty, we see it as independent
    double c1, c2;
    if (count1 == 0) {
        c1 = 0.0;
    } else {
        c1 = pairwise_cor(z1.block(0, 0, count1 - 1, 2));
    }
    if (count2 == 0) {
        c2 = 0.0;
    } else {
        c2 = pairwise_cor(z2.block(0, 0, count2 - 1, 2));
    }

    return {c1, c2};
}

bool preselect_family(std::vector<double> c, double tau, 
                      vinecopulib::BicopFamily family,  int rotation, 
                      bool is_rotationless)
{
    using namespace vinecopulib;
    using namespace tools_stl;

    bool preselect = false;
    if (is_rotationless) {
        preselect = true;
        if ((std::fabs(c[0] - c[1]) > 0.3) & (family == BicopFamily::frank))
            preselect = false;
    } else {
        if (is_member(family, bicop_families::BB)) {
            if ((tau > 0) && is_member(rotation, {0, 180})) {
                preselect = true;
            }
            if ((tau < 0) && is_member(rotation, {90, 270})) {
                preselect = true;
            }
        }
        bool is_90or180 = is_member(rotation, {90, 180});
        if (c[0] - c[1] > 0.05) {
            if (is_member(family, bicop_families::lt) & is_90or180) { 
                preselect = true;
            }
            if (is_member(family, bicop_families::ut) & !is_90or180) {
                preselect = true;
            }
        } else if (c[0] - c[1] < -0.05) {
            if (is_member(family, bicop_families::lt) & !is_90or180) { 
                preselect = true;
            }
            if (is_member(family, bicop_families::ut) & is_90or180) {
                preselect = true;
            }
        } else {
            if ((tau > 0) && is_member(rotation, {0, 180})) {
                preselect = true;
            }
            if ((tau < 0) && is_member(rotation, {90, 270})) {
                preselect = true;
            }
        }
    }
    
    return preselect;
}
