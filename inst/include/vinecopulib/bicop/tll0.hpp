// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/bicop/kernel.hpp>

namespace vinecopulib
{
    //! @brief The transformation local-constant likelihood estimator
    //!
    //! This class is used in the implementation underlying the Bicop class. 
    //! Users should not use AbstractBicop or derived classes directly, but 
    //! always work with the Bicop interface.
    //! 
    //! @literature
    //! Nagler, Thomas. *kdecopula: An R Package for the Kernel Estimation of 
    //! Copula Densities*. arXiv:1603.04229 [stat.CO], 2016
    class Tll0Bicop : public KernelBicop
    {
    public:
        Tll0Bicop();
    private:
        Eigen::VectorXd gaussian_kernel_2d(
            const Eigen::Matrix<double, Eigen::Dynamic, 2>& x);
        void fit(const Eigen::Matrix<double, Eigen::Dynamic, 2>& data, 
                 std::string, double mult);
    };
}
