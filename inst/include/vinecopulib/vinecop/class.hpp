// Copyright © 2018 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/bicop/class.hpp>
#include <vinecopulib/vinecop/rvine_matrix.hpp>
#include <vinecopulib/vinecop/tools_select.hpp>

namespace vinecopulib {
//! @brief A class for vine copula models
//!
//! A vine copula model is characterized by the structure matrix (see
//! RVineMatrix) and the pair-copulas.
class Vinecop
{
public:
    // Constructors
    Vinecop()
    {
    }

    Vinecop(size_t d);

    Vinecop(const Eigen::Matrix <size_t, Eigen::Dynamic, Eigen::Dynamic> &matrix,
            bool check_matrix = true);

    Vinecop(const std::vector <std::vector<Bicop>> &pair_copulas,
            const Eigen::Matrix <size_t, Eigen::Dynamic, Eigen::Dynamic> &matrix,
            bool check_matrix = true);

    Vinecop(const Eigen::MatrixXd &data,
            FitControlsVinecop controls = FitControlsVinecop());

    Vinecop(const Eigen::MatrixXd &data,
            const Eigen::Matrix <size_t, Eigen::Dynamic, Eigen::Dynamic> &matrix,
            FitControlsVinecop controls = FitControlsVinecop(),
            bool check_matrix = true);

    Vinecop(const char *filename, bool check_matrix = true);

    Vinecop(boost::property_tree::ptree input, bool check_matrix = true);

    // Serialize
    boost::property_tree::ptree to_ptree();

    void to_json(const char *filename);

    // Methods modifying structure and/or families and parameters
    void select_all(const Eigen::MatrixXd &data,
                    FitControlsVinecop controls = FitControlsVinecop());

    void select_families(const Eigen::MatrixXd &data,
                         FitControlsVinecop controls = FitControlsVinecop());

    // Getters for a single pair copula
    Bicop get_pair_copula(size_t tree, size_t edge) const;

    BicopFamily get_family(size_t tree, size_t edge) const;

    int get_rotation(size_t tree, size_t edge) const;

    Eigen::VectorXd get_parameters(size_t tree, size_t edge) const;

    Eigen::Matrix <size_t, Eigen::Dynamic, Eigen::Dynamic> get_matrix() const;

    // Getters for all pair copulas
    std::vector <std::vector<Bicop>> get_all_pair_copulas() const;

    std::vector <std::vector<BicopFamily>> get_all_families() const;

    std::vector <std::vector<int>> get_all_rotations() const;

    std::vector <std::vector<Eigen::VectorXd>> get_all_parameters() const;
    
    // getter for the threshold
    double get_threshold() const;

    // Stats methods
    Eigen::VectorXd pdf(const Eigen::MatrixXd &u) const;

    Eigen::VectorXd cdf(const Eigen::MatrixXd &u, const size_t N = 1e4) const;

    Eigen::MatrixXd simulate(size_t n) const;

    Eigen::MatrixXd inverse_rosenblatt(const Eigen::MatrixXd &u) const;

    // Fit statistics
    double calculate_npars() const;

    double loglik(const Eigen::MatrixXd &u) const;

    double aic(const Eigen::MatrixXd &u) const;

    double bic(const Eigen::MatrixXd &u) const;

    // Misc methods
    static std::vector <std::vector<Bicop>> make_pair_copula_store(size_t d,
                                                                   size_t truncation_level = std::numeric_limits<size_t>::max());

private:
    size_t d_;
    RVineMatrix vine_matrix_;
    std::vector <std::vector<Bicop>> pair_copulas_;
    double threshold_;

    void check_data_dim(const Eigen::MatrixXd &data);

    Eigen::Matrix<size_t, Eigen::Dynamic, 1> inverse_permutation(
        const Eigen::Matrix<size_t, Eigen::Dynamic, 1> &order) const;
};

}

#include <vinecopulib/vinecop/implementation/class.ipp>
