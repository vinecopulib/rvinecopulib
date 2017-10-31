// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/misc/tools_eigen.hpp>

namespace vinecopulib {
using namespace tools_eigen;

//! @brief A class for regular vine matrices.
//!
//! A regular vine (R-vine) matrix encodes the structure of a vine.
//!     An examplary matrix is
//! ```
//! 1 1 1 1
//! 2 2 2 0
//! 3 3 0 0
//! 4 0 0 0
//! ```
//! which encodes the following pair-copulas:
//! ```
//! | tree | edge | pair-copulas   |
//! |------|------|----------------|
//! | 0    | 0    | `(4, 1)`       |
//! |      | 1    | `(3, 1)`       |
//! |      | 2    | `(2, 1)`       |
//! | 1    | 0    | `(4, 2; 1)`    |
//! |      | 1    | `(3, 2; 1)`    |
//! | 2    | 0    | `(4, 3; 2, 1)` |
//! ```
//! Denoting by `M[i, j]` the matrix entry in row `i` and column `j`,
//! the pair-copula index for edge `e` in tree `t` of a `d` dimensional vine
//! is `(M[d - 1 - t, e], M[t, e]; M[t - 1, e], ..., M[0, e])`. Less
//! formally,
//! 1. Start with the counter-diagonal element of column `e` (first conditioned
//!    variable).
//! 2. Jump up to the element in row `t` (second conditioned variable).
//! 3. Gather all entries further up in column `e` (conditioning set).
//!
//! A valid R-vine matrix must satisfy several conditions which are checked
//! when `RVineMatrix()` is called:
//! 1. The lower right triangle must only contain zeros.
//! 2. The upper left triangle can only contain numbers between 1 and d.
//! 3. The antidiagonal must contain the numbers 1, ..., d.
//! 4. The antidiagonal entry of a column must not be contained in any
//!    column further to the right.
//! 5. The entries of a column must be contained in all columns to the left.
//! 6. The proximity condition must hold: For all t = 1, ..., d - 2 and
//!    e = 0, ..., d - t - 1 there must exist an index j > d, such that
//!    `(M[t, e], {M[0, e], ..., M[t-1, e]})` equals either
//!    `(M[d-j-1, j], {M[0, j], ..., M[t-1, j]})` or
//!    `(M[t-1, j], {M[d-j-1, j], M[0, j], ..., M[t-2, j]})`.
//!
//! Condition 6 already implies conditions 2-5, but is more difficult to
//! check by hand.
class RVineMatrix
{
public:
    RVineMatrix()
    {
    }

    RVineMatrix(
        const Eigen::Matrix <size_t, Eigen::Dynamic, Eigen::Dynamic> &matrix,
        bool check = true
    );

    size_t get_element(size_t row, size_t col) const;

    Eigen::Matrix <size_t, Eigen::Dynamic, Eigen::Dynamic> get_matrix() const;

    Eigen::Matrix<size_t, Eigen::Dynamic, 1> get_order() const;

    Eigen::Matrix <size_t, Eigen::Dynamic, Eigen::Dynamic>
    in_natural_order() const;

    Eigen::Matrix <size_t, Eigen::Dynamic, Eigen::Dynamic>
    get_max_matrix() const;

    MatrixXb get_needed_hfunc1() const;

    MatrixXb get_needed_hfunc2() const;

    static Eigen::Matrix <size_t, Eigen::Dynamic, Eigen::Dynamic>
    construct_d_vine_matrix(
        const Eigen::Matrix<size_t, Eigen::Dynamic, 1> &order);

    bool belongs_to_structure(const std::vector <size_t> conditioned,
                              const std::vector <size_t> conditioning);

    static void complete_matrix(
        Eigen::Matrix <size_t, Eigen::Dynamic, Eigen::Dynamic> &mat,
        size_t t_start);

private:
    void check_if_quadratic() const;

    void check_lower_tri() const;

    void check_upper_tri() const;

    void check_antidiagonal() const;

    void check_columns() const;

    void check_proximity_condition() const;

    size_t d_;
    Eigen::Matrix <size_t, Eigen::Dynamic, Eigen::Dynamic> matrix_;
};

int relabel_one(const int &x,
                const Eigen::Matrix<size_t, Eigen::Dynamic, 1> &old_labels,
                const Eigen::Matrix<size_t, Eigen::Dynamic, 1> &new_labels);

Eigen::Matrix <size_t, Eigen::Dynamic, Eigen::Dynamic> relabel_elements(
    const Eigen::Matrix <size_t, Eigen::Dynamic, Eigen::Dynamic> &matrix,
    const Eigen::Matrix<size_t, Eigen::Dynamic, 1> &new_labels);
}

#include <vinecopulib/vinecop/implementation/rvine_matrix.ipp>
