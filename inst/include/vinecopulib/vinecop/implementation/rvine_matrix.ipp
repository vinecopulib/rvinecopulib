// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_stl.hpp>
#include <vinecopulib/misc/tools_parallel.hpp>

namespace vinecopulib {
//! instantiates an RVineMatrix object.
//! @param matrix a valid R-vine matrix.
//! @param check whether the matrix shall be checked for validity.
inline RVineMatrix::RVineMatrix(
    const Eigen::Matrix <size_t, Eigen::Dynamic, Eigen::Dynamic> &matrix,
    bool check
)
{
    d_ = matrix.rows();
    matrix_ = matrix;
    if (check) {
        check_if_quadratic();
        check_lower_tri();
        check_upper_tri();
        check_antidiagonal();
        check_columns();
        check_proximity_condition();
    }
}

//! extract matrix(row, col)
//!
//! \param row
//! \param col
//! \return matrix(row, col)
inline size_t RVineMatrix::get_element(size_t row, size_t col) const
{
    if (row >= d_ || col >= d_) {
        throw std::runtime_error("row and col should be < d");
    }
    return matrix_(row, col);
}

//! extract the matrix.
//!
inline Eigen::Matrix <size_t, Eigen::Dynamic, Eigen::Dynamic>
RVineMatrix::get_matrix() const
{
    return matrix_;
}

//! extracts the variable order in the R-vine.
//!
inline Eigen::Matrix<size_t, Eigen::Dynamic, 1> RVineMatrix::get_order() const
{
    return matrix_.colwise().reverse().diagonal().reverse();
}

//! constructs a D-vine matrix, where each tree is a path.
//!
//! @param order order of the variables.
inline Eigen::Matrix <size_t, Eigen::Dynamic, Eigen::Dynamic>
RVineMatrix::construct_d_vine_matrix(
    const Eigen::Matrix<size_t, Eigen::Dynamic, 1> &order)
{
    size_t d = order.size();
    Eigen::Matrix <size_t, Eigen::Dynamic, Eigen::Dynamic> vine_matrix(d, d);
    vine_matrix.fill(0);

    for (size_t i = 0; i < d; ++i) {
        vine_matrix(d - 1 - i, i) = order(d - 1 - i);  // diagonal
    }

    for (size_t i = 1; i < d; ++i) {
        for (size_t j = 0; j < i; ++j) {
            vine_matrix(d - 1 - i, j) = order(i - j - 1);  // below diagonal
        }
    }

    return vine_matrix;
}

//! check whether an edge belong to the implied structure
//!
//! @param conditioned the conditioned set.
//! @param conditioning the conditioning set.
inline bool
RVineMatrix::belongs_to_structure(const std::vector <size_t> conditioned,
                                  const std::vector <size_t> conditioning)
{
    if (conditioned.size() != 2) {
        throw std::runtime_error("conditioned should have size 2 ");
    }

    size_t tree = conditioning.size();
    std::vector <size_t> conditioning_mat(tree);
    std::vector <size_t> conditioned_mat(2);
    if (tree + 2 <= d_) {
        for (size_t i = 0; i < d_ - tree - 1; ++i) {
            conditioned_mat = {matrix_(d_ - 1 - i, i), matrix_(tree, i)};
            if (conditioned == conditioned_mat) {
                // conditioned sets equal, need to check conditioning set
                if (tree == 0) {
                    return true;  // conditioning set is empty for tree == 0
                }
                auto cond = matrix_.col(i).head(tree);
                Eigen::Matrix<size_t, Eigen::Dynamic, 1>::Map(
                    &conditioning_mat[0], tree) = cond;
                if (tools_stl::is_same_set(conditioning, conditioning_mat)) {
                    return true;
                }
            }
        }
    }

    // edge is not contained in the structure implied by the matrix
    return false;
}


//! extracts the R-vine matrix in natural order, which means that the
//! counter-diagonal has entries (d, ..., 1). We
//! convert to natural order by relabeling the variables. Most algorithms for
//! estimation and evaluation assume that the R-vine matrix is in natural order.
inline Eigen::Matrix <size_t, Eigen::Dynamic, Eigen::Dynamic>
RVineMatrix::in_natural_order() const
{
    // create vector of new variable labels: d, ..., 1
    std::vector <size_t> ivec = tools_stl::seq_int(1, d_);
    tools_stl::reverse(ivec);
    Eigen::Map <Eigen::Matrix<size_t, Eigen::Dynamic, 1>> new_labels(&ivec[0],
                                                                     d_);

    return relabel_elements(matrix_, new_labels);
}

//! extracts the maximum matrix, which is derived from an R-vine matrix by
//! iteratively computing
//! the (elementwise) maximum of a row and the row below (starting from the
//! bottom). It is used in estimation and evaluation algorithms to find the right
//! pseudo observations for an edge.
inline Eigen::Matrix <size_t, Eigen::Dynamic, Eigen::Dynamic>
RVineMatrix::get_max_matrix() const
{
    auto max_matrix = this->in_natural_order();
    for (size_t i = 0; i < d_ - 1; ++i) {
        for (size_t j = 0; j < d_ - i - 1; ++j) {
            max_matrix(i + 1, j) = max_matrix.block(i, j, 2, 1).maxCoeff();
        }
    }
    return max_matrix;
}

//! extracts a matrix indicating which of the first h-functions are needed
//! (it is usually not necessary to apply both h-functions for each
//! pair-copula).
inline MatrixXb RVineMatrix::get_needed_hfunc1() const
{
    MatrixXb needed_hfunc1 = MatrixXb::Constant(d_, d_, false);

    auto no_matrix = this->in_natural_order();
    auto max_matrix = this->get_max_matrix();
    for (size_t i = 1; i < d_ - 1; ++i) {
        size_t j = d_ - i;
        MatrixXb isnt_mat_j = (no_matrix.block(0, 0, j, i).array() != j);
        MatrixXb is_max_j = (max_matrix.block(0, 0, j, i).array() == j);
        MatrixXb is_different = (isnt_mat_j.array() && is_max_j.array());
        needed_hfunc1.block(0, i, j, 1) = is_different.rowwise().any();
    }
    return needed_hfunc1;
}

//! extracts a matrix indicating which of the second h-functions are needed
//! (it is usually not necessary to apply both h-functions for each
//! pair-copula).
inline MatrixXb RVineMatrix::get_needed_hfunc2() const
{
    MatrixXb needed_hfunc2 = MatrixXb::Constant(d_, d_, false);
    needed_hfunc2.block(0, 0, d_ - 1, 1) = MatrixXb::Constant(d_ - 1, 1, true);
    auto no_matrix = this->in_natural_order();
    auto max_matrix = this->get_max_matrix();
    for (size_t i = 1; i < d_ - 1; ++i) {
        size_t j = d_ - i;
        // fill column i with true above the diagonal
        needed_hfunc2.block(0, i, d_ - i, 1) = MatrixXb::Constant(d_ - i, 1,
                                                                  true);
        // for diagonal, check whether matrix and maximum matrix coincide
        MatrixXb is_mat_j = (no_matrix.block(j - 1, 0, 1, i).array() == j);
        MatrixXb is_max_j = (max_matrix.block(j - 1, 0, 1, i).array() == j);
        needed_hfunc2(j - 1, i) = (is_mat_j.array() && is_max_j.array()).any();
    }

    return needed_hfunc2;
}
//! @}


inline void RVineMatrix::check_if_quadratic() const
{
    std::string problem = "must be quadratic.";
    if (matrix_.rows() != matrix_.cols()) {
        throw std::runtime_error("not a valid R-vine matrix: " + problem);
    }
}

inline void RVineMatrix::check_lower_tri() const
{
    std::string problem = "the lower right triangle must only contain zeros";
    size_t sum_lwr = 0;
    for (size_t j = 1; j < d_; ++j) {
        sum_lwr += matrix_.block(d_ - j, j, j, 1).array().sum();
        if (sum_lwr != 0) {
            throw std::runtime_error("not a valid R-vine matrix: " + problem);
        }
    }
}

inline void RVineMatrix::check_upper_tri() const
{
    std::string problem;
    problem += "the upper left triangle can only contain numbers ";
    problem += "between 1 and d (number of variables).";
    size_t min_upr = d_;
    size_t max_upr = 0;
    for (size_t j = 0; j < d_; ++j) {
        min_upr = std::min(min_upr, matrix_.col(j).head(d_ - j).minCoeff());
        max_upr = std::max(max_upr, matrix_.col(j).head(d_ - j).maxCoeff());
        if ((max_upr > d_) | (min_upr < 1)) {
            throw std::runtime_error("not a valid R-vine matrix: " + problem);
        }
    }
}

inline void RVineMatrix::check_antidiagonal() const
{
    std::string problem;
    problem += "the antidiagonal must contain the numbers ";
    problem += "1, ..., d (the number of variables)";
    auto diag = matrix_.colwise().reverse().diagonal();
    std::vector <size_t> diag_vec(d_);
    Eigen::Matrix<size_t, Eigen::Dynamic, 1>::Map(&diag_vec[0], d_) = diag;
    if (!tools_stl::is_same_set(diag_vec, tools_stl::seq_int(1, d_))) {
        throw std::runtime_error("not a valid R-vine matrix: " + problem);
    }
}

inline void RVineMatrix::check_columns() const
{
    using namespace tools_stl;
    auto no_matrix = in_natural_order();
    std::string problem;
    problem += "the antidiagonal entry of a column must not be ";
    problem += "contained in any column further to the right; ";
    problem += "the entries of a column must be contained ";
    problem += "in all columns to the left.";

    // In natural order: column j only contains indices 1:(d - j).
    bool ok = true;
    for (size_t j = 0; j < d_; ++j) {
        std::vector <size_t> col_vec(d_ - j);
        Eigen::Matrix<size_t, Eigen::Dynamic, 1>::Map(&col_vec[0], d_ - j) =
            no_matrix.col(j).head(d_ - j);
        ok = ok & is_same_set(col_vec, seq_int(1, d_ - j));
        if (!ok) {
            throw std::runtime_error("not a valid R-vine matrix: " + problem);
        }
    }
}

inline void RVineMatrix::check_proximity_condition() const
{
    using namespace tools_stl;
    for (size_t t = 1; t < d_ - 1; ++t) {
        for (size_t e = 0; e < d_ - t - 1; ++e) {
            // non-diagonal conditioned variable
            double v0 = matrix_(t, e);
            // conditioning set
            std::vector <size_t> D0(t);
            Eigen::Matrix<size_t, Eigen::Dynamic, 1>::Map(&D0[0], t) =
                matrix_.col(e).head(t);

            // the pair (v0, D0) has to be matched in another column
            bool found = false;

            // search antidiagonal right of column e for v0,
            // columns with less then t entries can be omitted
            for (size_t j = e + 1; j < d_ - t; ++j) {
                if (matrix_(d_ - j - 1, j) != v0) {
                    continue;
                }
                // check if conditioning sets coincide
                std::vector <size_t> D(t);
                Eigen::Matrix<size_t, Eigen::Dynamic, 1>::Map(&D[0], t) =
                    matrix_.col(j).head(t);
                if (is_same_set(D0, D)) {
                    found = true;
                    break;
                }
            }

            // if (v0, D0) was matched, continue with next pair-copula
            if (found) {
                continue;
            }

            // otherwise search row above t,
            // columns with less then t entries can be omitted
            for (size_t j = e + 1; j < d_ - t; ++j) {
                if (matrix_(t - 1, j) != v0) {
                    continue;
                }
                // check if conditioning sets coincide
                std::vector <size_t> D(t);
                Eigen::Matrix<size_t, Eigen::Dynamic, 1>::Map(&D[0], t - 1) =
                    matrix_.col(j).head(t - 1);
                D[t - 1] = matrix_(d_ - j - 1, j); // add antidiagonal element
                if (is_same_set(D0, D)) {
                    found = true;
                    break;
                }
            }

            if (!found) {
                std::stringstream problem;
                problem <<
                        "not a valid R-vine matrix: " <<
                        "proximity contition violated; " <<
                        "cannot extract conditional distribution (" <<
                        v0 << " | ";
                for (size_t i = 0; i < D0.size() - 1; ++i) {
                    problem << D0[i] << ", ";
                }
                problem << D0[D0.size() - 1] << ") from pair-copulas.";
                throw std::runtime_error(problem.str().c_str());
            }
        }
    }
}


// translates matrix_entry from old to new labels
inline size_t relabel_one(size_t x,
                          const Eigen::Matrix<size_t, Eigen::Dynamic, 1> &old_labels,
                          const Eigen::Matrix<size_t, Eigen::Dynamic, 1> &new_labels)
{
    for (int i = 0; i < old_labels.size(); ++i) {
        if (x == old_labels[i]) {
            return new_labels[i];
        }
    }
    return 0;
}

// relabels all elements of the matrix (upper triangle assumed to be 0)
inline Eigen::Matrix <size_t, Eigen::Dynamic, Eigen::Dynamic> relabel_elements(
    const Eigen::Matrix <size_t, Eigen::Dynamic, Eigen::Dynamic> &matrix,
    const Eigen::Matrix<size_t, Eigen::Dynamic, 1> &new_labels)
{
    size_t d = matrix.rows();
    auto old_labels = matrix.colwise().reverse().diagonal();
    auto new_matrix = matrix;
    for (size_t i = 0; i < d; ++i) {
        for (size_t j = 0; j < d - i; ++j) {
            new_matrix(i, j) = relabel_one(matrix(i, j), old_labels,
                                           new_labels);
        }
    }

    return new_matrix;
}

//! completes an R-vine matrix that has only been filled for the first couple
//! of trees.
//! @param mat R-vine matrix where counter-diagonal and upper block have been 
//!     filled; entries in between may still be 0.
//! @param t_start the tree level from where the matrix needs to be filled.
//! @return a valid R-vine structure matrix (= passing all 
//!     checks in `RVineMatrix()`).
inline void RVineMatrix::complete_matrix(
    Eigen::Matrix <size_t, Eigen::Dynamic, Eigen::Dynamic> &mat,
    size_t t_start,
    size_t num_threads)
{
    using namespace tools_stl;
    size_t d = mat.cols();
    // there are max d-2 off-diagonal entries in each column
    for (size_t t = t_start; t < d - 1; ++t) {
        // for each column, get all indices it already contains
        std::vector <std::vector<size_t>> all_indices(d - t);
        std::vector <size_t> tmp(t);
        for (size_t e = 0; e < d - t; e++) {
            Eigen::Matrix<size_t, Eigen::Dynamic, 1>::Map(&tmp[0], t) =
                mat.col(e).head(t);
            all_indices[e] = cat(mat(d - 1 - e, e), tmp);
        }

        // fill row t column by column
        auto complete_column = [&] (size_t e0) {
            // check all columns to the right of e0
            for (size_t e1 = e0 + 1; e1 < d - t; e1++) {
                // get possible conditioned sets joining edges e0 and e1
                auto ned_set = set_sym_diff(all_indices[e0], all_indices[e1]);
                // a valid conditioned set has two elements
                if (ned_set.size() == 2) {
                    // allowed entry in column e0 is the element of the
                    // conditioned set that is not contained in column e0
                    auto new_ind = set_diff(ned_set, all_indices[e0]);
                    if (new_ind.size() == 1) {
                        mat(t, e0) = new_ind[0];
                        break;  // continue with next column
                    }
                }
            }
        };
        
        tools_parallel::map_on_pool(complete_column, 
                                    seq_int(0, d - t), 
                                    num_threads);
    }
}
}
