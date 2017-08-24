// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/vinecop/class.hpp>
#include <vinecopulib/misc/tools_stl.hpp>
#include <vinecopulib/misc/tools_stats.hpp>

#include <stdexcept>
#include <vector>

namespace vinecopulib
{
    //! creates a D-vine on `d` variables with all pair-copulas set to 
    //! independence.
    //! @param d the dimension (= number of variables) of the model.
    Vinecop::Vinecop(size_t d)
    {
        d_ = d;
        // D-vine with variable order (1, ..., d)
        Eigen::Matrix<size_t, Eigen::Dynamic, 1> order(d);
        for (size_t i = 0; i < d; ++i) {
            order(i) = i + 1;
        }
        vine_matrix_ = RVineMatrix(RVineMatrix::construct_d_vine_matrix(order),
                                   false);  // don't check if R-vine matrix

        // all pair-copulas are independence
        pair_copulas_ = make_pair_copula_store(d);
        for (auto& tree : pair_copulas_) {
            for (auto& pc : tree) {
                pc = Bicop(BicopFamily::indep);
            }
        }
    }
    
    //! creates a vine copula with structure specified by an R-vine matrix; all
    //! pair-copulas are set to independence.
    //! @param matrix an R-vine matrix.
    //! @param check_matrix whether to check if `matrix` is a valid R-vine 
    //!     matrix.
    Vinecop::Vinecop(const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& matrix,
        bool check_matrix)
    {
        d_ = matrix.rows();
        vine_matrix_ = RVineMatrix(matrix, check_matrix);

        // all pair-copulas are independence
        pair_copulas_ = make_pair_copula_store(d_);
        for (auto& tree : pair_copulas_) {
            for (auto& pc : tree) {
                pc = Bicop(BicopFamily::indep);
            }
        }
    }
    
    //! creates an arbitrary vine copula model.
    //! @param pair_copulas Bicop objects specifying the pair-copulas, see 
    //!     make_pair_copula_store().
    //! @param matrix an R-vine matrix specifying the vine structure.
    //! @param check_matrix whether to check if `matrix` is a valid R-vine 
    //!     matrix.
    Vinecop::Vinecop(const std::vector<std::vector<Bicop>>& pair_copulas,
        const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& matrix,
        bool check_matrix)
    {
        d_ = matrix.rows();
        if (pair_copulas.size() != d_ - 1) {
            std::stringstream message;
            message <<
                    "size of pair_copulas does not match dimension of matrix (" <<
                    d_ << "); " <<
                    "expected size: " << d_ - 1 << ", "<<
                    "actual size: " << pair_copulas.size() << std::endl;
            throw std::runtime_error(message.str().c_str());
        }
        for (size_t t = 0; t < d_ - 1; ++t) {
            if (pair_copulas[t].size() != d_ - 1 - t) {
                std::stringstream message;
                message <<
                        "size of pair_copulas[" << t << "] " <<
                        "does not match dimension of matrix (" << d_ << "); " <<
                        "expected size: " << d_ - 1 - t << ", "<<
                        "actual size: " << pair_copulas[t].size() << std::endl;
                throw std::runtime_error(message.str().c_str());
            }
        }

        vine_matrix_ = RVineMatrix(matrix, check_matrix);
        pair_copulas_ = pair_copulas;
    }

    //! creates from a boost::property_tree::ptree object
    //! @param input the boost::property_tree::ptree object to convert from
    //! (see to_ptree() for the structure of the input).
    //! @param check_matrix whether to check if the `"matrix"` node represents
    //!      a valid R-vine matrix.
    Vinecop::Vinecop(boost::property_tree::ptree input, bool check_matrix) {

        auto matrix = tools_serialization::ptree_to_matrix<size_t>(
                input.get_child("matrix"));
        vine_matrix_ = RVineMatrix(matrix, check_matrix);
        d_ = (size_t) matrix.rows();

        pair_copulas_ = make_pair_copula_store(d_);
        boost::property_tree::ptree pcs_node = input.get_child("pair copulas");
        for (size_t tree = 0; tree < d_ - 1; ++tree) {
            boost::property_tree::ptree tree_node = pcs_node.get_child(
                    "tree" + std::to_string(tree));
            for (size_t edge = 0; edge < d_ - tree - 1; ++edge) {
                boost::property_tree::ptree pc_node = tree_node.get_child(
                        "pc" + std::to_string(edge));
                pair_copulas_[tree][edge] = Bicop(pc_node);
            }
        }
    }

    //! creates from a JSON file
    //! @param filename the name of the JSON file to read (see to_ptree() for the
    //! structure of the file).
    //! @param check_matrix whether to check if the `"matrix"` node represents
    //!      a valid R-vine matrix.
    Vinecop::Vinecop(const char *filename, bool check_matrix) :
            Vinecop(tools_serialization::json_to_ptree(filename),
                    check_matrix) {}

    //! constructs a vine copula model from data. 
    //! 
    //! The function creates a model and calls select_family().
    //! 
    //! @param data an \f$ n \times d \f$ matrix of observations.
    //! @param matrix either an empty matrix (default) or an R-vine structure 
    //!     matrix, see select_families().
    //! @param controls see FitControlsVinecop.
    //! @param check_matrix whether to check if `matrix` is a valid R-vine 
    //!     matrix.
    Vinecop::Vinecop(const Eigen::MatrixXd& data,
        const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& matrix,
        FitControlsVinecop controls,
        bool check_matrix)
    {
        d_ = data.cols();
        pair_copulas_ = make_pair_copula_store(d_);
        vine_matrix_ = RVineMatrix(matrix, check_matrix);
        select_families(data, controls);
    }

    //! constructs a vine copula model from data.
    //!
    //! The function creates a model and calls select_all().
    //!
    //! @param data an \f$ n \times d \f$ matrix of observations.
    //! @param controls see FitControlsVinecop.
    Vinecop::Vinecop(const Eigen::MatrixXd& data, FitControlsVinecop controls)
    {
        d_ = data.cols();
        pair_copulas_ = make_pair_copula_store(d_);
        select_all(data, controls);
    }

    //! Convert the copula into a boost::property_tree::ptree object
    //!
    //! The ptree object contains two nodes : `"matrix"` and `"pair copulas"`.
    //! The former is encodes the R-Vine structure and the latter is a list of
    //! child nodes for the trees (`"tree1"`, `"tree2"`, etc), each containing
    //! a list of child nodes for the edges (`"pc1"`, `"pc2"`, etc).
    //! See Bicop::to_ptree() for the encoding of pair-copulas.
    //!
    //! @return the boost::property_tree::ptree object containing the copula.
    boost::property_tree::ptree Vinecop::to_ptree()
    {
        boost::property_tree::ptree output;

        boost::property_tree::ptree pair_copulas;
        for (size_t tree = 0; tree < d_ - 1; ++tree) {
            boost::property_tree::ptree tree_node;
            for (size_t edge = 0; edge < d_ - tree - 1; ++edge) {
                tree_node.add_child("pc"+std::to_string(edge),
                                    pair_copulas_[tree][edge].to_ptree());
            }
            pair_copulas.add_child("tree"+std::to_string(tree), tree_node);
        }

        output.add_child("pair copulas", pair_copulas);
        auto matrix_node = tools_serialization::matrix_to_ptree(get_matrix());
        output.add_child("matrix", matrix_node);

        return output;
    }

    //! Write the copula object into a JSON file
    //!
    //! See to_ptree() for the structure of the file.
    //!
    //! @param filename the name of the file to write.
    void Vinecop::to_json(const char *filename)
    {
        boost::property_tree::write_json(filename, to_ptree());
    }

    //! Initialize object for storing pair copulas
    //!
    //! @param d dimension of the vine copula.
    //! @return A nested vector such that `pc_store[t][e]` contains a Bicop 
    //!     object for the pair copula corresponding to tree `t` and edge `e`.
    std::vector<std::vector<Bicop>> Vinecop::make_pair_copula_store(size_t d)
    {
        if (d < 2) {
            throw std::runtime_error("the dimension should be larger than 1");
        }

        std::vector<std::vector<Bicop>> pc_store(d - 1);
        for (size_t t = 0; t < d - 1; ++t) {
            pc_store[t].resize(d - 1 - t);
        }
    
        return pc_store;
    }
    
    
    //! automatically fits and selects a vine copula model
    //! 
    //! Selection of the structure is performed using the algorithm of
    //! Dissmann, J. F., E. C. Brechmann, C. Czado, and D. Kurowicka (2013).
    //! *Selecting and estimating regular vine copulae and application to 
    //! financial returns.* Computational Statistics & Data Analysis, 59 (1), 
    //! 52-69.
    //!
    //! @param data nxd matrix of copula data.
    //! @param controls the controls to the algorithm (see FitControlsVinecop).
    void Vinecop::select_all(const Eigen::MatrixXd& data,
                             FitControlsVinecop controls)
    {
        check_data_dim(data);
        tools_select::StructureSelector selector(data, controls);
        if (controls.needs_sparse_select()) {
            selector.sparse_select_all_trees(data);
        } else {
            selector.select_all_trees(data);
        }
        vine_matrix_ = selector.get_rvine_matrix();
        pair_copulas_ = selector.get_pair_copulas();
    }
        
    //! automatically selects all pair-copula families and fits all parameters.
    //!
    //! @param data nxd matrix of copula data.
    //! @param controls the controls to the algorithm (see FitControlsVinecop).
    void Vinecop::select_families(const Eigen::MatrixXd& data,
                                  FitControlsVinecop controls)
    {
        check_data_dim(data);
        tools_select::FamilySelector selector(data, vine_matrix_, controls);
        if (controls.needs_sparse_select()) {
            selector.sparse_select_all_trees(data);
        } else {
            selector.select_all_trees(data);
        }
        pair_copulas_ = selector.get_pair_copulas();
    }

    //! @name Getters
    //! @{

    //! extracts a pair copula.
    //!
    //! @param tree tree index (starting with 0).
    //! @param edge edge index (starting with 0).
    Bicop Vinecop::get_pair_copula(size_t tree, size_t edge) const
    {
        if (tree > d_ - 2) {
            std::stringstream message;
            message <<
                    "tree index out of bounds" << std::endl <<
                    "allowed: 0, ..., " << d_ - 2 << std::endl <<
                    "actual: " << tree << std::endl;
            throw std::runtime_error(message.str().c_str());
        }
        if (edge > d_ - tree - 2) {
            std::stringstream message;
            message <<
                    "edge index out of bounds" << std::endl <<
                    "allowed: 0, ..., " << d_ - tree - 2 << std::endl <<
                    "actual: " << edge << std::endl <<
                    "tree level: " <<  tree  << std::endl;
            throw std::runtime_error(message.str().c_str());
        }
        return pair_copulas_[tree][edge];
    }
    
    //! extracts all pair copulas.
    //! 
    //! @return a nested std::vector with entry `[t][e]` corresponding to
    //! edge `e` in tree `t`.
    std::vector<std::vector<Bicop>> Vinecop::get_all_pair_copulas() const
    {
        return pair_copulas_;
    }
    
    //! extracts the family of a pair copula.
    //!
    //! @param tree tree index (starting with 0).
    //! @param edge edge index (starting with 0).
    BicopFamily Vinecop::get_family(size_t tree, size_t edge) const
    {
        return get_pair_copula(tree, edge).get_family();
    }
    
    //! extracts the families of all pair copulas.
    //!
    //! @return a nested std::vector with entry `[t][e]` corresponding to
    //! edge `e` in tree `t`.
    std::vector<std::vector<BicopFamily>> Vinecop::get_all_families() const
    {
        std::vector<std::vector<BicopFamily>> families(d_ - 1);
        for (size_t t = 0; t < d_ - 1; ++t)
            families[t].resize(d_ - 1 - t);
        for (size_t tree = 0; tree < d_ - 1; ++tree) {
            for (size_t edge = 0; edge < d_ - 1 - tree; ++edge) {
                families[tree][edge] = get_family(tree, edge);
            }
        }
    
        return families;
    }
    
    //! extracts the rotation of a pair copula.
    //!
    //! @param tree tree index (starting with 0).
    //! @param edge edge index (starting with 0).
    int Vinecop::get_rotation(size_t tree, size_t edge) const
    {
        return get_pair_copula(tree, edge).get_rotation();
    }
    
    //! extracts the rotations of all pair copulas.
    //!
    //! @return a nested std::vector with entry `[t][e]` corresponding to
    //! edge `e` in tree `t`.
    std::vector<std::vector<int>> Vinecop::get_all_rotations() const
    {
        std::vector<std::vector<int>> rotations(d_ - 1);
        for (size_t t = 0; t < d_ - 1; ++t)
            rotations[t].resize(d_ - 1 - t);
        for (size_t tree = 0; tree < d_ - 1; ++tree) {
            for (size_t edge = 0; edge < d_ - 1 - tree; ++edge) {
                rotations[tree][edge] = get_rotation(tree, edge);
            }
        }
    
        return rotations;
    }
    
    //! extracts the parameters of a pair copula.
    //!
    //! @param tree tree index (starting with 0).
    //! @param edge edge index (starting with 0).
    Eigen::VectorXd Vinecop::get_parameters(size_t tree, size_t edge) const
    {
        return get_pair_copula(tree, edge).get_parameters();
    }
    
    //! extracts the parameters of all pair copulas.
    //!
    //! @return a nested std::vector with entry `[t][e]` corresponding to
    //! edge `e` in tree `t`.
    std::vector<std::vector<Eigen::VectorXd>> Vinecop::get_all_parameters() const
    {
        std::vector<std::vector<Eigen::VectorXd>> parameters(d_ - 1);
        for (size_t t = 0; t < d_ - 1; ++t)
            parameters[t].resize(d_ - 1 - t);
        for (size_t tree = 0; tree < d_ - 1; ++tree) {
            for (size_t edge = 0; edge < d_ - 1 - tree; ++edge) {
                parameters[tree][edge] = get_parameters(tree, edge);
            }
        }
    
        return parameters;
    }
    
    //! extracts the structure matrix of the vine copula model.
    Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> Vinecop::get_matrix() const
    {
        return vine_matrix_.get_matrix();
    } 
    
    //! @}
    
    //! calculates the density function of the vine copula model.
    //!
    //! @param u \f$ n \times d \f$ matrix of evaluation points.
    Eigen::VectorXd Vinecop::pdf(const Eigen::MatrixXd& u)
    {
        size_t d = u.cols();
        size_t n = u.rows();
        if (d != d_) {
            std::stringstream message;
            message << "u has wrong number of columns. " <<
                    "expected: " << d_ <<
                    ", actual: " << d << std::endl;
            throw std::runtime_error(message.str().c_str());
        }
    
        // info about the vine structure (reverse rows (!) for more natural indexing)
        Eigen::Matrix<size_t, Eigen::Dynamic, 1> revorder = vine_matrix_.get_order().reverse();
        auto no_matrix  = vine_matrix_.in_natural_order();
        auto max_matrix = vine_matrix_.get_max_matrix();
        MatrixXb needed_hfunc1 = vine_matrix_.get_needed_hfunc1();
        MatrixXb needed_hfunc2 = vine_matrix_.get_needed_hfunc2();
    
        // initial value must be 1.0 for multiplication
        Eigen::VectorXd vine_density = Eigen::VectorXd::Constant(u.rows(), 1.0);
    
        // temporary storage objects for h-functions
        Eigen::MatrixXd hfunc1(n, d);
        Eigen::MatrixXd hfunc2(n, d);
        Eigen::MatrixXd u_e(n, 2);
    
        // fill first row of hfunc2 matrix with evaluation points;
        // points have to be reordered to correspond to natural order
        for (size_t j = 0; j < d; ++j)
            hfunc2.col(j) = u.col(revorder(j) - 1);
    
        for (size_t tree = 0; tree < d - 1; ++tree) {
            for (size_t edge = 0; edge < d - tree - 1; ++edge) {
                // extract evaluation point from hfunction matrices (have been
                // computed in previous tree level)
                size_t m = max_matrix(tree, edge);
                u_e.col(0) = hfunc2.col(edge);
                if (m == no_matrix(tree, edge)) {
                    u_e.col(1) = hfunc2.col(d - m);
                } else {
                    u_e.col(1) = hfunc1.col(d - m);
                }
                
                Bicop edge_copula = get_pair_copula(tree, edge);
                vine_density = vine_density.cwiseProduct(edge_copula.pdf(u_e));
                
                // h-functions are only evaluated if needed in next step
                if (needed_hfunc1(tree + 1, edge)) {
                    hfunc1.col(edge) = edge_copula.hfunc1(u_e);
                }
                if (needed_hfunc2(tree + 1, edge)) {
                    hfunc2.col(edge) = edge_copula.hfunc2(u_e);
                }
            }
        }
    
        return vine_density;
    }

    //! calculates the cumulative distribution of the vine copula model.
    //!
    //! @param u \f$ n \times d \f$ matrix of evaluation points.
    //! @param N integer for the number of quasi-random numbers to draw
    //! to evaluate the distribution (default: 1e4).
    Eigen::VectorXd Vinecop::cdf(const Eigen::MatrixXd& u, const size_t N)
    {
        if (d_ > 360) {
            std::stringstream message;
            message << "cumulative distribution available for models of " <<
                    "dimension 360 or less. This model's dimension: " << d_
                    << std::endl;
            throw std::runtime_error(message.str().c_str());
        }

        size_t d = u.cols();
        size_t n = u.rows();
        if (d != d_) {
            std::stringstream message;
            message << "u has wrong number of columns. " <<
                    "expected: " << d_ <<
                    ", actual: " << d << std::endl;
            throw std::runtime_error(message.str().c_str());
        }

        // Simulate N quasi-random numbers from the vine model
        auto U = tools_stats::ghalton(N, d);
        U = inverse_rosenblatt(U);

        // Alternative: simulate N pseudo-random numbers from the vine model
        //auto U = simulate(N);

        Eigen::VectorXd vine_distribution(n);
        Eigen::ArrayXXd x(N,1);
        Eigen::RowVectorXd temp(d);
        for (size_t i = 0; i < n; i++) {
            temp = u.block(i,0,1,d);
            x = (U.rowwise() - temp).rowwise().maxCoeff().array();
            vine_distribution(i) = (x <= 0.0).count();
        }

        return vine_distribution/((double) N);
    }
    
    //! simulates from a vine copula model, see inverse_rosenblatt().
    //! 
    //! @param n number of observations.
    Eigen::MatrixXd Vinecop::simulate(size_t n)
    {
        Eigen::MatrixXd U = tools_stats::simulate_uniform(n, d_);
        return inverse_rosenblatt(U);
    }

    //! calculates the log-likelihood.
    //!
    //! The log-likelihood is defined as
    //! \f[ \mathrm{loglik} = \sum_{i = 1}^n \ln c(U_{1, i}, ..., U_{d, i}), \f]
    //! where \f$ c \f$ is the copula density pdf().
    //!
    //! @param u \f$n \times d\f$ matrix of observations.
    double Vinecop::loglik(const Eigen::MatrixXd& u)
    {
        return pdf(u).array().log().sum();
    }

    //! calculates the Akaike information criterion (AIC).
    //!
    //! The AIC is defined as
    //! \f[ \mathrm{AIC} = -2\, \mathrm{loglik} + 2 p, \f]
    //! where \f$ \mathrm{loglik} \f$ is the log-liklihood and \f$ p \f$ is the
    //! (effective) number of parameters of the model, see loglik() and
    //! calculate_npars(). The AIC is a consistent model selection criterion
    //! for nonparametric models.
    //!
    //! @param u \f$n \times 2\f$ matrix of observations.
    double Vinecop::aic(const Eigen::MatrixXd& u)
    {
        return -2 * loglik(u) + 2 * calculate_npars();
    }

    //! calculates the Bayesian information criterion (BIC).
    //!
    //! The BIC is defined as
    //! \f[ \mathrm{BIC} = -2\, \mathrm{loglik} +  \ln(n) p, \f]
    //! where \f$ \mathrm{loglik} \f$ is the log-liklihood and \f$ p \f$ is the
    //! (effective) number of parameters of the model, see loglik() and
    //! calculate_npars(). The BIC is a consistent model selection criterion
    //! for nonparametric models.
    //!
    //! @param u \f$n \times 2\f$ matrix of observations.
    double Vinecop::bic(const Eigen::MatrixXd& u)
    {
        return -2 * loglik(u) + calculate_npars() * log(u.rows());
    }

    //! calculates the effective number of parameters.
    //!
    //! Returns sum of the number of parameters for all pair copulas (see
    //! Bicop::calculate_npars()).
    double Vinecop::calculate_npars()
    {
        double npars = 0.0;
        for (auto& tree : pair_copulas_) {
            for (auto& pc : tree) {
                npars += pc.calculate_npars();
            }
        }
        return npars;
    }
    
    
    //! calculates the inverse Rosenblatt transform for a vine copula model.
    //! 
    //! The inverse Rosenblatt transform can be used for simulation: the 
    //! function applied to independent uniform variates resembles simulated
    //! data from the vine copula model.
    //! 
    //! If the problem is too large, it is split recursively into halves (w.r.t.
    //! n, the number of observations).
    //! "Too large" means that the required memory will exceed 1 GB. An 
    //! examplary configuration requiring less than 1 GB is \f$ n = 1000 \f$, 
    //! \f$d = 200\f$.
    //! 
    //! @param u \f$ n \times d \f$ matrix of evaluation points.
    Eigen::MatrixXd Vinecop::inverse_rosenblatt(const Eigen::MatrixXd& u)
    {
        size_t n = u.rows();
        if (n < 1) {
            throw std::runtime_error("n must be at least one");
        }
        size_t d = u.cols();
        if (d != d_) {
            std::stringstream message;
            message << "U has wrong number of columns; " <<
                    "expected: " << d_ <<
                    ", actual: " << d << std::endl;
            throw std::runtime_error(message.str().c_str());
        }
    
        Eigen::MatrixXd U_vine = u;  // output matrix
    
        //                   (direct + indirect)    (U_vine)       (info matrices)
        size_t bytes_required = (8 * 2 * n * d * d) +  (8 * n * d)  + (4 * 4 * d * d);
        // if the problem is too large (requires more than 1 GB memory), split
        // the data into two halves and call simulate on the reduced data.
        if ((n > 1) & (bytes_required > 1e9)) {
            size_t n_half = n / 2;
            size_t n_left = n - n_half;
            U_vine.block(0, 0, n_half, d) = 
                inverse_rosenblatt(u.block(0, 0, n_half, d));
            U_vine.block(n_half, 0, n_left, d) = 
                inverse_rosenblatt(u.block(n_half, 0, n_left, d));
            return U_vine;
        }

        if (d > 2) {
            // info about the vine structure (in upper triangular matrix notation)
            Eigen::Matrix<size_t, Eigen::Dynamic, 1> revorder = vine_matrix_.get_order().reverse();
            auto no_matrix  = vine_matrix_.in_natural_order();
            auto max_matrix = vine_matrix_.get_max_matrix();
            MatrixXb needed_hfunc1 = vine_matrix_.get_needed_hfunc1();
            MatrixXb needed_hfunc2 = vine_matrix_.get_needed_hfunc2();

            // temporary storage objects for (inverse) h-functions
            Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic> hinv2(d, d);
            Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic> hfunc1(d, d);

            // initialize with independent uniforms (corresponding to natural order)
            for (size_t j = 0; j < d; ++j)
                hinv2(d - j - 1, j) = u.col(revorder(j) - 1);
            hfunc1(0, d - 1) = hinv2(0, d - 1);

            // loop through variables (0 is just the inital uniform)
            for (ptrdiff_t var = d - 2; var >= 0; --var) {
                for (ptrdiff_t tree = d - var - 2; tree >= 0; --tree) {
                    Bicop edge_copula = get_pair_copula(tree, var);

                    // extract data for conditional pair
                    Eigen::MatrixXd U_e(n, 2);
                    size_t m = max_matrix(tree, var);
                    U_e.col(0) = hinv2(tree + 1, var);
                    if (m == no_matrix(tree, var)) {
                        U_e.col(1) = hinv2(tree, d - m);
                    } else {
                        U_e.col(1) = hfunc1(tree, d - m);
                    }

                    // inverse Rosenblatt transform simulates data for conditional pair
                    hinv2(tree, var) = edge_copula.hinv2(U_e);

                    // if required at later stage, also calculate hfunc2
                    if (var < (ptrdiff_t) d_ - 1) {
                        if (needed_hfunc1(tree + 1, var)) {
                            U_e.col(0) = hinv2(tree, var);
                            hfunc1(tree + 1, var) = edge_copula.hfunc1(U_e);
                        }
                    }
                }
            }

            // go back to original order
            auto inverse_order = inverse_permutation(revorder);
            for (size_t j = 0; j < d; ++j)
                U_vine.col(j) = hinv2(0, inverse_order(j));
        } else {
            U_vine.col(1) = get_pair_copula(0, 0).hinv1(u);
        }
    
        return U_vine;
    }
    
    //! checks if dimension d of the data matches the dimension of the vine.
    void Vinecop::check_data_dim(const Eigen::MatrixXd& data)
    {
        size_t d_data = data.cols();
        if (d_data != d_) {
            std::stringstream message;
            message << "data has wrong number of columns; " <<
                    "expected: " << d_ <<
                    ", actual: " << d_data << std::endl;
            throw std::runtime_error(message.str().c_str());
        }
    }  
    
          
    // get indexes for reverting back to old order in simulation routine
    Eigen::Matrix<size_t, Eigen::Dynamic, 1> Vinecop::inverse_permutation(
        const Eigen::Matrix<size_t, Eigen::Dynamic, 1>& order) {
        // start with (0, 1, .., k)
        auto indexes = tools_stl::seq_int(0, order.size());
    
        // get sort indexes by comparing values in order
        std::sort(indexes.begin(), indexes.end(),
                  [&order](size_t i1, size_t i2) {return order(i1) < order(i2);});
    
        // convert to Eigen::Matrix<size_t, Eigen::Dynamic, 1>;
        return Eigen::Map<Eigen::Matrix<size_t, Eigen::Dynamic, 1>>(&indexes[0], order.size());
    }

}
