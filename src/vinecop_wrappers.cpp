#include <RcppEigen.h>
#include "bicop_wrappers.hpp"

using namespace vinecopulib;

// [[Rcpp::export()]]
void rvine_matrix_check_cpp(Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> matrix) {
    auto rvm = RVineMatrix(matrix);
}

Vinecop vinecop_wrap(const Rcpp::List& vinecop_r)
{
    Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> matrix = vinecop_r["matrix"];
    size_t d = matrix.cols();
    
    Rcpp::List pair_copulas = vinecop_r["pair_copulas"];
    auto pc_store = Vinecop::make_pair_copula_store(d, pair_copulas.size());
    Rcpp::List tree_pcs, pc;
    for (size_t t = 0; t < pc_store.size(); ++t) {
        tree_pcs = pair_copulas[t];
        if (tree_pcs.size() != d - 1 - t) {
            throw std::runtime_error("length(pair_copulas[[t]]) must be d-t");
        }
        for(size_t e = 0; e < d - 1 - t; ++e) {
            pc = tree_pcs[e];
            pc_store[t][e] = bicop_wrap(pc);
        }
    }
    
    // omit R-vine matrix check, already done in R
    return Vinecop(pc_store, matrix, false);
}


Rcpp::List vinecop_wrap(const Vinecop& vinecop_cpp) {
    auto matrix = vinecop_cpp.get_matrix();
    auto pcs = vinecop_cpp.get_all_pair_copulas();
    size_t d = matrix.cols();
    
    Rcpp::List pair_copulas(pcs.size());

    for (size_t t = 0; t < pcs.size(); ++t) {
        Rcpp::List tree_pcs(d - 1 - t);
        for(size_t e = 0; e < d - 1 - t; ++e) {
            tree_pcs[e] = bicop_wrap(pcs[t][e]);
        }
        pair_copulas[t] = tree_pcs;
    }
    double npars = vinecop_cpp.calculate_npars();
    double threshold = vinecop_cpp.get_threshold();
    return Rcpp::List::create(
        Rcpp::Named("pair_copulas") = pair_copulas,
        Rcpp::Named("matrix") = matrix,
        Rcpp::Named("npars") = npars,
        Rcpp::Named("threshold") = threshold
    );
}

// [[Rcpp::export()]]
void vinecop_check_cpp(Rcpp::List vinecop_r) {
    vinecop_wrap(vinecop_r);
}


// [[Rcpp::export()]]
Eigen::MatrixXd vinecop_inverse_rosenblatt_cpp(const Eigen::MatrixXd& U,
                                               const Rcpp::List& vinecop_r)
{
    return vinecop_wrap(vinecop_r).inverse_rosenblatt(U);
}

// [[Rcpp::export()]]
Eigen::VectorXd vinecop_pdf_cpp(const Eigen::MatrixXd& u, const Rcpp::List& vinecop_r)
{
    return vinecop_wrap(vinecop_r).pdf(u);
}

// [[Rcpp::export()]]
Eigen::VectorXd vinecop_cdf_cpp(const Eigen::MatrixXd& u, const Rcpp::List& vinecop_r, size_t N)
{
    return vinecop_wrap(vinecop_r).cdf(u, N);
}

// [[Rcpp::export()]]
double vinecop_loglik_cpp(const Eigen::MatrixXd& u, const Rcpp::List& vinecop_r)
{
    return vinecop_wrap(vinecop_r).loglik(u);
}

// [[Rcpp::export()]]
Rcpp::List vinecop_select_cpp(
        const Eigen::MatrixXd& data,
        Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> matrix,
        std::vector<std::string> family_set,
        std::string par_method,
        std::string nonpar_method,
        double mult,
        int truncation_level,
        std::string tree_criterion,
        double threshold,
        std::string selection_criterion,
        bool select_truncation_level,
        bool select_threshold,
        bool preselect_families,
        bool show_trace,
        size_t num_threads
)
{
    std::vector<BicopFamily> fam_set(family_set.size());
    for (unsigned int fam = 0; fam < fam_set.size(); ++fam) {
        fam_set[fam] = to_cpp_family(family_set[fam]);
    }

    FitControlsVinecop fit_controls(
            fam_set,
            par_method,
            nonpar_method,
            mult,
            truncation_level,
            tree_criterion,
            threshold,
            selection_criterion,
            preselect_families,
            select_truncation_level,
            select_threshold,
            show_trace,
            num_threads
    );
    Vinecop vinecop_cpp(data.cols());
    if (matrix.cols() > 1) {
        vinecop_cpp = Vinecop(matrix);
        vinecop_cpp.select_families(data, fit_controls);
    } else {
        vinecop_cpp.select_all(data, fit_controls);
    }

    return vinecop_wrap(vinecop_cpp);
}
