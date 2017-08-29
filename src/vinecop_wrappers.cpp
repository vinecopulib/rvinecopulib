#include <RcppEigen.h>
#include "bicop_wrappers.hpp"

using namespace vinecopulib;

// [[Rcpp::export()]]
void rvine_matrix_check_cpp(Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> matrix) {
    auto rvm = RVineMatrix(matrix);
}

Vinecop vinecop_wrap(const Rcpp::List& vinecop_r)
{
    Rcpp::List pair_copulas = vinecop_r["pair_copulas"];
    size_t d = pair_copulas.size() + 1;
    auto pc_store = Vinecop::make_pair_copula_store(d);
    Rcpp::List tree_pcs, pc;
    for (size_t t = 0; t < d - 1; ++t) {
        tree_pcs = pair_copulas[t];
        for(size_t e = 0; e < d - 1 - t; ++e) {
            pc = tree_pcs[e];
            pc_store[t][e] = bicop_wrap(pc);
        }
    }

    Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> matrix = vinecop_r["matrix"];
    
    // omit R-vine matrix check, already done in R
    return Vinecop(pc_store, matrix, false);
}


Rcpp::List vinecop_wrap(const Vinecop& vinecop_cpp) {
    auto matrix = vinecop_cpp.get_matrix();
    size_t d = matrix.cols();
    Rcpp::List pair_copulas(d - 1);

    for (size_t t = 0; t < d - 1; ++t) {
        Rcpp::List tree_pcs(d - 1 - t);
        for(size_t e = 0; e < d - 1 - t; ++e) {
            tree_pcs[e] = bicop_wrap(vinecop_cpp.get_pair_copula(t, e));
        }
        pair_copulas[t] = tree_pcs;
    }
    double npars = vinecop_cpp.calculate_npars();
    return Rcpp::List::create(
        Rcpp::Named("pair_copulas") = pair_copulas,
        Rcpp::Named("matrix") = matrix,
        Rcpp::Named("npars") = npars
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
    Progress p(0, false);
    return vinecop_wrap(vinecop_r).inverse_rosenblatt(U);
}

// [[Rcpp::export()]]
Eigen::VectorXd vinecop_pdf_cpp(const Eigen::MatrixXd& u, const Rcpp::List& vinecop_r)
{
    Progress p(0, false);
    return vinecop_wrap(vinecop_r).pdf(u);
}

// [[Rcpp::export()]]
Eigen::VectorXd vinecop_cdf_cpp(const Eigen::MatrixXd& u, const Rcpp::List& vinecop_r, size_t N)
{
    Progress p(0, false);
    return vinecop_wrap(vinecop_r).cdf(u, N);
}

// [[Rcpp::export()]]
double vinecop_loglik_cpp(const Eigen::MatrixXd& u, const Rcpp::List& vinecop_r)
{
    Progress p(0, false);
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
        bool show_trace
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
            show_trace
    );
    Vinecop vinecop_cpp(data.cols());
    Progress p(0, false);
    if (matrix.cols() > 1) {
        vinecop_cpp = Vinecop(matrix);
        vinecop_cpp.select_families(data, fit_controls);
    } else {
        vinecop_cpp.select_all(data, fit_controls);
    }

    return vinecop_wrap(vinecop_cpp);
}
