#include <RcppEigen.h>
#include "Bicop_wrappers.hpp"

using namespace vinecopulib;


Vinecop vinecop_wrap(const Rcpp::List& vinecop_r)
{
    Rcpp::List pair_copulas = vinecop_r["pair_copulas"];
    int d = pair_copulas.size() + 1;
    auto pc_store = Vinecop::make_pair_copula_store(d);
    Rcpp::List tree_pcs, pc;
    for (int t = 0; t < d - 1; ++t) {
        for(int e = 0; e < d - 1 - t; ++e) {
            tree_pcs = pair_copulas[t];
            pc = tree_pcs[e];
            pc_store[t][e] = bicop_wrap(pc);
        }
    }
    
    Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> matrix = vinecop_r["matrix"];
    
    return Vinecop(pc_store, matrix);
}


Rcpp::List vinecop_wrap(Vinecop& vinecop_cpp) {
    auto matrix = vinecop_cpp.get_matrix();
    int d = matrix.cols();
    Rcpp::List pair_copulas(d - 1);
    
    for (int t = 0; t < d - 1; ++t) {
        Rcpp::List tree_pcs(d - 1 - t);
        for(int e = 0; e < d - 1 - t; ++e) {
            tree_pcs[e] = bicop_wrap(vinecop_cpp.get_pair_copula(t, e));
        }
        pair_copulas[t] = tree_pcs;
    }
    
    return Rcpp::List::create(
        Rcpp::Named("pair_copulas") = pair_copulas,
        Rcpp::Named("matrix") = matrix
    );
}

// [[Rcpp::export()]]
void vinecop_check_cpp(Rcpp::List vinecop_r) {
    vinecop_wrap(vinecop_r);
}


// [[Rcpp::export()]]
Eigen::MatrixXd vinecop_sim_cpp(int n, const Rcpp::List& vinecop_r)
{
    return vinecop_wrap(vinecop_r).simulate(n);
}

// [[Rcpp::export()]]
Eigen::VectorXd vinecop_pdf_cpp(const Eigen::MatrixXd& u, const Rcpp::List& vinecop_r)
{
    return vinecop_wrap(vinecop_r).pdf(u);
}

// [[Rcpp::export()]]
Rcpp::List vinecop_select_cpp(
        const Eigen::MatrixXd& data,
        Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> matrix,
        std::vector<std::string> family_set,
        std::string method,
        double mult,
        int truncation_level,
        std::string tree_criterion,
        double threshold,
        std::string selection_criterion,
        bool preselect_families
)
{
    std::vector<BicopFamily> fam_set(family_set.size());
    for (unsigned int fam = 0; fam < fam_set.size(); ++fam) {
        fam_set[fam] = to_cpp_family(family_set[fam]);
    }
    
    FitControlsVinecop fit_controls(
            fam_set,
            method,
            mult,
            truncation_level,
            tree_criterion,
            threshold,
            selection_criterion,
            preselect_families,
            false  // show_trace
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
