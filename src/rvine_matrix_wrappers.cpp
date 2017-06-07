#include <RcppEigen.h>
#include <vinecopulib/vinecop/rvine_matrix.hpp>

using namespace vinecopulib;

// Vinecop rvine_matrix_wrap(const Rcpp::List& vinecop_r)
// {
//     Rcpp::List pair_copulas = vinecop_r["pair_copulas"];
//     size_t d = pair_copulas.size() + 1;
//     auto pc_store = Vinecop::make_pair_copula_store(d);
//     Rcpp::List tree_pcs, pc;
//     for (size_t t = 0; t < d - 1; ++t) {
//         tree_pcs = pair_copulas[t];
//         for(size_t e = 0; e < d - 1 - t; ++e) {
//             pc = tree_pcs[e];
//             pc_store[t][e] = bicop_wrap(pc);
//         }
//     }
// 
//     Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> matrix = vinecop_r["matrix"];
// 
//     return Vinecop(pc_store, matrix);
// }
// 
// 
// Rcpp::List vinecop_wrap(const Vinecop& vinecop_cpp) {
//     auto matrix = vinecop_cpp.get_matrix();
//     size_t d = matrix.cols();
//     Rcpp::List pair_copulas(d - 1);
// 
//     for (size_t t = 0; t < d - 1; ++t) {
//         Rcpp::List tree_pcs(d - 1 - t);
//         for(size_t e = 0; e < d - 1 - t; ++e) {
//             tree_pcs[e] = bicop_wrap(vinecop_cpp.get_pair_copula(t, e));
//         }
//         pair_copulas[t] = tree_pcs;
//     }
// 
//     return Rcpp::List::create(
//         Rcpp::Named("pair_copulas") = pair_copulas,
//         Rcpp::Named("matrix") = matrix
//     );
// }

// [[Rcpp::export()]]
void rvine_matrix_check_cpp(Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> matrix) {
    auto rvm = RVineMatrix(matrix);
}

// [[Rcpp::export()]]
Eigen::Matrix<size_t, Eigen::Dynamic, 1> 
rvine_matrix_get_order_cpp(Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> matrix) {
    auto rvm = RVineMatrix(matrix, false);  // don't check validity again
    return rvm.get_order();
}
