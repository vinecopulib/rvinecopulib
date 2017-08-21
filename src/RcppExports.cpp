// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// bicop_check_cpp
void bicop_check_cpp(const Rcpp::List& bicop_r);
RcppExport SEXP _rvinecopulib_bicop_check_cpp(SEXP bicop_rSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type bicop_r(bicop_rSEXP);
    bicop_check_cpp(bicop_r);
    return R_NilValue;
END_RCPP
}
// bicop_select_cpp
Rcpp::List bicop_select_cpp(Eigen::MatrixXd& data, std::vector<std::string> family_set, std::string par_method, std::string nonpar_method, double mult, std::string selcrit, bool presel);
RcppExport SEXP _rvinecopulib_bicop_select_cpp(SEXP dataSEXP, SEXP family_setSEXP, SEXP par_methodSEXP, SEXP nonpar_methodSEXP, SEXP multSEXP, SEXP selcritSEXP, SEXP preselSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type family_set(family_setSEXP);
    Rcpp::traits::input_parameter< std::string >::type par_method(par_methodSEXP);
    Rcpp::traits::input_parameter< std::string >::type nonpar_method(nonpar_methodSEXP);
    Rcpp::traits::input_parameter< double >::type mult(multSEXP);
    Rcpp::traits::input_parameter< std::string >::type selcrit(selcritSEXP);
    Rcpp::traits::input_parameter< bool >::type presel(preselSEXP);
    rcpp_result_gen = Rcpp::wrap(bicop_select_cpp(data, family_set, par_method, nonpar_method, mult, selcrit, presel));
    return rcpp_result_gen;
END_RCPP
}
// bicop_pdf_cpp
Eigen::VectorXd bicop_pdf_cpp(const Eigen::MatrixXd& u, const Rcpp::List& bicop_r);
RcppExport SEXP _rvinecopulib_bicop_pdf_cpp(SEXP uSEXP, SEXP bicop_rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type bicop_r(bicop_rSEXP);
    rcpp_result_gen = Rcpp::wrap(bicop_pdf_cpp(u, bicop_r));
    return rcpp_result_gen;
END_RCPP
}
// bicop_cdf_cpp
Eigen::VectorXd bicop_cdf_cpp(const Eigen::MatrixXd& u, const Rcpp::List& bicop_r);
RcppExport SEXP _rvinecopulib_bicop_cdf_cpp(SEXP uSEXP, SEXP bicop_rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type bicop_r(bicop_rSEXP);
    rcpp_result_gen = Rcpp::wrap(bicop_cdf_cpp(u, bicop_r));
    return rcpp_result_gen;
END_RCPP
}
// bicop_hfunc1_cpp
Eigen::VectorXd bicop_hfunc1_cpp(const Eigen::MatrixXd& u, const Rcpp::List& bicop_r);
RcppExport SEXP _rvinecopulib_bicop_hfunc1_cpp(SEXP uSEXP, SEXP bicop_rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type bicop_r(bicop_rSEXP);
    rcpp_result_gen = Rcpp::wrap(bicop_hfunc1_cpp(u, bicop_r));
    return rcpp_result_gen;
END_RCPP
}
// bicop_hfunc2_cpp
Eigen::VectorXd bicop_hfunc2_cpp(const Eigen::MatrixXd& u, const Rcpp::List& bicop_r);
RcppExport SEXP _rvinecopulib_bicop_hfunc2_cpp(SEXP uSEXP, SEXP bicop_rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type bicop_r(bicop_rSEXP);
    rcpp_result_gen = Rcpp::wrap(bicop_hfunc2_cpp(u, bicop_r));
    return rcpp_result_gen;
END_RCPP
}
// bicop_hinv1_cpp
Eigen::VectorXd bicop_hinv1_cpp(const Eigen::MatrixXd& u, const Rcpp::List& bicop_r);
RcppExport SEXP _rvinecopulib_bicop_hinv1_cpp(SEXP uSEXP, SEXP bicop_rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type bicop_r(bicop_rSEXP);
    rcpp_result_gen = Rcpp::wrap(bicop_hinv1_cpp(u, bicop_r));
    return rcpp_result_gen;
END_RCPP
}
// bicop_hinv2_cpp
Eigen::VectorXd bicop_hinv2_cpp(const Eigen::MatrixXd& u, const Rcpp::List& bicop_r);
RcppExport SEXP _rvinecopulib_bicop_hinv2_cpp(SEXP uSEXP, SEXP bicop_rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type bicop_r(bicop_rSEXP);
    rcpp_result_gen = Rcpp::wrap(bicop_hinv2_cpp(u, bicop_r));
    return rcpp_result_gen;
END_RCPP
}
// bicop_simulate_cpp
Eigen::MatrixXd bicop_simulate_cpp(int n, const Rcpp::List& bicop_r);
RcppExport SEXP _rvinecopulib_bicop_simulate_cpp(SEXP nSEXP, SEXP bicop_rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type bicop_r(bicop_rSEXP);
    rcpp_result_gen = Rcpp::wrap(bicop_simulate_cpp(n, bicop_r));
    return rcpp_result_gen;
END_RCPP
}
// bicop_loglik_cpp
double bicop_loglik_cpp(Eigen::MatrixXd& u, const Rcpp::List& bicop_r);
RcppExport SEXP _rvinecopulib_bicop_loglik_cpp(SEXP uSEXP, SEXP bicop_rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type bicop_r(bicop_rSEXP);
    rcpp_result_gen = Rcpp::wrap(bicop_loglik_cpp(u, bicop_r));
    return rcpp_result_gen;
END_RCPP
}
// bicop_aic_cpp
double bicop_aic_cpp(Eigen::MatrixXd& u, const Rcpp::List& bicop_r);
RcppExport SEXP _rvinecopulib_bicop_aic_cpp(SEXP uSEXP, SEXP bicop_rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type bicop_r(bicop_rSEXP);
    rcpp_result_gen = Rcpp::wrap(bicop_aic_cpp(u, bicop_r));
    return rcpp_result_gen;
END_RCPP
}
// bicop_bic_cpp
double bicop_bic_cpp(Eigen::MatrixXd& u, const Rcpp::List& bicop_r);
RcppExport SEXP _rvinecopulib_bicop_bic_cpp(SEXP uSEXP, SEXP bicop_rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type bicop_r(bicop_rSEXP);
    rcpp_result_gen = Rcpp::wrap(bicop_bic_cpp(u, bicop_r));
    return rcpp_result_gen;
END_RCPP
}
// bicop_par_to_tau_cpp
double bicop_par_to_tau_cpp(const Rcpp::List& bicop_r);
RcppExport SEXP _rvinecopulib_bicop_par_to_tau_cpp(SEXP bicop_rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type bicop_r(bicop_rSEXP);
    rcpp_result_gen = Rcpp::wrap(bicop_par_to_tau_cpp(bicop_r));
    return rcpp_result_gen;
END_RCPP
}
// bicop_tau_to_par_cpp
Eigen::MatrixXd bicop_tau_to_par_cpp(const Rcpp::List& bicop_r, const double& tau);
RcppExport SEXP _rvinecopulib_bicop_tau_to_par_cpp(SEXP bicop_rSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type bicop_r(bicop_rSEXP);
    Rcpp::traits::input_parameter< const double& >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(bicop_tau_to_par_cpp(bicop_r, tau));
    return rcpp_result_gen;
END_RCPP
}
// rvine_matrix_check_cpp
void rvine_matrix_check_cpp(Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> matrix);
RcppExport SEXP _rvinecopulib_rvine_matrix_check_cpp(SEXP matrixSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> >::type matrix(matrixSEXP);
    rvine_matrix_check_cpp(matrix);
    return R_NilValue;
END_RCPP
}
// vinecop_check_cpp
void vinecop_check_cpp(Rcpp::List vinecop_r);
RcppExport SEXP _rvinecopulib_vinecop_check_cpp(SEXP vinecop_rSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type vinecop_r(vinecop_rSEXP);
    vinecop_check_cpp(vinecop_r);
    return R_NilValue;
END_RCPP
}
// vinecop_sim_cpp
Eigen::MatrixXd vinecop_sim_cpp(int n, const Rcpp::List& vinecop_r);
RcppExport SEXP _rvinecopulib_vinecop_sim_cpp(SEXP nSEXP, SEXP vinecop_rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type vinecop_r(vinecop_rSEXP);
    rcpp_result_gen = Rcpp::wrap(vinecop_sim_cpp(n, vinecop_r));
    return rcpp_result_gen;
END_RCPP
}
// vinecop_pdf_cpp
Eigen::VectorXd vinecop_pdf_cpp(const Eigen::MatrixXd& u, const Rcpp::List& vinecop_r);
RcppExport SEXP _rvinecopulib_vinecop_pdf_cpp(SEXP uSEXP, SEXP vinecop_rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type vinecop_r(vinecop_rSEXP);
    rcpp_result_gen = Rcpp::wrap(vinecop_pdf_cpp(u, vinecop_r));
    return rcpp_result_gen;
END_RCPP
}
// vinecop_cdf_cpp
Eigen::VectorXd vinecop_cdf_cpp(const Eigen::MatrixXd& u, const Rcpp::List& vinecop_r, size_t N);
RcppExport SEXP _rvinecopulib_vinecop_cdf_cpp(SEXP uSEXP, SEXP vinecop_rSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type vinecop_r(vinecop_rSEXP);
    Rcpp::traits::input_parameter< size_t >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(vinecop_cdf_cpp(u, vinecop_r, N));
    return rcpp_result_gen;
END_RCPP
}
// vinecop_loglik_cpp
double vinecop_loglik_cpp(const Eigen::MatrixXd& u, const Rcpp::List& vinecop_r);
RcppExport SEXP _rvinecopulib_vinecop_loglik_cpp(SEXP uSEXP, SEXP vinecop_rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type vinecop_r(vinecop_rSEXP);
    rcpp_result_gen = Rcpp::wrap(vinecop_loglik_cpp(u, vinecop_r));
    return rcpp_result_gen;
END_RCPP
}
// vinecop_aic_cpp
double vinecop_aic_cpp(Eigen::MatrixXd& u, const Rcpp::List& vinecop_r);
RcppExport SEXP _rvinecopulib_vinecop_aic_cpp(SEXP uSEXP, SEXP vinecop_rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type vinecop_r(vinecop_rSEXP);
    rcpp_result_gen = Rcpp::wrap(vinecop_aic_cpp(u, vinecop_r));
    return rcpp_result_gen;
END_RCPP
}
// vinecop_bic_cpp
double vinecop_bic_cpp(Eigen::MatrixXd& u, const Rcpp::List& vinecop_r);
RcppExport SEXP _rvinecopulib_vinecop_bic_cpp(SEXP uSEXP, SEXP vinecop_rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type vinecop_r(vinecop_rSEXP);
    rcpp_result_gen = Rcpp::wrap(vinecop_bic_cpp(u, vinecop_r));
    return rcpp_result_gen;
END_RCPP
}
// vinecop_select_cpp
Rcpp::List vinecop_select_cpp(const Eigen::MatrixXd& data, Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> matrix, std::vector<std::string> family_set, std::string par_method, std::string nonpar_method, double mult, int truncation_level, std::string tree_criterion, double threshold, std::string selection_criterion, bool preselect_families);
RcppExport SEXP _rvinecopulib_vinecop_select_cpp(SEXP dataSEXP, SEXP matrixSEXP, SEXP family_setSEXP, SEXP par_methodSEXP, SEXP nonpar_methodSEXP, SEXP multSEXP, SEXP truncation_levelSEXP, SEXP tree_criterionSEXP, SEXP thresholdSEXP, SEXP selection_criterionSEXP, SEXP preselect_familiesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> >::type matrix(matrixSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type family_set(family_setSEXP);
    Rcpp::traits::input_parameter< std::string >::type par_method(par_methodSEXP);
    Rcpp::traits::input_parameter< std::string >::type nonpar_method(nonpar_methodSEXP);
    Rcpp::traits::input_parameter< double >::type mult(multSEXP);
    Rcpp::traits::input_parameter< int >::type truncation_level(truncation_levelSEXP);
    Rcpp::traits::input_parameter< std::string >::type tree_criterion(tree_criterionSEXP);
    Rcpp::traits::input_parameter< double >::type threshold(thresholdSEXP);
    Rcpp::traits::input_parameter< std::string >::type selection_criterion(selection_criterionSEXP);
    Rcpp::traits::input_parameter< bool >::type preselect_families(preselect_familiesSEXP);
    rcpp_result_gen = Rcpp::wrap(vinecop_select_cpp(data, matrix, family_set, par_method, nonpar_method, mult, truncation_level, tree_criterion, threshold, selection_criterion, preselect_families));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rvinecopulib_bicop_check_cpp", (DL_FUNC) &_rvinecopulib_bicop_check_cpp, 1},
    {"_rvinecopulib_bicop_select_cpp", (DL_FUNC) &_rvinecopulib_bicop_select_cpp, 7},
    {"_rvinecopulib_bicop_pdf_cpp", (DL_FUNC) &_rvinecopulib_bicop_pdf_cpp, 2},
    {"_rvinecopulib_bicop_cdf_cpp", (DL_FUNC) &_rvinecopulib_bicop_cdf_cpp, 2},
    {"_rvinecopulib_bicop_hfunc1_cpp", (DL_FUNC) &_rvinecopulib_bicop_hfunc1_cpp, 2},
    {"_rvinecopulib_bicop_hfunc2_cpp", (DL_FUNC) &_rvinecopulib_bicop_hfunc2_cpp, 2},
    {"_rvinecopulib_bicop_hinv1_cpp", (DL_FUNC) &_rvinecopulib_bicop_hinv1_cpp, 2},
    {"_rvinecopulib_bicop_hinv2_cpp", (DL_FUNC) &_rvinecopulib_bicop_hinv2_cpp, 2},
    {"_rvinecopulib_bicop_simulate_cpp", (DL_FUNC) &_rvinecopulib_bicop_simulate_cpp, 2},
    {"_rvinecopulib_bicop_loglik_cpp", (DL_FUNC) &_rvinecopulib_bicop_loglik_cpp, 2},
    {"_rvinecopulib_bicop_aic_cpp", (DL_FUNC) &_rvinecopulib_bicop_aic_cpp, 2},
    {"_rvinecopulib_bicop_bic_cpp", (DL_FUNC) &_rvinecopulib_bicop_bic_cpp, 2},
    {"_rvinecopulib_bicop_par_to_tau_cpp", (DL_FUNC) &_rvinecopulib_bicop_par_to_tau_cpp, 1},
    {"_rvinecopulib_bicop_tau_to_par_cpp", (DL_FUNC) &_rvinecopulib_bicop_tau_to_par_cpp, 2},
    {"_rvinecopulib_rvine_matrix_check_cpp", (DL_FUNC) &_rvinecopulib_rvine_matrix_check_cpp, 1},
    {"_rvinecopulib_vinecop_check_cpp", (DL_FUNC) &_rvinecopulib_vinecop_check_cpp, 1},
    {"_rvinecopulib_vinecop_sim_cpp", (DL_FUNC) &_rvinecopulib_vinecop_sim_cpp, 2},
    {"_rvinecopulib_vinecop_pdf_cpp", (DL_FUNC) &_rvinecopulib_vinecop_pdf_cpp, 2},
    {"_rvinecopulib_vinecop_cdf_cpp", (DL_FUNC) &_rvinecopulib_vinecop_cdf_cpp, 3},
    {"_rvinecopulib_vinecop_loglik_cpp", (DL_FUNC) &_rvinecopulib_vinecop_loglik_cpp, 2},
    {"_rvinecopulib_vinecop_aic_cpp", (DL_FUNC) &_rvinecopulib_vinecop_aic_cpp, 2},
    {"_rvinecopulib_vinecop_bic_cpp", (DL_FUNC) &_rvinecopulib_vinecop_bic_cpp, 2},
    {"_rvinecopulib_vinecop_select_cpp", (DL_FUNC) &_rvinecopulib_vinecop_select_cpp, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_rvinecopulib(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
