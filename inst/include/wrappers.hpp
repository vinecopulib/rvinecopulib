#define BOOST_NO_AUTO_PTR
#define BOOST_MATH_PROMOTE_DOUBLE_POLICY false

#include <RcppEigen.h>
#include "vinecopulib.hpp"

using namespace vinecopulib;

// tools exports

Eigen::MatrixXd pseudo_obs_cpp(Eigen::MatrixXd x, std::string ties_method);

// bicop wrappers

BicopFamily to_cpp_family(const std::string& fam);

std::string to_r_family(const BicopFamily& fam);

Bicop bicop_wrap(const Rcpp::List& bicop_r);

Rcpp::List bicop_wrap(Bicop bicop_cpp, bool is_fitted);

// bicop exports

void bicop_check_cpp(const Rcpp::List& bicop_r);

Rcpp::List bicop_select_cpp(const Eigen::MatrixXd& data,
                            std::vector<std::string> family_set,
                            std::string par_method,
                            std::string nonpar_method,
                            double mult,
                            std::string selcrit,
                            const Eigen::VectorXd& weights,
                            double psi0,
                            bool presel,
                            size_t num_threads);

Eigen::VectorXd bicop_pdf_cpp(const Eigen::MatrixXd& u, 
                              const Rcpp::List& bicop_r);

Eigen::VectorXd bicop_cdf_cpp(const Eigen::MatrixXd& u, 
                              const Rcpp::List& bicop_r);

Eigen::VectorXd bicop_hfunc1_cpp(const Eigen::MatrixXd& u, 
                                 const Rcpp::List& bicop_r);

Eigen::VectorXd bicop_hfunc2_cpp(const Eigen::MatrixXd& u, 
                                 const Rcpp::List& bicop_r);

Eigen::VectorXd bicop_hinv1_cpp(const Eigen::MatrixXd& u, 
                                const Rcpp::List& bicop_r);

Eigen::VectorXd bicop_hinv2_cpp(const Eigen::MatrixXd& u, 
                                const Rcpp::List& bicop_r);

Eigen::MatrixXd bicop_sim_cpp(const Rcpp::List& bicop_r, const size_t &n, 
                              const bool qrng,
                              std::vector<int> seeds);

double bicop_loglik_cpp(Eigen::MatrixXd& u, const Rcpp::List& bicop_r);

double bicop_par_to_tau_cpp(const Rcpp::List& bicop_r);

Eigen::MatrixXd bicop_tau_to_par_cpp(const Rcpp::List& bicop_r, 
                                     const double& tau);

// structure wrappers

TriangularArray<size_t> struct_array_wrap(const Rcpp::List& struct_array_r,
                                          size_t trunc_lvl);

Rcpp::List struct_array_wrap(const TriangularArray<size_t>& struct_array);

RVineStructure rvine_structure_wrap(const Rcpp::List& rvine_structure_r,
                                    bool check = false,
                                    bool is_natural_order = true);

Rcpp::List rvine_structure_wrap(const RVineStructure& rvine_struct);

// structure exports

void rvine_structure_check_cpp(const Rcpp::List& rvine_struct,
                               bool is_natural_order = true);

void rvine_matrix_check_cpp(Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> matrix);

// vinecop wrappers

std::vector<std::vector<Bicop>> pair_copulas_wrap(const Rcpp::List& pair_copulas_r,
                                                  size_t d);

Rcpp::List pair_copulas_wrap(std::vector<std::vector<Bicop>> pair_copulas,
                             size_t d,
                             bool is_fitted);

Vinecop vinecop_wrap(const Rcpp::List& vinecop_r);


Rcpp::List vinecop_wrap(const Vinecop& vinecop_cpp, bool is_fitted = FALSE);

// vinecop exports

void vinecop_check_cpp(Rcpp::List vinecop_r);

Eigen::MatrixXd vinecop_inverse_rosenblatt_cpp(const Eigen::MatrixXd& U,
                                               const Rcpp::List& vinecop_r,
                                               size_t cores);

Eigen::MatrixXd vinecop_rosenblatt_cpp(const Eigen::MatrixXd& U,
                                       const Rcpp::List& vinecop_r,
                                       size_t cores);

Eigen::MatrixXd vinecop_sim_cpp(const Rcpp::List& vinecop_r, 
                                const size_t n, 
                                const bool qrng,
                                size_t cores,
                                std::vector<int> seeds);

Eigen::VectorXd vinecop_pdf_cpp(const Eigen::MatrixXd& u, 
                                const Rcpp::List& vinecop_r,
                                size_t cores);

Eigen::VectorXd vinecop_cdf_cpp(const Eigen::MatrixXd& u, 
                                const Rcpp::List& vinecop_r, 
                                size_t N,
                                size_t cores,
                                std::vector<int> seeds);

double vinecop_loglik_cpp(const Eigen::MatrixXd& u, 
                          const Rcpp::List& vinecop_r,
                          size_t cores);

double vinecop_mbicv_cpp(const Eigen::MatrixXd& u, 
                         const Rcpp::List& vinecop_r,
                         double psi0,
                         size_t cores);

Rcpp::List vinecop_select_cpp(const Eigen::MatrixXd& data,
                              bool is_structure_provided,
                              Rcpp::List& structure,
                              std::vector<std::string> family_set,
                              std::string par_method,
                              std::string nonpar_method,
                              double mult,
                              int truncation_level,
                              std::string tree_criterion,
                              double threshold,
                              std::string selection_criterion,
                              const Eigen::VectorXd& weights,
                              double psi0,
                              bool select_truncation_level,
                              bool select_threshold,
                              bool preselect_families,
                              bool show_trace,
                              size_t num_threads);
