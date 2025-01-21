#include "vinecopulib-wrappers.hpp"
#include "kde1d-wrappers.hpp"

using namespace vinecopulib;

// tools exports -------------------------------------------

// [[Rcpp::export]]
Eigen::MatrixXd pseudo_obs_cpp(Eigen::MatrixXd x, std::string ties_method)
{
  return vinecopulib::tools_stats::to_pseudo_obs(x, ties_method);
}

// bicop exports -------------------------------------------

// [[Rcpp::export]]
void bicop_check_cpp(const Rcpp::List& bicop_r)
{
  bicop_wrap(bicop_r);
}

// [[Rcpp::export()]]
Rcpp::List bicop_select_cpp(const Eigen::MatrixXd& data,
                            std::vector<std::string> family_set,
                            std::string par_method,
                            std::string nonpar_method,
                            double mult,
                            std::string selcrit,
                            const Eigen::VectorXd& weights,
                            double psi0,
                            bool presel,
                            size_t num_threads,
                            bool allow_rotations,
                            std::vector<std::string> var_types)
{
  std::vector<BicopFamily> fam_set(family_set.size());
  for (unsigned int fam = 0; fam < fam_set.size(); ++fam) {
    fam_set[fam] = to_cpp_family(family_set[fam]);
  }
  FitControlsBicop controls;
  controls.set_family_set(fam_set);
  controls.set_parametric_method(par_method);
  controls.set_nonparametric_method(nonpar_method);
  controls.set_nonparametric_mult(mult);
  controls.set_selection_criterion(selcrit);
  controls.set_weights(weights);
  controls.set_psi0(psi0);
  controls.set_preselect_families(presel);
  controls.set_allow_rotations(allow_rotations);
  controls.set_num_threads(num_threads);

  Bicop bicop_cpp;
  bicop_cpp.set_var_types(var_types);
  bicop_cpp.select(data, controls);

  return bicop_wrap(bicop_cpp, TRUE);
}

// [[Rcpp::export()]]
Eigen::VectorXd bicop_pdf_cpp(const Eigen::MatrixXd& u,
                              const Rcpp::List& bicop_r)
{
  return bicop_wrap(bicop_r).pdf(u);
}

// [[Rcpp::export()]]
Eigen::VectorXd bicop_cdf_cpp(const Eigen::MatrixXd& u,
                              const Rcpp::List& bicop_r)
{
  return bicop_wrap(bicop_r).cdf(u);
}

// [[Rcpp::export()]]
Eigen::VectorXd bicop_hfunc1_cpp(const Eigen::MatrixXd& u,
                                 const Rcpp::List& bicop_r)
{
  return bicop_wrap(bicop_r).hfunc1(u);
}

// [[Rcpp::export()]]
Eigen::VectorXd bicop_hfunc2_cpp(const Eigen::MatrixXd& u,
                                 const Rcpp::List& bicop_r)
{
  return bicop_wrap(bicop_r).hfunc2(u);
}

// [[Rcpp::export()]]
Eigen::VectorXd bicop_hinv1_cpp(const Eigen::MatrixXd& u,
                                const Rcpp::List& bicop_r)
{
  return bicop_wrap(bicop_r).hinv1(u);
}

// [[Rcpp::export()]]
Eigen::VectorXd bicop_hinv2_cpp(const Eigen::MatrixXd& u,
                                const Rcpp::List& bicop_r)
{
  return bicop_wrap(bicop_r).hinv2(u);
}

// [[Rcpp::export()]]
Eigen::MatrixXd bicop_sim_cpp(const Rcpp::List& bicop_r,
                              const size_t &n,
                              const bool qrng,
                              std::vector<int> seeds)
{
  return bicop_wrap(bicop_r).simulate(n, qrng, seeds);
}

// [[Rcpp::export()]]
double bicop_par_to_tau_cpp(const Rcpp::List& bicop_r)
{
  Bicop bicop_cpp = bicop_wrap(bicop_r);
  return bicop_cpp.parameters_to_tau(bicop_cpp.get_parameters());
}

// [[Rcpp::export()]]
Eigen::MatrixXd bicop_tau_to_par_cpp(const Rcpp::List& bicop_r,
                                     const double& tau)
{
  Bicop bicop_cpp = bicop_wrap(bicop_r);
  return bicop_cpp.tau_to_parameters(tau);
}

// structure exports ---------------------------------------------
// [[Rcpp::export()]]
Rcpp::List rvine_structure_cpp(const Rcpp::List& rvine_structure_r,
                               bool check,
                               bool is_natural_order)
{
  auto rvine_structure = rvine_structure_wrap(rvine_structure_r, check,
                                              is_natural_order);
  return rvine_structure_wrap(rvine_structure);
}

// [[Rcpp::export()]]
void rvine_structure_check_cpp(const Rcpp::List& rvine_struct,
                               bool is_natural_order) {
  auto rvine_structure = rvine_structure_wrap(rvine_struct, true,
                                              is_natural_order);
}

// [[Rcpp::export()]]
Rcpp::List rvine_structure_sim_cpp(size_t d,
                                   bool natural_order,
                                   const std::vector<int>& seeds) {

  auto rvs = RVineStructure::simulate(d, natural_order, seeds);
  return rvine_structure_wrap(rvs);
}


// [[Rcpp::export()]]
void rvine_matrix_check_cpp(
    Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> matrix) {
  auto rvine_structure = RVineStructure(matrix);
}

// vinecop exports --------------------------------------------

// [[Rcpp::export()]]
void vinecop_check_cpp(Rcpp::List vinecop_r)
{
  vinecop_wrap(vinecop_r, true);
}

// [[Rcpp::export()]]
Eigen::MatrixXd vinecop_inverse_rosenblatt_cpp(const Eigen::MatrixXd& U,
                                               const Rcpp::List& vinecop_r,
                                               size_t cores)
{
  return vinecop_wrap(vinecop_r).inverse_rosenblatt(U, cores);
}

// [[Rcpp::export()]]
Eigen::MatrixXd vinecop_rosenblatt_cpp(const Eigen::MatrixXd& U,
                                       const Rcpp::List& vinecop_r,
                                       size_t cores,
                                       bool randomize_discrete,
                                       std::vector<int> seeds)
{
  return vinecop_wrap(vinecop_r).rosenblatt(U, cores, randomize_discrete, seeds);
}

// [[Rcpp::export()]]
Eigen::MatrixXd vinecop_sim_cpp(const Rcpp::List& vinecop_r,
                                const size_t n,
                                const bool qrng,
                                size_t cores,
                                std::vector<int> seeds)
{
  return vinecop_wrap(vinecop_r).simulate(n, qrng, cores, seeds);
}

// [[Rcpp::export()]]
Eigen::VectorXd vinecop_pdf_cpp(const Eigen::MatrixXd& u,
                                const Rcpp::List& vinecop_r,
                                size_t cores)
{
  return vinecop_wrap(vinecop_r).pdf(u, cores);
}

// [[Rcpp::export()]]
Eigen::VectorXd vinecop_cdf_cpp(const Eigen::MatrixXd& u,
                                const Rcpp::List& vinecop_r,
                                size_t N,
                                size_t cores,
                                std::vector<int> seeds)
{
  return vinecop_wrap(vinecop_r).cdf(u, N, cores, seeds);
}

// [[Rcpp::export()]]
Rcpp::List vinecop_select_cpp(const Eigen::MatrixXd& data,
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
                              bool select_families,
                              bool allow_rotations,
                              bool show_trace,
                              size_t num_threads,
                              std::vector<std::string> var_types)
{
  std::vector<BicopFamily> fam_set(family_set.size());
  for (unsigned int fam = 0; fam < fam_set.size(); ++fam) {
    fam_set[fam] = to_cpp_family(family_set[fam]);
  }

  FitControlsVinecop fit_controls;
  fit_controls.set_family_set(fam_set);
  fit_controls.set_parametric_method(par_method);
  fit_controls.set_nonparametric_method(nonpar_method);
  fit_controls.set_nonparametric_mult(mult);
  fit_controls.set_selection_criterion(selection_criterion);
  fit_controls.set_weights(weights);
  fit_controls.set_psi0(psi0);
  fit_controls.set_preselect_families(preselect_families);
  fit_controls.set_allow_rotations(allow_rotations);
  fit_controls.set_num_threads(num_threads);
  fit_controls.set_trunc_lvl(truncation_level);
  fit_controls.set_tree_criterion(tree_criterion);
  fit_controls.set_threshold(threshold);
  fit_controls.set_select_threshold(select_threshold);
  fit_controls.set_select_trunc_lvl(select_truncation_level);
  fit_controls.set_select_families(select_families);
  fit_controls.set_show_trace(show_trace);

  Vinecop vinecop_cpp(rvine_structure_wrap(structure, false));
  vinecop_cpp.set_var_types(var_types);
  vinecop_cpp.select(data, fit_controls);

  return vinecop_wrap(vinecop_cpp, TRUE);
}

// [[Rcpp::export()]]
Rcpp::List vinecop_fit_cpp(const Eigen::MatrixXd& data,
                              Rcpp::List& vinecop_r,
                              std::string par_method,
                              std::string nonpar_method,
                              double mult,
                              const Eigen::VectorXd& weights,
                              bool show_trace,
                              size_t num_threads)
{
  FitControlsVinecop fit_controls;
  fit_controls.set_parametric_method(par_method);
  fit_controls.set_nonparametric_method(nonpar_method);
  fit_controls.set_nonparametric_mult(mult);
  fit_controls.set_weights(weights);
  fit_controls.set_show_trace(show_trace);
  fit_controls.set_num_threads(num_threads);

  Vinecop vinecop_cpp = vinecop_wrap(vinecop_r, false);
  vinecop_cpp.fit(data, fit_controls);

  return vinecop_wrap(vinecop_cpp, true);
}

// [[Rcpp::export()]]
std::vector<Rcpp::List> fit_margins_cpp(const Eigen::MatrixXd& data,
                                        const Eigen::VectorXd& xmin,
                                        const Eigen::VectorXd& xmax,
                                        const std::vector<std::string>& type,
                                        const Eigen::VectorXd& mult,
                                        const Eigen::VectorXd& bw,
                                        const Eigen::VectorXi& deg,
                                        const Eigen::VectorXd& weights,
                                        size_t num_threads)
{
  size_t d = data.cols();
  std::vector<kde1d::Kde1d> fits_cpp(d);
  num_threads = (num_threads > 1) ? num_threads : 0;
  RcppThread::parallelFor(0,
                          d,
                          [&](const size_t& k) {
                            fits_cpp[k] = kde1d::Kde1d(
                              xmin(k),
                              xmax(k),
                              type.at(k),
                              mult(k),
                              bw(k),
                              deg(k)
                            );
                            fits_cpp[k].fit(data.col(k), weights);
                          },
                          num_threads);

  // we can't do the following in parallel because it calls R API
  std::vector<Rcpp::List> fits_r(d);
  for (size_t k = 0; k < d; ++k) {
    fits_r[k] = kde1d::kde1d_wrap(fits_cpp[k]);
  }
  return fits_r;
}




