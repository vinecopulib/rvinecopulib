#include "vinecopulib-wrappers.hpp"
#include "kde1d-wrappers.hpp"

using namespace vinecopulib;

inline Eigen::MatrixXd get_bicop_parameters(const Rcpp::List& bicop_r)
{
  return bicop_r["parameters"];
}

inline bool has_vectorized_bicop_parameters(const Bicop& bicop_cpp,
                                            const Eigen::MatrixXd& parameters)
{
  if (!tools_stl::is_member(bicop_cpp.get_family(), bicop_families::parametric) ||
      (parameters.size() == 0)) {
    return false;
  }

  const Eigen::Index p = bicop_cpp.get_parameters().rows();
  if (parameters.cols() == 1) {
    // n x 1 is vectorized only for one-parameter families.
    return (p == 1) && (parameters.rows() > 1);
  }

  return (parameters.cols() == p) && (parameters.rows() > 1);
}


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
  Bicop bicop_cpp = bicop_wrap(bicop_r);
  Eigen::MatrixXd parameters = get_bicop_parameters(bicop_r);
  if (has_vectorized_bicop_parameters(bicop_cpp, parameters)) {
    return bicop_cpp.pdf(u, parameters);
  }
  return bicop_cpp.pdf(u);
}

// [[Rcpp::export()]]
Eigen::VectorXd bicop_cdf_cpp(const Eigen::MatrixXd& u,
                              const Rcpp::List& bicop_r)
{
  Bicop bicop_cpp = bicop_wrap(bicop_r);
  Eigen::MatrixXd parameters = get_bicop_parameters(bicop_r);
  if (has_vectorized_bicop_parameters(bicop_cpp, parameters)) {
    return bicop_cpp.cdf(u, parameters);
  }
  return bicop_cpp.cdf(u);
}

// [[Rcpp::export()]]
Eigen::VectorXd bicop_hfunc1_cpp(const Eigen::MatrixXd& u,
                                 const Rcpp::List& bicop_r)
{
  Bicop bicop_cpp = bicop_wrap(bicop_r);
  Eigen::MatrixXd parameters = get_bicop_parameters(bicop_r);
  if (has_vectorized_bicop_parameters(bicop_cpp, parameters)) {
    return bicop_cpp.hfunc1(u, parameters);
  }
  return bicop_cpp.hfunc1(u);
}

// [[Rcpp::export()]]
Eigen::VectorXd bicop_hfunc2_cpp(const Eigen::MatrixXd& u,
                                 const Rcpp::List& bicop_r)
{
  Bicop bicop_cpp = bicop_wrap(bicop_r);
  Eigen::MatrixXd parameters = get_bicop_parameters(bicop_r);
  if (has_vectorized_bicop_parameters(bicop_cpp, parameters)) {
    return bicop_cpp.hfunc2(u, parameters);
  }
  return bicop_cpp.hfunc2(u);
}

// [[Rcpp::export()]]
Eigen::VectorXd bicop_hinv1_cpp(const Eigen::MatrixXd& u,
                                const Rcpp::List& bicop_r)
{
  Bicop bicop_cpp = bicop_wrap(bicop_r);
  Eigen::MatrixXd parameters = get_bicop_parameters(bicop_r);
  if (has_vectorized_bicop_parameters(bicop_cpp, parameters)) {
    return bicop_cpp.hinv1(u, parameters);
  }
  return bicop_cpp.hinv1(u);
}

// [[Rcpp::export()]]
Eigen::VectorXd bicop_hinv2_cpp(const Eigen::MatrixXd& u,
                                const Rcpp::List& bicop_r)
{
  Bicop bicop_cpp = bicop_wrap(bicop_r);
  Eigen::MatrixXd parameters = get_bicop_parameters(bicop_r);
  if (has_vectorized_bicop_parameters(bicop_cpp, parameters)) {
    return bicop_cpp.hinv2(u, parameters);
  }
  return bicop_cpp.hinv2(u);
}

// [[Rcpp::export()]]
Eigen::MatrixXd bicop_sim_cpp(const Rcpp::List& bicop_r,
                              const size_t &n,
                              const bool qrng,
                              std::vector<int> seeds)
{
  Bicop bicop_cpp = bicop_wrap(bicop_r);
  return bicop_cpp.simulate(n, qrng, seeds);
}

// [[Rcpp::export()]]
double bicop_par_to_tau_cpp(const Rcpp::List& bicop_r)
{
  Bicop bicop_cpp = bicop_wrap(bicop_r);
  return bicop_cpp.parameters_to_tau(bicop_cpp.get_parameters());
}

// [[Rcpp::export()]]
Eigen::MatrixXd bicop_tail_dep_cpp(const Rcpp::List& bicop_r)
{
  return bicop_wrap(bicop_r).get_taildep();
}

// [[Rcpp::export()]]
double bicop_beta_cpp(const Rcpp::List& bicop_r)
{
  return bicop_wrap(bicop_r).get_beta();
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

inline Rcpp::NumericVector eigen_vector_wrap(const Eigen::VectorXd& x)
{
  Rcpp::NumericVector result(x.size());
  std::copy(x.data(), x.data() + x.size(), result.begin());
  return result;
}

inline Rcpp::List vector_triangular_array_wrap(
    const TriangularArray<Eigen::VectorXd>& x)
{
  size_t trunc_lvl = x.get_trunc_lvl();
  size_t d = x.get_dim();
  Rcpp::List result(trunc_lvl);
  for (size_t t = 0; t < trunc_lvl; ++t) {
    Rcpp::List row(d - t - 1);
    for (size_t e = 0; e < d - t - 1; ++e) {
      row[e] = eigen_vector_wrap(x(t, e));
    }
    result[t] = row;
  }
  return result;
}

// [[Rcpp::export()]]
Rcpp::List vinecop_pdf_full_cpp(const Eigen::MatrixXd& u,
                                const Rcpp::List& vinecop_r,
                                size_t cores)
{
  auto result = vinecop_wrap(vinecop_r).pdf_full(u, cores, true);
  return Rcpp::List::create(
    Rcpp::Named("pdf") = eigen_vector_wrap(result.pdf),
    Rcpp::Named("pdf_edges") = vector_triangular_array_wrap(result.pdf_edges),
    Rcpp::Named("hfunc1") = vector_triangular_array_wrap(result.hfunc1),
    Rcpp::Named("hfunc2") = vector_triangular_array_wrap(result.hfunc2),
    Rcpp::Named("hfunc1_sub") = vector_triangular_array_wrap(result.hfunc1_sub),
    Rcpp::Named("hfunc2_sub") = vector_triangular_array_wrap(result.hfunc2_sub)
  );
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
Eigen::MatrixXd vinecop_scores_cpp(const Eigen::MatrixXd& u,
                                   const Rcpp::List& vinecop_r,
                                   bool step_wise,
                                   size_t cores)
{
  return vinecop_wrap(vinecop_r).scores(u, step_wise, cores);
}

// [[Rcpp::export()]]
Eigen::MatrixXd vinecop_hessian_cpp(const Eigen::MatrixXd& u,
                                    const Rcpp::List& vinecop_r,
                                    bool step_wise,
                                    size_t cores)
{
  return vinecop_wrap(vinecop_r).hessian(u, step_wise, cores);
}

// [[Rcpp::export()]]
Rcpp::List vinecop_select_cpp(const Eigen::MatrixXd &data,
                              Rcpp::List &structure,
                              std::vector<std::string> family_set,
                              std::string par_method,
                              std::string nonpar_method,
                              double mult,
                              int truncation_level,
                              std::string tree_criterion,
                              double threshold,
                              std::string selection_criterion,
                              const Eigen::VectorXd &weights,
                              double psi0,
                              bool select_truncation_level,
                              bool select_threshold,
                              bool preselect_families,
                              bool select_families,
                              bool allow_rotations,
                              bool show_trace,
                              size_t num_threads,
                              std::vector<std::string> var_types,
                              std::string tree_algorithm,
                              std::vector<int> seeds)
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
  fit_controls.set_tree_algorithm(tree_algorithm);
  fit_controls.set_seeds(seeds);

  Vinecop vinecop_cpp(rvine_structure_wrap(structure, false));
  vinecop_cpp.set_var_types(var_types);
  vinecop_cpp.select(data, fit_controls);

  return vinecop_wrap(vinecop_cpp, TRUE);
}

// [[Rcpp::export()]]
Rcpp::List vinecop_fit_cpp(const Eigen::MatrixXd &data,
                           Rcpp::List &vinecop_r,
                           std::string par_method,
                           std::string nonpar_method,
                           double mult,
                           const Eigen::VectorXd &weights,
                           bool show_trace,
                           size_t num_threads,
                           std::string tree_algorithm,
                           std::vector<int> seeds)
{
  FitControlsVinecop fit_controls;
  fit_controls.set_parametric_method(par_method);
  fit_controls.set_nonparametric_method(nonpar_method);
  fit_controls.set_nonparametric_mult(mult);
  fit_controls.set_weights(weights);
  fit_controls.set_show_trace(show_trace);
  fit_controls.set_num_threads(num_threads);
  fit_controls.set_tree_algorithm(tree_algorithm);
  fit_controls.set_seeds(seeds);

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

