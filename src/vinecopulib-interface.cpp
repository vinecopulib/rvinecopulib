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
                            std::vector<std::string> var_types)
{
  std::vector<BicopFamily> fam_set(family_set.size());
  for (unsigned int fam = 0; fam < fam_set.size(); ++fam) {
    fam_set[fam] = to_cpp_family(family_set[fam]);
  }
  FitControlsBicop controls(
      fam_set,
      par_method,
      nonpar_method,
      mult,
      selcrit,
      weights,
      psi0,
      presel,
      num_threads
  );
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
double bicop_loglik_cpp(Eigen::MatrixXd& u,
                        const Rcpp::List& bicop_r)
{
  return bicop_wrap(bicop_r).loglik(u);
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
                                       size_t cores)
{
  return vinecop_wrap(vinecop_r).rosenblatt(U, cores);
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
double vinecop_loglik_cpp(const Eigen::MatrixXd& u,
                          const Rcpp::List& vinecop_r,
                          size_t cores)
{
  return vinecop_wrap(vinecop_r).loglik(u, cores);
}

// [[Rcpp::export()]]
double vinecop_mbicv_cpp(const Eigen::MatrixXd& u,
                         const Rcpp::List& vinecop_r,
                         double psi0,
                         size_t cores)
{
  return vinecop_wrap(vinecop_r).mbicv(u, psi0, cores);
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
                              bool show_trace,
                              size_t num_threads,
                              std::vector<std::string> var_types)
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
      weights,
      psi0,
      preselect_families,
      select_truncation_level,
      select_threshold,
      show_trace,
      num_threads
  );

  Vinecop vinecop_cpp(rvine_structure_wrap(structure, false));
  vinecop_cpp.set_var_types(var_types);
  vinecop_cpp.select(data, fit_controls);

  return vinecop_wrap(vinecop_cpp, TRUE);
}

// [[Rcpp::export()]]
std::vector<Rcpp::List> fit_margins_cpp(const Eigen::MatrixXd& data,
                                        const Eigen::VectorXi& nlevels,
                                        const Eigen::VectorXd& mult,
                                        const Eigen::VectorXd& xmin,
                                        const Eigen::VectorXd& xmax,
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
                            fits_cpp[k] = kde1d::Kde1d(data.col(k),
                                                       nlevels(k),
                                                       bw(k),
                                                       mult(k),
                                                       xmin(k),
                                                       xmax(k),
                                                       deg(k),
                                                       weights);
                          },
                          num_threads);

  // we can't do the following in parallel because it calls R API
  std::vector<Rcpp::List> fits_r(d);
  for (size_t k = 0; k < d; ++k) {
    fits_r[k] = kde1d::kde1d_wrap(fits_cpp[k]);
  }
  return fits_r;
}

//' Rosenblatt transform for discrete variables
//'
//' @param u data
//' @param vinecop_r model
//' @param seed seed
//' @param num_threads cores
//'
//' @export
//' @examples
//' a <- 1
// [[Rcpp::export]]
Eigen::MatrixXd rosenblatt_discrete(const Eigen::MatrixXd& u,
                                    const Rcpp::List& vinecop_r,
                                    int seed = 5,
                                    size_t num_threads = 1)
{
  auto d = vinecop_wrap(vinecop_r).get_dim();
  auto R = tools_stats::simulate_uniform(u.rows(), d, true, {seed});
  auto w = vinecop_wrap(vinecop_r).rosenblatt_discrete(u, num_threads);
  return w.leftCols(d).array() * R.array() +
    w.rightCols(d).array() * (1 - R.array());
}



//' asd
 //' @param u data
 //'
 //' @export
 //' @examples
 //' a <- 1
 // [[Rcpp::export]]
 std::vector<size_t> which_in_box(const std::vector<size_t>& ord0,
                                  const std::vector<size_t>& ord1,
                                  const Eigen::MatrixXd& lower,
                                  const Eigen::MatrixXd& upper) {
   auto n = ord1.size();
   std::vector<size_t> indices;
   indices.reserve(n + 1);

   auto l0 = static_cast<size_t>(std::ceil(lower(0) * n));
   auto l1 = static_cast<size_t>(std::ceil(lower(1) * n));
   auto u0 = static_cast<size_t>(std::floor(upper(0) * n));
   auto u1 = static_cast<size_t>(std::floor(upper(1) * n));

   for (size_t k = l0; k < std::min(u0, n); k++) {
     if (l1 <= ord1.at(ord0.at(k)) & ord1.at(ord0.at(k)) <= u1)
       indices.push_back(ord0.at(k));
   }

   return indices;
 }

//' asd
//' @param u data
//'
//' @export
//' @examples
//' a <- 1
// [[Rcpp::export]]
Eigen::MatrixXd find_latent_sample(const Eigen::MatrixXd& u, double b, size_t niter = 5)
{
  size_t n = u.rows();

  auto w = tools_stats::simulate_uniform(n, 2);
  Eigen::MatrixXd uu = w.array() * u.leftCols(2).array() +
    (1 - w.array()) * u.rightCols(2).array();
  // auto x = tools_stats::qnorm(uu);
  // auto xu = tools_stats::safe_qnorm(u.leftCols(2));
  // auto xl = tools_stats::safe_qnorm(u.rightCols(2));
  // xu.array() += b;
  // xl.array() -= b;

  Eigen::MatrixXd x(n, 2), norm_sim(n, 2);

  std::vector<size_t> ord0(n), ord1(n);

  for (size_t it = 0; it < niter; it++) {
    uu = tools_stats::to_pseudo_obs(uu);
    x = tools_stats::qnorm(uu);
    norm_sim = tools_stats::simulate_normal(n, 2).array() * b;
    w = tools_stats::simulate_uniform(n, 1);
    ord0 = tools_eigen::get_order(uu.col(0));
    ord1 = tools_eigen::get_order(uu.col(1));


    for (size_t i = 0; i < n; i++) {
      // indices = which_in_box(ord0, ord1,
      //                             u.row(i).rightCols(2).array(),
      //                             u.row(i).leftCols(2).array());
      std::vector<size_t> indices;
      indices.reserve(n);
      for (size_t k = 0; k < n; k++) {
        if (u(i, 2) <= uu(k, 0) && uu(k, 0) <= u(i, 0) && u(i, 3) <= uu(k, 1) && uu(k, 1) <= u(i, 1))
          indices.push_back(k);
      }
      if (indices.size() > 0) {
        int j = indices.at(static_cast<size_t>(w(i) * indices.size()));
        x.row(i) = x.row(j) + norm_sim.row(i);
      }
    }
  }

  return tools_stats::to_pseudo_obs(x);
}


//' asd
 //' @param u data
 //'
 //' @export
 //' @examples
 //' a <- 1
 // [[Rcpp::export]]
 Eigen::MatrixXd find_latent_sample2(const Eigen::MatrixXd& u, double b, size_t niter = 5)
 {
   size_t n = u.rows();

   auto w = tools_stats::simulate_uniform(n, 2);
   Eigen::MatrixXd uu = w.array() * u.leftCols(2).array() +
     (1 - w.array()) * u.rightCols(2).array();

   auto covering = BoxCovering(uu);
   std::vector<size_t> indices;

   Eigen::MatrixXd x(n, 2), norm_sim(n, 2);

   for (size_t it = 0; it < niter; it++) {
     uu = tools_stats::to_pseudo_obs(uu);
     x = tools_stats::qnorm(uu);
     norm_sim = tools_stats::simulate_normal(n, 2).array() * b;
     w = tools_stats::simulate_uniform(n, 1);

     for (size_t i = 0; i < n; i++) {
       indices = covering.get_box_indices(u.row(i).rightCols(2), u.row(i).leftCols(2));
       if (indices.size() > 0) {
         int j = indices.at(static_cast<size_t>(w(i) * indices.size()));
         x.row(i) = x.row(j) + norm_sim.row(i);
         uu.row(i) = tools_stats::pnorm(x.row(i));
         covering.swap_sample(i, uu.row(i));
       }
     }
   }

   return tools_stats::to_pseudo_obs(x);
 }




