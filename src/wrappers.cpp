#include "wrappers.hpp"

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
                            size_t num_threads)
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
    Bicop bicop_cpp(data, controls);

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
    return(rvine_structure_wrap(rvine_structure));
}

// [[Rcpp::export()]]
void rvine_structure_check_cpp(const Rcpp::List& rvine_struct,
                               bool is_natural_order) {

    auto rvine_structure = rvine_structure_wrap(rvine_struct, true,
                                                is_natural_order);
}

// [[Rcpp::export()]]
void rvine_matrix_check_cpp(Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> matrix) {
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
                              size_t num_threads)
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
    Vinecop vinecop_cpp(data.cols());
    if (is_structure_provided) {
        vinecop_cpp = Vinecop(rvine_structure_wrap(structure, false));
        vinecop_cpp.select_families(data, fit_controls);
    } else {
        vinecop_cpp.select_all(data, fit_controls);
    }
    return vinecop_wrap(vinecop_cpp, TRUE);
}
