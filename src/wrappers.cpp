#include "wrappers.hpp"

// tools exports

// [[Rcpp::export]]
Eigen::MatrixXd pseudo_obs_cpp(Eigen::MatrixXd x, std::string ties_method) 
{
    return vinecopulib::tools_stats::to_pseudo_obs(x, ties_method);
}

// bicop wrappers

BicopFamily to_cpp_family(const std::string& fam)
{
    BicopFamily bicop_fam;
    if (fam == "indep") {
        bicop_fam = BicopFamily::indep;
    } else if (fam == "gaussian") {
        bicop_fam = BicopFamily::gaussian;
    } else if (fam == "t") {
        bicop_fam = BicopFamily::student;
    } else if (fam == "clayton") {
        bicop_fam = BicopFamily::clayton;
    } else if (fam == "gumbel") {
        bicop_fam = BicopFamily::gumbel;
    } else if (fam == "frank") {
        bicop_fam = BicopFamily::frank;
    } else if (fam == "joe") {
        bicop_fam = BicopFamily::joe;
    } else if (fam == "bb1") {
        bicop_fam = BicopFamily::bb1;
    } else if (fam == "bb6") {
        bicop_fam = BicopFamily::bb6;
    } else if (fam == "bb7") {
        bicop_fam = BicopFamily::bb7;
    } else if (fam == "bb8") {
        bicop_fam = BicopFamily::bb8;
    } else if (fam == "tll") {
        bicop_fam = BicopFamily::tll;
    } else {
        throw std::runtime_error("family not implemented");
    }

    return bicop_fam;
}

std::string to_r_family(const BicopFamily& fam)
{
    std::string bicop_fam;
    if (fam == BicopFamily::indep) {
        bicop_fam = "indep";
    } else if (fam == BicopFamily::gaussian) {
        bicop_fam = "gaussian";
    } else if (fam == BicopFamily::student) {
        bicop_fam = "t";
    } else if (fam == BicopFamily::clayton) {
        bicop_fam = "clayton";
    } else if (fam == BicopFamily::gumbel) {
        bicop_fam = "gumbel";
    } else if (fam == BicopFamily::frank) {
        bicop_fam = "frank";
    } else if (fam == BicopFamily::joe) {
        bicop_fam = "joe";
    } else if (fam == BicopFamily::bb1) {
        bicop_fam = "bb1";
    } else if (fam == BicopFamily::bb6) {
        bicop_fam = "bb6";
    } else if (fam == BicopFamily::bb7) {
        bicop_fam = "bb7";
    } else if (fam == BicopFamily::bb8) {
        bicop_fam = "bb8";
    } else if (fam == BicopFamily::tll) {
        bicop_fam = "tll";
    } else {
        throw std::runtime_error("family not implemented");
    }

    return bicop_fam;
}

Bicop bicop_wrap(const Rcpp::List& bicop_r)
{
    Eigen::MatrixXd par = bicop_r["parameters"];
    Bicop bicop_cpp;
    if (par.size() == 0) {
        bicop_cpp = Bicop(
            to_cpp_family(bicop_r["family"]),
            bicop_r["rotation"]
        );
    } else {
        Eigen::MatrixXd pars = bicop_r["parameters"];
        bicop_cpp = Bicop(
            to_cpp_family(bicop_r["family"]),
            bicop_r["rotation"],
            pars
        );
    }

    return bicop_cpp;
}

Rcpp::List bicop_wrap(Bicop bicop_cpp, bool is_fitted)
{
    double loglik = NAN;
    if (is_fitted)
        loglik = bicop_cpp.get_loglik();
    return Rcpp::List::create(
        Rcpp::Named("family")     = to_r_family(bicop_cpp.get_family()),
        Rcpp::Named("rotation")   = bicop_cpp.get_rotation(),
        Rcpp::Named("parameters") = bicop_cpp.get_parameters(),
        Rcpp::Named("npars")      = bicop_cpp.calculate_npars(),
        Rcpp::Named("loglik")     = loglik
    );
}

// bicop exports

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

// structure wrappers

TriangularArray<size_t> struct_array_wrap(const Rcpp::List& struct_array_r,
                                          size_t trunc_lvl)
{
    size_t d = struct_array_r.size() + 1;
    auto struct_array = TriangularArray<size_t>(d, trunc_lvl);
    for (size_t i = 0; i < d - 1; i++) {
        struct_array.set_column(i, struct_array_r[i]);
    }
    
    return struct_array;
}

Rcpp::List struct_array_wrap(const TriangularArray<size_t>& struct_array)
{
    size_t d = struct_array.get_dim();
    
    Rcpp::List struct_array_r(d - 1);
    for (size_t i = 0; i < d - 1; i++) {
        struct_array_r[i] = struct_array[i];
    }
    
    return struct_array_r;
}

RVineStructure rvine_structure_wrap(const Rcpp::List& rvine_structure_r,
                                    bool check,
                                    bool is_natural_order)
{
    size_t trunc_lvl = rvine_structure_r["trunc_lvl"];
    std::vector<size_t> order = rvine_structure_r["order"];
    TriangularArray<size_t> struct_array = struct_array_wrap(
        rvine_structure_r["struct_array"], trunc_lvl);
    
    return RVineStructure(order, struct_array, is_natural_order, check);
}

Rcpp::List rvine_structure_wrap(const RVineStructure& rvine_struct)
{
    auto struct_array = struct_array_wrap(rvine_struct.get_struct_array());
    return Rcpp::List::create(
        Rcpp::Named("order")        = rvine_struct.get_order(),
        Rcpp::Named("struct_array") = struct_array,
        Rcpp::Named("d")            = rvine_struct.get_dim(),
        Rcpp::Named("trunc_lvl")    = rvine_struct.get_trunc_lvl()
    );
}


// structure exports
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

// vinecop wrappers

std::vector<std::vector<Bicop>> pair_copulas_wrap(const Rcpp::List& pair_copulas_r,
                                                  size_t d)
{
    size_t trunc_lvl = pair_copulas_r.size();
    auto pair_copulas = Vinecop::make_pair_copula_store(d, trunc_lvl);
    Rcpp::List tree_pcs;
    for (size_t t = 0; t < trunc_lvl; ++t) {
        tree_pcs = pair_copulas_r[t];
        if (static_cast<size_t>(tree_pcs.size()) != d - 1 - t) {
            throw std::runtime_error("length(pair_copulas[[t]]) must be d-t");
        }
        for(size_t e = 0; e < d - 1 - t; ++e) {
            pair_copulas[t][e] = bicop_wrap(tree_pcs[e]);
        }
    }
    return pair_copulas;
}

Rcpp::List pair_copulas_wrap(std::vector<std::vector<Bicop>> pair_copulas,
                             size_t d,
                             bool is_fitted)
{
    Rcpp::List pair_copulas_r(pair_copulas.size());
    for (size_t t = 0; t < pair_copulas.size(); ++t) {
        Rcpp::List tree_pcs(d - 1 - t);
        for(size_t e = 0; e < d - 1 - t; ++e) {
            tree_pcs[e] = bicop_wrap(pair_copulas[t][e], is_fitted);
        }
        pair_copulas_r[t] = tree_pcs;
    }
    
    return pair_copulas_r;
}

Vinecop vinecop_wrap(const Rcpp::List& vinecop_r, bool check)
{
    // omit R-vine matrix check, already done in R
    auto structure = rvine_structure_wrap(vinecop_r["structure"], check);
    
    // extract pair-copulas
    auto pair_copulas = pair_copulas_wrap(vinecop_r["pair_copulas"],
                                          structure.get_dim());

    return Vinecop(pair_copulas, structure);
}

Rcpp::List vinecop_wrap(const Vinecop& vinecop_cpp, 
                        bool is_fitted) 
{
    
    auto vine_structure = rvine_structure_wrap(vinecop_cpp.get_rvine_structure());
    auto pair_copulas = pair_copulas_wrap(vinecop_cpp.get_all_pair_copulas(), 
                                          vinecop_cpp.get_dim(), is_fitted);
    
    double npars = vinecop_cpp.calculate_npars();
    double threshold = vinecop_cpp.get_threshold();
    double loglik = NAN;
    if (is_fitted)
        loglik = vinecop_cpp.get_loglik();
    return Rcpp::List::create(
        Rcpp::Named("pair_copulas")      = pair_copulas,
        Rcpp::Named("structure")         = vine_structure,
        Rcpp::Named("npars")             = npars,
        Rcpp::Named("loglik")            = loglik,
        Rcpp::Named("threshold")         = threshold
    );
}

// vinecop exports

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
