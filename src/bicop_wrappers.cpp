#include <RcppEigen.h>
#include "bicop_wrappers.hpp"
using namespace vinecopulib;


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

Rcpp::List bicop_wrap(Bicop bicop_cpp)
{
    return Rcpp::List::create(
        Rcpp::Named("family")     = to_r_family(bicop_cpp.get_family()),
        Rcpp::Named("rotation")   = bicop_cpp.get_rotation(),
        Rcpp::Named("parameters") = bicop_cpp.get_parameters(),
        Rcpp::Named("npars")      = bicop_cpp.calculate_npars()
    );
}

// [[Rcpp::export]]
void bicop_check_cpp(const Rcpp::List& bicop_r)
{
    bicop_wrap(bicop_r);
}

// [[Rcpp::export()]]
Rcpp::List bicop_select_cpp(
        Eigen::MatrixXd& data,
        std::vector<std::string> family_set,
        std::string par_method,
        std::string nonpar_method,
        double mult,
        std::string selcrit,
        bool presel,
        size_t num_threads
)
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
            presel,
            num_threads
    );
    Bicop bicop_cpp(data, controls);

    return bicop_wrap(bicop_cpp);
}

// [[Rcpp::export()]]
Eigen::VectorXd bicop_pdf_cpp(const Eigen::MatrixXd& u, const Rcpp::List& bicop_r)
{
    return bicop_wrap(bicop_r).pdf(u);
}

// [[Rcpp::export()]]
Eigen::VectorXd bicop_cdf_cpp(const Eigen::MatrixXd& u, const Rcpp::List& bicop_r)
{
    return bicop_wrap(bicop_r).cdf(u);
}

// [[Rcpp::export()]]
Eigen::VectorXd bicop_hfunc1_cpp(const Eigen::MatrixXd& u, const Rcpp::List& bicop_r)
{
    return bicop_wrap(bicop_r).hfunc1(u);
}

// [[Rcpp::export()]]
Eigen::VectorXd bicop_hfunc2_cpp(const Eigen::MatrixXd& u, const Rcpp::List& bicop_r)
{
    return bicop_wrap(bicop_r).hfunc2(u);
}

// [[Rcpp::export()]]
Eigen::VectorXd bicop_hinv1_cpp(const Eigen::MatrixXd& u, const Rcpp::List& bicop_r)
{
    return bicop_wrap(bicop_r).hinv1(u);
}

// [[Rcpp::export()]]
Eigen::VectorXd bicop_hinv2_cpp(const Eigen::MatrixXd& u, const Rcpp::List& bicop_r)
{
    return bicop_wrap(bicop_r).hinv2(u);
}

// [[Rcpp::export()]]
double bicop_loglik_cpp(Eigen::MatrixXd& u, const Rcpp::List& bicop_r)
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
Eigen::MatrixXd bicop_tau_to_par_cpp(const Rcpp::List& bicop_r, const double& tau)
{
    Bicop bicop_cpp = bicop_wrap(bicop_r);
    return bicop_cpp.tau_to_parameters(tau);
}
