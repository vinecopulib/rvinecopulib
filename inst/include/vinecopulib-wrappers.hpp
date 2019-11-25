#pragma once

#ifndef BOOST_NO_AUTO_PTR
#define BOOST_NO_AUTO_PTR
#endif

#ifndef BOOST_MATH_PROMOTE_DOUBLE_POLICY
#define BOOST_MATH_PROMOTE_DOUBLE_POLICY false
#else
#undef BOOST_MATH_PROMOTE_DOUBLE_POLICY
#define BOOST_MATH_PROMOTE_DOUBLE_POLICY false
#endif


#include <RcppEigen.h>
#include "vinecopulib.hpp"

namespace vinecopulib {

// bicop wrapppers -----------------------------------

inline BicopFamily to_cpp_family(const std::string& fam)
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

inline std::string to_r_family(const BicopFamily& fam)
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

inline Bicop bicop_wrap(const Rcpp::List& bicop_r)
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
  bicop_cpp.set_var_types(bicop_r["var_types"]);

  return bicop_cpp;
}

inline Rcpp::List bicop_wrap(Bicop bicop_cpp, bool is_fitted)
{
  double loglik = NAN;
  if (is_fitted)
    loglik = bicop_cpp.get_loglik();
  auto bc = Rcpp::List::create(
    Rcpp::Named("family")     = to_r_family(bicop_cpp.get_family()),
    Rcpp::Named("rotation")   = bicop_cpp.get_rotation(),
    Rcpp::Named("parameters") = bicop_cpp.get_parameters(),
    Rcpp::Named("var_types")  = bicop_cpp.get_var_types(),
    Rcpp::Named("npars")      = bicop_cpp.get_npars(),
    Rcpp::Named("loglik")     = loglik
  );
  bc.attr("class") = "bicop_dist";
  return bc;
}

// structure wrappers ---------------------------------------

inline TriangularArray<size_t> struct_array_wrap(const Rcpp::List& struct_array_r,
                                                 size_t trunc_lvl)
{
  std::vector<std::vector<size_t>> rows(trunc_lvl);
  for (size_t i = 0; i < trunc_lvl; i++) {
    rows.at(i) = Rcpp::as<std::vector<size_t>>(struct_array_r[i]);
  }
  return TriangularArray<size_t>(rows);
}

inline Rcpp::List struct_array_wrap(const TriangularArray<size_t>& struct_array)
{
  size_t d = struct_array.get_dim();
  size_t trunc_lvl = struct_array.get_trunc_lvl();
  Rcpp::List struct_array_r(trunc_lvl);
  for (size_t i = 0; i < trunc_lvl; i++) {
    std::vector<size_t> newrow(d - 1 - i);
    for (size_t j = 0; j < d - 1 - i; j++) {
      newrow[j] = struct_array(i, j);
    }
    struct_array_r[i] = newrow;
  }

  return struct_array_r;
}

inline RVineStructure rvine_structure_wrap(const Rcpp::List& rvine_structure_r,
                                           bool check = false,
                                           bool is_natural_order = true)
{
  size_t trunc_lvl = rvine_structure_r["trunc_lvl"];
  std::vector<size_t> order = rvine_structure_r["order"];
  TriangularArray<size_t> struct_array = struct_array_wrap(
    rvine_structure_r["struct_array"], trunc_lvl);

  return RVineStructure(order, struct_array, is_natural_order, check);
}

inline Rcpp::List rvine_structure_wrap(const RVineStructure& rvine_struct)
{
  auto struct_array = struct_array_wrap(rvine_struct.get_struct_array(true));
  auto rvs = Rcpp::List::create(
    Rcpp::Named("order")        = rvine_struct.get_order(),
    Rcpp::Named("struct_array") = struct_array,
    Rcpp::Named("d")            = rvine_struct.get_dim(),
    Rcpp::Named("trunc_lvl")    = rvine_struct.get_trunc_lvl()
  );
  rvs.attr("class") = Rcpp::CharacterVector{"rvine_structure", "list"};
  return rvs;
}


// vinecop wrappers --------------------------------------------------------

inline std::vector<std::vector<Bicop>> pair_copulas_wrap(
    const Rcpp::List& pair_copulas_r,
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

inline Rcpp::List pair_copulas_wrap(
    std::vector<std::vector<Bicop>> pair_copulas,
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

inline Vinecop vinecop_wrap(const Rcpp::List& vinecop_r, bool check = false)
{
  // omit R-vine matrix check, already done in R
  auto structure = rvine_structure_wrap(vinecop_r["structure"], check);
  auto pair_copulas = pair_copulas_wrap(vinecop_r["pair_copulas"],
                                        structure.get_dim());
  Vinecop vc(structure, pair_copulas);
  vc.set_var_types(vinecop_r["var_types"]);
  return vc;
}

inline Rcpp::List vinecop_wrap(const Vinecop& vinecop_cpp,
                               bool is_fitted = false)
{

  auto vine_structure = rvine_structure_wrap(
    vinecop_cpp.get_rvine_structure());
  auto pair_copulas = pair_copulas_wrap(vinecop_cpp.get_all_pair_copulas(),
                                        vinecop_cpp.get_dim(), is_fitted);

  double npars = vinecop_cpp.get_npars();
  double threshold = vinecop_cpp.get_threshold();
  double loglik = NAN;
  if (is_fitted)
    loglik = vinecop_cpp.get_loglik();
  auto vc = Rcpp::List::create(
    Rcpp::Named("pair_copulas")      = pair_copulas,
    Rcpp::Named("structure")         = vine_structure,
    Rcpp::Named("var_types")         = vinecop_cpp.get_var_types(),
    Rcpp::Named("npars")             = npars,
    Rcpp::Named("loglik")            = loglik,
    Rcpp::Named("threshold")         = threshold
  );
  vc.attr("class") = Rcpp::CharacterVector{"vinecop", "vinecop_dist"};
  return vc;
}

}
