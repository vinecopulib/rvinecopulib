// Copyright © 2018 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_stl.hpp>
#include <stdexcept>
#include <thread>

//! Tools for bivariate and vine copula modeling
namespace vinecopulib {
//! creates the controls for fitting bivariate copula models.
//! @param family_set the set of copula families to consider (if empty, then
//!     all families are included).
//! @param parametric_method the fit method for parametric families;
//!     possible choices: `"mle"`, `"itau"`.
//! @param nonparametric_method the fit method for the local-likelihood
//!     nonparametric family (TLLs); possible choices: `"constant"`,
//!     `"linear"`, `"quadratic"`.
//! @param nonparametric_mult a factor with which the smoothing parameters
//!     are multiplied.
//! @param selection_criterion the selection criterion (`"loglik"`, `"aic"` 
//!     or `"bic"`).
//! @param preselect_families whether to exclude families before fitting
//!     based on symmetry properties of the data.
//! @param num_threads number of concurrent threads to use while fitting
//!     copulas for different families; never uses more than the number
//!     returned by `std::thread::hardware_concurrency()`.
inline FitControlsBicop::FitControlsBicop(std::vector <BicopFamily> family_set,
                                          std::string parametric_method,
                                          std::string nonparametric_method,
                                          double nonparametric_mult,
                                          std::string selection_criterion,
                                          bool preselect_families,
                                          size_t num_threads)
{
    set_family_set(family_set);
    set_parametric_method(parametric_method);
    set_nonparametric_method(nonparametric_method);
    set_nonparametric_mult(nonparametric_mult);
    set_selection_criterion(selection_criterion);
    set_preselect_families(preselect_families);
    set_num_threads(num_threads);
}

//! @param parametric_method the fit method for parametric families;
//!     possible choices: `"mle"`, `"itau"`.
inline FitControlsBicop::FitControlsBicop(std::string parametric_method) :
    FitControlsBicop()
{
    set_parametric_method(parametric_method);
}

//! @param nonparametric_method the fit method for the local-likelihood
//!     nonparametric family (TLLs); possible choices: `"constant"`,
//!     `"linear"`, `"quadratic"`.
//! @param nonparametric_mult a factor with which the smoothing parameters
//!     are multiplied.
inline FitControlsBicop::FitControlsBicop(std::string nonparametric_method,
                                          double nonparametric_mult) :
    FitControlsBicop()
{
    set_nonparametric_method(nonparametric_method);
    set_nonparametric_mult(nonparametric_mult);
}

//! Sanity checks
//! @{
inline void
FitControlsBicop::check_parametric_method(std::string parametric_method)
{
    if (!tools_stl::is_member(parametric_method, {"itau", "mle"})) {
        throw std::runtime_error("parametric_method should be mle or itau");
    }
}

inline void
FitControlsBicop::check_nonparametric_method(std::string nonparametric_method)
{
    if (!tools_stl::is_member(nonparametric_method,
                              {"constant", "linear", "quadratic"})) {
        throw std::runtime_error(
            "parametric_method should be constant, linear or quadratic");
    }
}

inline void
FitControlsBicop::check_nonparametric_mult(double nonparametric_mult)
{
    if (nonparametric_mult <= 0.0) {
        throw std::runtime_error("nonparametric_mult must be positive");
    }
}

inline void
FitControlsBicop::check_selection_criterion(std::string selection_criterion)
{
    if (!tools_stl::is_member(selection_criterion, {"loglik", "aic", "bic"})) {
        throw std::runtime_error(
            "selection_criterion should be 'loglik', 'aic', or 'bic'");
    }
}
//! @}

//! Getters and setters.
//! @{
inline std::vector <BicopFamily> FitControlsBicop::get_family_set() const
{
    return family_set_;
}

inline std::string FitControlsBicop::get_parametric_method() const
{
    return parametric_method_;
}

inline std::string FitControlsBicop::get_nonparametric_method() const
{
    return nonparametric_method_;
}

inline double FitControlsBicop::get_nonparametric_mult() const
{
    return nonparametric_mult_;
}

inline size_t FitControlsBicop::get_num_threads()
{
    return num_threads_;
}

inline std::string FitControlsBicop::get_selection_criterion() const
{
    return selection_criterion_;
}

inline bool FitControlsBicop::get_preselect_families() const
{
    return preselect_families_;
}

inline void
FitControlsBicop::set_family_set(std::vector <BicopFamily> family_set)
{
    family_set_ = family_set;
}

inline void
FitControlsBicop::set_parametric_method(std::string parametric_method)
{
    check_parametric_method(parametric_method);
    parametric_method_ = parametric_method;
}

inline void
FitControlsBicop::set_nonparametric_method(std::string nonparametric_method)
{
    check_nonparametric_method(nonparametric_method);
    nonparametric_method_ = nonparametric_method;
}

inline void
FitControlsBicop::set_nonparametric_mult(double nonparametric_mult)
{
    check_nonparametric_mult(nonparametric_mult);
    nonparametric_mult_ = nonparametric_mult;
}

inline void
FitControlsBicop::set_selection_criterion(std::string selection_criterion)
{
    check_selection_criterion(selection_criterion);
    selection_criterion_ = selection_criterion;
}

inline void FitControlsBicop::set_preselect_families(bool preselect_families)
{
    preselect_families_ = preselect_families;
}

inline void FitControlsBicop::set_num_threads(size_t num_threads)
{
    num_threads_ = process_num_threads(num_threads);
}

inline size_t FitControlsBicop::process_num_threads(size_t num_threads)
{
    // use at least one thread
    num_threads = std::max(num_threads, static_cast<size_t>(1));
    // don't use more threads than supported by the system
    size_t max_threads = std::thread::hardware_concurrency();
    num_threads = std::min(num_threads, max_threads);

    return num_threads;
}
//! @}
}
