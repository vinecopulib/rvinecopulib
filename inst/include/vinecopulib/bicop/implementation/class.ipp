// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/bicop/tools_select.hpp>
#include <vinecopulib/misc/tools_stats.hpp>
#include <vinecopulib/misc/tools_stl.hpp>
#include <vinecopulib/misc/tools_interface.hpp>
#include <vinecopulib/misc/tools_parallel.hpp>
#include <mutex>

//! Tools for bivariate and vine copula modeling
namespace vinecopulib {

//! creates the independence copula.
inline Bicop::Bicop()
{
    bicop_ = AbstractBicop::create();
    rotation_ = 0;
}

//! creates a specific bivariate copula model.
//! @param family the copula family.
//! @param rotation the rotation of the copula; one of 0, 90, 180, or 270
//!     (for Independence, Gaussian, Student, Frank, and nonparametric
//!     families, only 0 is allowed).
//! @param parameters the copula parameters.
inline Bicop::Bicop(BicopFamily family, int rotation,
                    const Eigen::MatrixXd &parameters)
{
    bicop_ = AbstractBicop::create(family, parameters);
    // family must be set before checking the rotation
    set_rotation(rotation);
}

//! create a copula model from the data,
//! equivalent to `Bicop cop; cop.select(data, controls)`.
//! @param data see select().
//! @param controls see select().
inline Bicop::Bicop(Eigen::Matrix<double, Eigen::Dynamic, 2> data,
                    FitControlsBicop controls)
{
    select(data, controls);
}

//! creates from a boost::property_tree::ptree object
//! @param input the boost::property_tree::ptree object to convert from
//! (see to_ptree() for the structure of the input).
inline Bicop::Bicop(boost::property_tree::ptree input) :
    Bicop(
        get_family_enum(input.get<std::string>("family")),
        input.get<int>("rotation"),
        tools_serialization::ptree_to_matrix<double>(
            input.get_child("parameters"))
    )
{
}

//! creates from a JSON file
//! @param filename the name of the JSON file to read (see to_ptree() for the
//! structure of the file).
inline Bicop::Bicop(const char *filename) :
    Bicop(tools_serialization::json_to_ptree(filename))
{
}

//! Convert the copula into a boost::property_tree::ptree object
//!
//! The boost::property_tree::ptree is contains of three values named
//! `"family"`, `"rotation"`, `"parameters"`, respectively a string
//! for the family name, an integer for the rotation, and an Eigen::MatrixXd
//! for the parameters.
//!
//! @return the boost::property_tree::ptree object containing the copula.
inline boost::property_tree::ptree Bicop::to_ptree() const
{
    boost::property_tree::ptree output;

    output.put("family", get_family_name());
    output.put("rotation", rotation_);
    auto mat_node = tools_serialization::matrix_to_ptree(get_parameters());
    output.add_child("parameters", mat_node);

    return output;
}

//! Write the copula object into a JSON file
//!
//! See to_ptree() for the structure of the file.
//!
//! @param filename the name of the file to write.
inline void Bicop::to_json(const char *filename) const
{
    boost::property_tree::write_json(filename, to_ptree());
}

//! evaluates the copula density.
//!
//! @param u \f$n \times 2\f$ matrix of evaluation points.
//! @return The copula density evaluated at \c u.
inline Eigen::VectorXd
Bicop::pdf(const Eigen::Matrix<double, Eigen::Dynamic, 2> &u)
const
{
    Eigen::VectorXd f = bicop_->pdf(cut_and_rotate(u));
    f = f.unaryExpr([](const double x) { return std::min(x, 1e16); });
    return f;
}

//! evaluates the copula distribution.
//!
//! @param u \f$n \times 2\f$ matrix of evaluation points.
//! @return The copula distribution evaluated at \c u.
inline Eigen::VectorXd
Bicop::cdf(const Eigen::Matrix<double, Eigen::Dynamic, 2> &u)
const
{
    Eigen::VectorXd p = bicop_->cdf(cut_and_rotate(u));
    switch (rotation_) {
        case 0:
            return p;

        case 90:
            return u.col(1) - p;

        case 180: {
            Eigen::VectorXd f = Eigen::VectorXd::Ones(p.rows());
            f = f - u.rowwise().sum();
            return p - f;
        }

        case 270:
            return u.col(0) - p;

        default:
            throw std::runtime_error(std::string(
                "rotation can only take values in {0, 90, 180, 270}"
            ));
    }
}

//! calculates the first h-function, i.e.,
//! \f$ h_1(u_1, u_2) = \int_0^{u_2} c(u_1, s) \f$.
//! @param u \f$m \times 2\f$ matrix of evaluation points.
inline Eigen::VectorXd
Bicop::hfunc1(const Eigen::Matrix<double, Eigen::Dynamic, 2> &u)
const
{
    switch (rotation_) {
        case 0:
            return bicop_->hfunc1(cut_and_rotate(u));

        case 90:
            return bicop_->hfunc2(cut_and_rotate(u));

        case 180:
            return 1.0 - bicop_->hfunc1(cut_and_rotate(u)).array();

        case 270:
            return 1.0 - bicop_->hfunc2(cut_and_rotate(u)).array();

        default:
            throw std::runtime_error(std::string(
                "rotation can only take values in {0, 90, 180, 270}"
            ));
    }
}

//! calculates the second h-function, i.e.,
//! \f$ h_2(u_1, u_2) = \int_0^{u_1} c(s, u_2) \f$.
//! @param u \f$m \times 2\f$ matrix of evaluation points.
inline Eigen::VectorXd
Bicop::hfunc2(const Eigen::Matrix<double, Eigen::Dynamic, 2> &u)
const
{
    switch (rotation_) {
        case 0:
            return bicop_->hfunc2(cut_and_rotate(u));

        case 90:
            return 1.0 - bicop_->hfunc1(cut_and_rotate(u)).array();

        case 180:
            return 1.0 - bicop_->hfunc2(cut_and_rotate(u)).array();

        case 270:
            return bicop_->hfunc1(cut_and_rotate(u));

        default:
            throw std::runtime_error(std::string(
                "rotation can only take values in {0, 90, 180, 270}"
            ));
    }
}

//! calculates the inverse of \f$ h_1 f\f$ (see hfunc1()) w.r.t. the second
//! argument.
//! @param u \f$m \times 2\f$ matrix of evaluation points.
inline Eigen::VectorXd
Bicop::hinv1(const Eigen::Matrix<double, Eigen::Dynamic, 2> &u)
const
{
    switch (rotation_) {
        case 0:
            return bicop_->hinv1(cut_and_rotate(u));

        case 90:
            return bicop_->hinv2(cut_and_rotate(u));

        case 180:
            return 1.0 - bicop_->hinv1(cut_and_rotate(u)).array();

        case 270:
            return 1.0 - bicop_->hinv2(cut_and_rotate(u)).array();

        default:
            throw std::runtime_error(std::string(
                "rotation only takes value in {0, 90, 180, 270}"
            ));
    }
}

//! calculates the inverse of \f$ h_2 f\f$ (see hfunc2()) w.r.t. the first
//! argument.
//! @param u \f$m \times 2\f$ matrix of evaluation points.
inline Eigen::VectorXd
Bicop::hinv2(const Eigen::Matrix<double, Eigen::Dynamic, 2> &u)
const
{
    switch (rotation_) {
        case 0:
            return bicop_->hinv2(cut_and_rotate(u));

        case 90:
            return 1.0 - bicop_->hinv1(cut_and_rotate(u)).array();

        case 180:
            return 1.0 - bicop_->hinv2(cut_and_rotate(u)).array();

        case 270:
            return bicop_->hinv1(cut_and_rotate(u));

        default:
            throw std::runtime_error(std::string(
                "rotation only takes value in {0, 90, 180, 270}"
            ));
    }
}
//! @}


//! simulates from a bivariate copula.
//!
//! @param n number of observations.
//! @return An \f$ n \times 2 \f$ matrix of samples from the copula model.
inline Eigen::Matrix<double, Eigen::Dynamic, 2>
Bicop::simulate(const int &n) const
{
    Eigen::Matrix<double, Eigen::Dynamic, 2> U =
        tools_stats::simulate_uniform(n, 2);
    // use inverse Rosenblatt transform to generate a sample from the copula
    U.col(1) = hinv1(U);
    return U;
}

//! calculates the log-likelihood, defined as
//! \f[ \mathrm{loglik} = \sum_{i = 1}^n \ln c(U_{1, i}, U_{2, i}), \f]
//! where \f$ c \f$ is the copula density pdf().
//!
//! @param u \f$n \times 2\f$ matrix of observations.
inline double
Bicop::loglik(const Eigen::Matrix<double, Eigen::Dynamic, 2> &u) const
{
    return pdf(tools_eigen::nan_omit(u)).array().log().sum();
}

//! calculates the Akaike information criterion (AIC), defined as
//! \f[ \mathrm{AIC} = -2\, \mathrm{loglik} + 2 p, \f]
//! where \f$ \mathrm{loglik} \f$ is the log-liklihood and \f$ p \f$ is the
//! (effective) number of parameters of the model, see loglik() and
//! calculate_npars(). The AIC is a consistent model selection criterion
//! for nonparametric models.
//!
//! @param u \f$n \times 2\f$ matrix of observations.
inline double
Bicop::aic(const Eigen::Matrix<double, Eigen::Dynamic, 2> &u) const
{
    return -2 * loglik(u) + 2 * calculate_npars();
}

//! calculates the Bayesian information criterion (BIC), defined as
//! \f[ \mathrm{BIC} = -2\, \mathrm{loglik} +  \ln(n) p, \f]
//! where \f$ \mathrm{loglik} \f$ is the log-liklihood and \f$ p \f$ is the
//! (effective) number of parameters of the model, see loglik() and
//! calculate_npars(). The BIC is a consistent model selection criterion
//! for nonparametric models.
//!
//! @param u \f$n \times 2\f$ matrix of observations.
inline double
Bicop::bic(const Eigen::Matrix<double, Eigen::Dynamic, 2> &u) const
{
    return -2 * loglik(u) + calculate_npars() * log(static_cast<double>(u.rows()));
}

//! Returns the actual number of parameters for parameteric families. For
//! nonparametric families, there is a conceptually similar definition in
//! the sense that it can be used in the calculation of fit statistics.
inline double Bicop::calculate_npars() const
{
    return bicop_->calculate_npars();
}

//! converts a Kendall's \f$ \tau \f$ to the copula parameters of the
//! current family (only works for one-parameter families).
//!
//! @param tau a value in \f$ (-1, 1) \f$.
inline Eigen::MatrixXd Bicop::tau_to_parameters(const double &tau) const
{
    return bicop_->tau_to_parameters(tau);
}

//! converts the parameters to the Kendall's \f$ tau \f$ for the current
//! family (works for all families but `BicopFamily::tll`).
//!
//! @param parameters the parameters (must be a valid parametrization of
//!     the current family).
inline double
Bicop::parameters_to_tau(const Eigen::MatrixXd &parameters) const
{
    double tau = bicop_->parameters_to_tau(parameters);
    if (tools_stl::is_member(rotation_, {90, 270})) {
        tau *= -1;
    }
    return tau;
}

//! @name Getters and setters
//!
//! @{
inline BicopFamily Bicop::get_family() const
{
    return bicop_->get_family();
}

inline std::string Bicop::get_family_name() const
{
    return bicop_->get_family_name();
}

inline int Bicop::get_rotation() const
{
    return rotation_;
}

inline Eigen::MatrixXd Bicop::get_parameters() const
{
    return bicop_->get_parameters();
}

inline void Bicop::set_rotation(int rotation)
{
    check_rotation(rotation);
    rotation_ = rotation;
}

inline void Bicop::set_parameters(const Eigen::MatrixXd &parameters)
{
    bicop_->set_parameters(parameters);
}
//! @}


//! @name Utilities
//! @{
//! useful functions for bivariate copulas

//! adjust's the copula model to a change in the variable order.
inline void Bicop::flip()
{
    BicopFamily family = bicop_->get_family();
    if (tools_stl::is_member(family, bicop_families::flip_by_rotation)) {
        if (rotation_ == 90) {
            set_rotation(270);
        } else if (rotation_ == 270) {
            set_rotation(90);
        }
    } else {
        bicop_->flip();
    }
}

//! summarizes the model into a string (can be used for printing).
inline std::string Bicop::str() const
{
    std::stringstream bicop_str;
    bicop_str << get_family_name();
    if (get_rotation() != 0) {
        bicop_str << " " << get_rotation() << "°";
    }
    if (get_family() != BicopFamily::indep) {
        bicop_str << ", parameters = " << get_parameters();
    }

    return bicop_str.str().c_str();
}
//! @}

inline BicopPtr Bicop::get_bicop() const
{
    return bicop_;
}

//! fits a bivariate copula (with fixed family) to data.
//!
//! For parametric models, two different methods are available. `"mle"` fits
//! the parameters by maximum-likelihood. `"itau"` uses inversion of
//! Kendall's \f$ \tau \f$, but is only available for one-parameter families
//! and the Student t copula. For the latter, there is a one-to-one
//! transformation for the first parameter, the second is found by profile
//! likelihood optimization (with accuracy of at least 0.5). Nonparametric
//! families have specialized methods, no specification is required.
//!
//! @param data an \f$ n \times 2 \f$ matrix of observations contained in
//!     \f$(0, 1)^2 \f$.
//! @param controls the controls (see FitControlsBicop).
inline void Bicop::fit(const Eigen::Matrix<double, Eigen::Dynamic, 2> &data,
                       FitControlsBicop controls)
{
    std::string method;
    if (tools_stl::is_member(bicop_->get_family(),
                             bicop_families::parametric)) {
        method = controls.get_parametric_method();
    } else {
        method = controls.get_nonparametric_method();
    }
    bicop_->fit(tools_eigen::nan_omit(cut_and_rotate(data)), method,
                controls.get_nonparametric_mult());
}

//! selects the best fitting model, by calling fit() for all families in
//! `family_set` and selecting the best fitting model by either BIC or AIC,
//! see bic() and aic().
//!
//! @param data an \f$ n \times 2 \f$ matrix of observations contained in
//!     \f$(0, 1)^2 \f$.
//! @param controls the controls (see FitControlsBicop).
inline void Bicop::select(Eigen::Matrix<double, Eigen::Dynamic, 2> data,
                          FitControlsBicop controls)
{
    using namespace tools_select;
    data = tools_eigen::nan_omit(data);

    if (data.rows() < 10) {
        bicop_ = AbstractBicop::create();
        rotation_ = 0;
    } else {
        rotation_ = 0;
        data = cut_and_rotate(data);
        std::vector <Bicop> bicops = create_candidate_bicops(data, controls);

        // Estimate all models and select the best one using the
        // selection_criterion
        double fitted_criterion = 1e6;
        std::mutex m;
        auto fit_and_compare = [&](Bicop cop) {
            tools_interface::check_user_interrupt();

            // Estimate the model
            cop.fit(data, controls);

            // Compute the selection criterion
            double new_criterion;
            if (controls.get_selection_criterion() == "loglik") {
                new_criterion = -cop.loglik(data);
            } else if (controls.get_selection_criterion() == "aic") {
                new_criterion = cop.aic(data);
            } else {
                new_criterion = cop.bic(data);
            }

            // the following block modifies thread-external variables
            // and is thus shielded by a mutex
            {
                std::lock_guard <std::mutex> lk(m);
                // If the new model is better than the current one,
                // then replace the current model by the new one
                if (new_criterion < fitted_criterion) {
                    fitted_criterion = new_criterion;
                    bicop_ = cop.get_bicop();
                    rotation_ = cop.get_rotation();
                }
            }
        };

        tools_parallel::map_on_pool(fit_and_compare,
                                    bicops,
                                    controls.get_num_threads());
    }
}

//! Data manipulations for rotated families
//!
//! @param u \f$m \times 2\f$ matrix of data.
//! @return The manipulated data.
inline Eigen::Matrix<double, Eigen::Dynamic, 2> Bicop::cut_and_rotate(
    const Eigen::Matrix<double, Eigen::Dynamic, 2> &u) const
{
    Eigen::Matrix<double, Eigen::Dynamic, 2> u_new(u.rows(), 2);

    // counter-clockwise rotations
    switch (rotation_) {
        case 0:
            u_new = u;
            break;

        case 90:
            u_new.col(0) = u.col(1);
            u_new.col(1) = 1.0 - u.col(0).array();
            break;

        case 180:
            u_new.col(0) = 1.0 - u.col(0).array();
            u_new.col(1) = 1.0 - u.col(1).array();
            break;

        case 270:
            u_new.col(0) = 1.0 - u.col(1).array();
            u_new.col(1) = u.col(0);
            break;
    }

    // truncate to interval [eps, 1 - eps]
    Eigen::Matrix<double, Eigen::Dynamic, 2> eps =
        Eigen::Matrix<double, Eigen::Dynamic, 2>::Constant(u.rows(), 2, 1e-10);
    u_new = (1.0 - eps.array()).min(u_new.array());
    u_new = eps.array().max(u_new.array());

    return u_new;
}

inline void Bicop::check_rotation(int rotation) const
{
    using namespace tools_stl;
    std::vector<int> allowed_rotations = {0, 90, 180, 270};
    if (!is_member(rotation, allowed_rotations)) {
        throw std::runtime_error("rotation must be one of {0, 90, 180, 270}");
    }
    if (is_member(bicop_->get_family(), bicop_families::rotationless)) {
        if (rotation != 0) {
            throw std::runtime_error("rotation must be 0 for the " +
                                     bicop_->get_family_name() + " copula");
        }
    }
}
}
