// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/bicop/class.hpp>
#include <vinecopulib/bicop/tools_select.hpp>
#include <vinecopulib/misc/tools_stats.hpp>
#include <vinecopulib/misc/tools_stl.hpp>

//! Tools for bivariate and vine copula modeling
namespace vinecopulib
{
    
    //! creates the independence copula.
    Bicop::Bicop()
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
    Bicop::Bicop(BicopFamily family, int rotation,
                 const Eigen::MatrixXd& parameters)
    {
        bicop_ = AbstractBicop::create(family, parameters);
        // family must be set before checking the rotation
        set_rotation(rotation);
    }

    //! equivalent to `Bicop cop; cop.select(data, controls)`.
    //! @param data see select().
    //! @param controls see select().
    Bicop::Bicop(Eigen::Matrix<double, Eigen::Dynamic, 2> data,
                 FitControlsBicop controls)
    {
        select(data, controls);
    }

    //! creates from a boost::property_tree::ptree object
    //! @param input the boost::property_tree::ptree object to convert from
    //! (see to_ptree() for the structure of the input).
    Bicop::Bicop(boost::property_tree::ptree input) :
            Bicop(
                    get_family_enum(input.get<std::string>("family")),
                    input.get<int>("rotation"),
                    tools_serialization::ptree_to_matrix<double>(input.get_child("parameters"))
            ) {}

    //! creates from a JSON file
    //! @param filename the name of the JSON file to read (see to_ptree() for the
    //! structure of the file).
    Bicop::Bicop(const char *filename) :
            Bicop(tools_serialization::json_to_ptree(filename)) {}

    //! Convert the copula into a boost::property_tree::ptree object
    //!
    //! The boost::property_tree::ptree is contains of three values named
    //! `"family"`, `"rotation"`, `"parameters"`, respectively a string
    //! for the family name, an integer for the rotation, and an Eigen::MatrixXd
    //! for the parameters.
    //!
    //! @return the boost::property_tree::ptree object containing the copula.
    boost::property_tree::ptree Bicop::to_ptree()
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
    void Bicop::to_json(const char *filename)
    {
        boost::property_tree::write_json(filename, to_ptree());
    }

    //! evaluates the copula density.
    //!
    //! @param u \f$n \times 2\f$ matrix of evaluation points.
    //! @return The copula density evaluated at \c u.
    Eigen::VectorXd Bicop::pdf(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
    {
        Eigen::VectorXd f = bicop_->pdf(cut_and_rotate(u));
        f = f.unaryExpr([](const double x){ return std::min(x,1e16);});
        return f;
    }

    //! evaluates the copula distribution.
    //!
    //! @param u \f$n \times 2\f$ matrix of evaluation points.
    //! @return The copula distribution evaluated at \c u.
    Eigen::VectorXd Bicop::cdf(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
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
    Eigen::VectorXd Bicop::hfunc1(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
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
    Eigen::VectorXd Bicop::hfunc2(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
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
    Eigen::VectorXd Bicop::hinv1(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
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
    Eigen::VectorXd Bicop::hinv2(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
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
    Eigen::Matrix<double, Eigen::Dynamic, 2> Bicop::simulate(const int& n)
    {
        Eigen::Matrix<double, Eigen::Dynamic, 2> U =
                tools_stats::simulate_uniform(n, 2);
        // use inverse Rosenblatt transform to generate a sample from the copula
        U.col(1) = hinv1(U);
        return U;
    }

    //! calculates the log-likelihood.
    //!
    //! The log-likelihood is defined as
    //! \f[ \mathrm{loglik} = \sum_{i = 1}^n \ln c(U_{1, i}, U_{2, i}), \f]
    //! where \f$ c \f$ is the copula density pdf().
    //!
    //! @param u \f$n \times 2\f$ matrix of observations.
    double Bicop::loglik(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
    {
        return pdf(u).array().log().sum();
    }

    //! calculates the Akaike information criterion (AIC).
    //!
    //! The AIC is defined as
    //! \f[ \mathrm{AIC} = -2\, \mathrm{loglik} + 2 p, \f]
    //! where \f$ \mathrm{loglik} \f$ is the log-liklihood and \f$ p \f$ is the
    //! (effective) number of parameters of the model, see loglik() and
    //! calculate_npars(). The AIC is a consistent model selection criterion
    //! for nonparametric models.
    //!
    //! @param u \f$n \times 2\f$ matrix of observations.
    double Bicop::aic(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
    {
        return -2 * loglik(u) + 2 * calculate_npars();
    }

    //! calculates the Bayesian information criterion (BIC).
    //!
    //! The BIC is defined as
    //! \f[ \mathrm{BIC} = -2\, \mathrm{loglik} +  \ln(n) p, \f]
    //! where \f$ \mathrm{loglik} \f$ is the log-liklihood and \f$ p \f$ is the
    //! (effective) number of parameters of the model, see loglik() and
    //! calculate_npars(). The BIC is a consistent model selection criterion
    //! for nonparametric models.
    //!
    //! @param u \f$n \times 2\f$ matrix of observations.
    double Bicop::bic(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
    {
        return -2 * loglik(u) + calculate_npars() * log(u.rows());
    }

    //! calculates the effective number of parameters.
    //!
    //! Returns the actual number of parameters for parameteric families. For
    //! nonparametric families, there is a conceptually similar definition in
    //! the sense that it can be used in the calculation of fit statistics.
    double Bicop::calculate_npars()
    {
        return bicop_->calculate_npars();
    }

    //! converts a Kendall's \f$ \tau \f$ to the copula parameters of the
    //! current family (only works for one-parameter families).
    //!
    //! @param tau a value in \f$ (-1, 1) \f$.
    Eigen::MatrixXd Bicop::tau_to_parameters(const double& tau)
    {
        return bicop_->tau_to_parameters(tau);
    }

    //! converts the parameters to the Kendall's \f$ tau \f$ for the current
    //! family (works for all families but `BicopFamily::tll`).
    //!
    //! @param parameters the parameters (must be a valid parametrization of
    //!     the current family).
    double Bicop::parameters_to_tau(const Eigen::MatrixXd& parameters)
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
    BicopFamily Bicop::get_family() const
    {
        return bicop_->get_family();
    }

    std::string Bicop::get_family_name() const
    {
        return bicop_->get_family_name();
    };

    int Bicop::get_rotation() const
    {
        return rotation_;
    }

    Eigen::MatrixXd Bicop::get_parameters() const
    {
        return bicop_->get_parameters();
    }

    //! @param rotation
    void Bicop::set_rotation(int rotation) {
        check_rotation(rotation);
        rotation_ = rotation;
    }

    //! @param parameters
    void Bicop::set_parameters(const Eigen::MatrixXd& parameters)
    {
        bicop_->set_parameters(parameters);
    }
    //! @}


    //! @name Utilities
    //! @{
    //! adjust's the copula model to a change in the variable order.
    void Bicop::flip()
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
    std::string Bicop::str()
    {
        std::stringstream bicop_str;
        bicop_str << "family = "    << get_family_name() <<
                  ", rotation = "   << get_rotation() <<
                  ", parameters = " << get_parameters();

        return bicop_str.str().c_str();
    }
    //! @}

    BicopPtr Bicop::get_bicop()
    {
        return bicop_;
    };



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
    void Bicop::fit(const Eigen::Matrix<double, Eigen::Dynamic, 2> &data,
                    FitControlsBicop controls)
    {
        std::string method;
        if (tools_stl::is_member(bicop_->get_family(),
                                 bicop_families::parametric)) {
            method = controls.get_parametric_method();
        } else {
            method = controls.get_nonparametric_method();
        }
        bicop_->fit(cut_and_rotate(data), method,
                    controls.get_nonparametric_mult());
    }

    //! selects the best fitting model.
    //!
    //! The function calls fit() for all families in `family_set`)  and selects
    //! the best fitting model by either BIC or AIC, see bic() and aic().
    //!
    //! @param data an \f$ n \times 2 \f$ matrix of observations contained in
    //!     \f$(0, 1)^2 \f$.
    //! @param controls the controls (see FitControlsBicop).
    void Bicop::select(Eigen::Matrix<double, Eigen::Dynamic, 2> data,
                       FitControlsBicop controls)
    {
        using namespace tools_stl;
        std::vector<BicopFamily> family_set = controls.get_family_set();
        std::string parametric_method = controls.get_parametric_method();
        std::string nonparametric_method = controls.get_nonparametric_method();
        std::string method;
        double mult = controls.get_nonparametric_mult();
        std::string selection_criterion = controls.get_selection_criterion();
        bool preselect_families = controls.get_preselect_families();

        // If the familyset is empty, use all families.
        // If the familyset is not empty, check that all included families are implemented.
        if (family_set.empty()) {
            if (parametric_method == "itau") {
                family_set = bicop_families::itau;
            } else {
                family_set = bicop_families::all;
            }
        } else {
            if (intersect(family_set, bicop_families::all).empty()) {
                throw std::runtime_error("One of the families is not implemented");
            }
            if (parametric_method == "itau") {
                family_set = intersect(family_set, bicop_families::itau);
                if (family_set.empty()) {
                    throw std::runtime_error("No family with method itau provided");
                }
            }
        }

        // When using rotations, add only the ones that yield the appropriate
        // association direction.
        auto tau = tools_stats::pairwise_tau(data);
        std::vector<int> which_rotations;
        if (tau > 0) {
            which_rotations = {0, 180};
        } else {
            which_rotations = {90, 270};
        }

        std::vector<double> c(2);
        if (preselect_families) {
            rotation_ = 0;
            c = get_c1c2(cut_and_rotate(data), tau);
        }

        // Create the combinations of families and rotations to estimate
        std::vector<BicopFamily> families;
        std::vector<int> rotations;
        for (auto family : family_set) {
            bool is_rotationless = is_member(family,
                                             bicop_families::rotationless);
            bool preselect = true;
            if (is_rotationless) {
                if (preselect_families) {
                    preselect = preselect_family(c, tau, family,
                                                 0, is_rotationless);
                }
                if (preselect) {
                    families.push_back(family);
                    rotations.push_back(0);
                }
            } else {
                for (auto rotation : which_rotations) {
                    if (preselect_families) {
                        preselect = preselect_family(c, tau, family,
                                                     rotation, is_rotationless);
                    }
                    if (preselect) {
                        families.push_back(family);
                        rotations.push_back(rotation);
                    }
                }
            }
        }

        // Estimate all models and select the best one using the selection_criterion
        BicopPtr fitted_bicop;
        int fitted_rotation = 0;
        double fitted_criterion = 1e6;
        for (unsigned int j = 0; j < families.size(); j++) {
            // Estimate the model
            bicop_ = AbstractBicop::create(families[j]);
            rotation_ = rotations[j];
            if (tools_stl::is_member(bicop_->get_family(),
                                     bicop_families::parametric)) {
                method = parametric_method;
            } else {
                method = nonparametric_method;
            }
            bicop_->fit(cut_and_rotate(data), method, mult);

            // Compute the selection criterion
            double new_criterion;
            if (selection_criterion == "aic") {
                new_criterion = aic(data);
            } else if (selection_criterion == "bic") {
                new_criterion = bic(data);
            } else {
                throw std::runtime_error("Selection criterion not implemented");
            }

            // If the new model is better than the current one,
            // then replace the current model by the new one
            if (new_criterion < fitted_criterion) {
                fitted_criterion = new_criterion;
                fitted_bicop = bicop_;
                fitted_rotation = rotation_;
            }
        }

        bicop_ = fitted_bicop;
        rotation_ = fitted_rotation;
    }

    //! Data manipulations for rotated families
    //!
    //! @param u \f$m \times 2\f$ matrix of data.
    //! @return The manipulated data.
    Eigen::Matrix<double, Eigen::Dynamic, 2> Bicop::cut_and_rotate(
        const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
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
        u_new = u_new.array().min(1.0 - eps.array());
        u_new = u_new.array().max(eps.array());

        return u_new;
    }

    void Bicop::check_rotation(int rotation)
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
