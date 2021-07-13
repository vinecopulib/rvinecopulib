// Copyright © 2016-2021 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <boost/math/special_functions/expm1.hpp>
#include <boost/math/special_functions/log1p.hpp>
#include <vinecopulib/misc/tools_eigen.hpp>

namespace vinecopulib {

inline FrankBicop::FrankBicop()
{
  family_ = BicopFamily::frank;
  parameters_ = Eigen::VectorXd(1);
  parameters_lower_bounds_ = Eigen::VectorXd(1);
  parameters_upper_bounds_ = Eigen::VectorXd(1);
  parameters_ << 0;
  parameters_lower_bounds_ << -35;
  parameters_upper_bounds_ << 35;
}

inline double
FrankBicop::generator(const double& u)
{
  double theta = double(this->parameters_(0));
  return -std::log(boost::math::expm1(-theta * u) / boost::math::expm1(-theta));
}

inline double
FrankBicop::generator_inv(const double& u)
{
  double theta = double(this->parameters_(0));
  return -boost::math::log1p(boost::math::expm1(-theta) * std::exp(-u)) / theta;
}

inline double
FrankBicop::generator_derivative(const double& u)
{
  double theta = double(this->parameters_(0));
  return -theta / boost::math::expm1(theta * u);
}

// inline double FrankBicop::generator_derivative2(const double &u)
//{
//    double theta = double(this->parameters_(0));
//    return std::pow(theta, 2) /
//           std::pow(boost::math::expm1(theta * u) * std::exp(-theta * u / 2),
//           2);
//}

inline Eigen::VectorXd
FrankBicop::pdf_raw(const Eigen::MatrixXd& u)
{
  double theta = static_cast<double>(parameters_(0));
  auto f = [theta](const double& u1, const double& u2) {
    return (theta * std::expm1(theta) *
            std::exp(theta * u2 + theta * u1 + theta)) /
           std::pow(std::exp(theta * u2 + theta * u1) -
                      std::exp(theta * u2 + theta) -
                      std::exp(theta * u1 + theta) + std::exp(theta),
                    2.0);
  };
  return tools_eigen::binaryExpr_or_nan(u, f);
}

inline Eigen::MatrixXd
FrankBicop::tau_to_parameters(const double& tau)
{
  Eigen::VectorXd tau0 = Eigen::VectorXd::Constant(1, tau);
  auto f = [&](const Eigen::VectorXd& par) {
    return Eigen::VectorXd::Constant(1, parameters_to_tau(par));
  };
  return tools_eigen::invert_f(tau0,
                               f,
                               parameters_lower_bounds_(0) + 1e-6,
                               parameters_upper_bounds_(0) - 1e-5);
}

inline double
FrankBicop::parameters_to_tau(const Eigen::MatrixXd& parameters)
{
  double par = std::fabs(parameters(0));
  if (par < 1e-5) {
    return 0.0;
  }
  double tau = 1 - 4 / par + (4 / par) * debye1(par) / par;
  return (parameters(0) >= 0) ? tau : -tau;
}

inline Eigen::VectorXd
FrankBicop::get_start_parameters(const double tau)
{
  Eigen::VectorXd par = tau_to_parameters(tau);
  par = par.cwiseMax(parameters_lower_bounds_);
  par = par.cwiseMin(parameters_upper_bounds_);
  return par;
}

//! @brief computes the Debye function of order 1.
//! @param x the argument and upper limit of the integral. x>=0.
//! @return the Debye function. Zero if x<=0.
//! @details The Debye function is defined as \f$ \int_0^x t  / [exp(t)-1] dt
//! \f$. Code modified from implementation in
//! https://github.com/openturns/openturns.
inline double
debye1(const double& x)
{
  if (x <= 0.)
    return 0.;
  double m_1_2pi = .159154943091895335768883763373; // 1/(2pi)

  if (x >= 3.) {
    double sum = 1.64493406684822643647241516665e+00;
    static constexpr std::array<short, 14> kLim{ 0, 0, 0, 13, 10, 8, 7,
                                                 6, 5, 5, 4,  4,  4, 3 };
    const short kmax = (x < 14) ? kLim[static_cast<int>(x)] : 3;
    for (short k = 1; k <= kmax; k++) {
      const double xk = x * k;
      double ksum = 1. / xk;
      ksum += ksum / xk;
      sum -= exp(-xk) * ksum * x * x;
    }
    return sum;
  } else {
    static constexpr std::array<double, 71> koeff = {
      0.,
      1.289868133696452872944830333292e+00,
      1.646464674222763830320073930823e-01,
      3.468612396889827942903585958184e-02,
      8.154712395888678757370477017305e-03,
      1.989150255636170674291917800638e-03,
      4.921731066160965972759960954793e-04,
      1.224962701174096585170902102707e-04,
      3.056451881730374346514297527344e-05,
      7.634586529999679712923289243879e-06,
      1.907924067745592226304077366899e-06,
      4.769010054554659800072963735060e-07,
      1.192163781025189592248804158716e-07,
      2.980310965673008246931701326140e-08,
      7.450668049576914109638408036805e-09,
      1.862654864839336365743529470042e-09,
      4.656623667353010984002911951881e-10,
      1.164154417580540177848737197821e-10,
      2.910384378208396847185926449064e-11,
      7.275959094757302380474472711747e-12,
      1.818989568052777856506623677390e-12,
      4.547473691649305030453643155957e-13,
      1.136868397525517121855436593505e-13,
      2.842170965606321353966861428348e-14,
      7.105427382674227346596939068119e-15,
      1.776356842186163180619218277278e-15,
      4.440892101596083967998640188409e-16,
      1.110223024969096248744747318102e-16,
      2.775557561945046552567818981300e-17,
      6.938893904331845249488542992219e-18,
      1.734723476023986745668411013469e-18,
      4.336808689994439570027820336642e-19,
      1.084202172491329082183740080878e-19,
      2.710505431220232916297046799365e-20,
      6.776263578041593636171406200902e-21,
      1.694065894509399669649398521836e-21,
      4.235164736272389463688418879636e-22,
      1.058791184067974064762782460584e-22,
      2.646977960169798160618902050189e-23,
      6.617444900424343177893912768629e-24,
      1.654361225106068880734221123349e-24,
      4.135903062765153408791935838694e-25,
      1.033975765691286264082026643327e-25,
      2.584939414228213340076225223666e-26,
      6.462348535570530772269628236053e-27,
      1.615587133892632406631747637268e-27,
      4.038967834731580698317525293132e-28,
      1.009741958682895139216954234507e-28,
      2.524354896707237808750799932127e-29,
      6.310887241768094478219682436680e-30,
      1.577721810442023614704107565240e-30,
      3.944304526105059031370476640000e-31,
      9.860761315262647572437533499000e-32,
      2.465190328815661892443976898000e-32,
      6.162975822039154730370601500000e-33,
      1.540743955509788682510501190000e-33,
      3.851859888774471706184973900000e-34,
      9.629649721936179265360991000000e-35,
      2.407412430484044816328953000000e-35,
      6.018531076210112040809600000000e-36,
      1.504632769052528010200750000000e-36,
      3.761581922631320025497600000000e-37,
      9.403954806578300063715000000000e-38,
      2.350988701644575015901000000000e-38,
      5.877471754111437539470000000000e-39,
      1.469367938527859384580000000000e-39,
      3.673419846319648458500000000000e-40,
      9.183549615799121117000000000000e-41,
      2.295887403949780249000000000000e-41,
      5.739718509874450320000000000000e-42,
      1.434929627468612270000000000000e-42
    };

    double sum = 0, sumold = 1;
    const double x2pi = x * m_1_2pi;
    for (short k = 1; (k < 70) && (sum != sumold); k++) {
      sumold = sum;
      sum += (2. + koeff[k]) * std::pow(x2pi, 2. * k) / (2 * k + 1.);
      k++;
      sum -= (2. + koeff[k]) * std::pow(x2pi, 2. * k) / (2 * k + 1.);
    }
    return x * (sum + 1. - x / 4.);
  }
}
}
