#include <RcppArmadillo.h>
#include <string.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

//' Fitted Values of CiAR model
//'
//' Fit a CiAR model to an irregularly observed time series.
//'
//' @param coef An array with the parameters of the CiAR model. The elements of the array are, in order, the real and the imaginary part of the phi parameter of the CiAR model.
//' @param series Array with the time series observations.
//' @param times Array with the irregular observational times.
//' @param standardized logical; if TRUE, the array series is standardized; if FALSE, series contains the raw time series
//' @param zero_mean logical; if TRUE, the array series has zero mean; if FALSE, series has a mean different from zero.
//' @param c Nuisance parameter corresponding to the variance of the imaginary part.
//'
//' @return A list with the following components:
//' \item{fitted}{ Fitted values of the observable part of CiAR model.}
//' \item{xhat}{ Fitted values of both observable part and imaginary part of CiAR model.}
//' \item{Lambda}{ Lambda value estimated by the CiAR model at the last time point.}
//' \item{Theta}{ Theta array estimated by the CiAR model at the last time point.}
//' \item{Sighat}{ Covariance matrix estimated by the CiAR model at the last time point.}
//' \item{Qt}{ Covariance matrix of the state equation estimated by the CiAR model at the last time point.}
//' @references
//' \insertRef{Elorrieta_2019}{iAR}
//' 
//' @keywords internal
//' @examples
//' \dontshow{
//' n=100
//' set.seed(6714)
//' }
//' @noRd
// [[Rcpp::export]]
List CiARfit(arma::vec coef, arma::vec series, arma::vec times, bool zero_mean=true, bool standardized=true, double c=1) {
  List output;
  
  arma::vec y=series;
  arma::cx_double phi(coef[0], coef[1]);
  double phiMod = std::sqrt(std::pow(phi.real(), 2.0) + std::pow(phi.imag(), 2.0));
  if(phiMod >= 1.0) {
    return output;
  }

  double sigmaY = 1;
  double mu=0;
  if(standardized == false) {
    sigmaY = arma::var(y);
  }
  if(zero_mean == false) {
    mu = arma::mean(y);
    y = y - mu;
  }
  if(standardized == false) {
    y = y/std::sqrt(sigmaY); 
  }
  int n = y.size();
  arma::mat auxSigmaHat = {{1, 0}, {0, c}};
  arma::mat sigmaHat = sigmaY * auxSigmaHat;

  arma::mat xhat(2, n, fill::zeros);
  arma::vec delta = arma::diff(times);

  arma::mat Q = sigmaHat;
  arma::mat F(2, 2, fill::zeros);
  arma::mat G = arma::eye(1,2);
  G.col(0) = 1;

  double psi = -acos(phi.real()/phiMod);

  arma::mat theta;
  arma::mat Lambda(2, 2, fill::zeros);
  arma::mat Qt;

  for (int i = 0; i < (n-1); ++i) {
    Lambda = G * sigmaHat * G.t();

    if(arma::as_scalar(Lambda) <= 0 || Lambda.has_nan()) {
      break;
    }

    double deltaPsi = delta[i] * psi;
    double phiModDelta = pow(phiMod, delta[i]);

    double phi2Real = phiModDelta * cos(deltaPsi);
    double phi2Imag = phiModDelta * sin(deltaPsi);

    arma::mat F = {{phi2Real, -phi2Imag}, {phi2Imag, phi2Real}};

    arma::cx_double phi2Inter = pow(phi, delta[i]);
    double phi2Mod = std::pow(std::sqrt(std::pow(phi2Inter.real(), 2.0) + std::pow(phi2Inter.imag(), 2.0)), 2.0);
    auto phi2 = 1 - phi2Mod;

    Qt = phi2 * Q;
    theta = F * sigmaHat * G.t();

    arma::mat aux = y[i] - (G * xhat.col(i));
    xhat.col(i + 1) = F * xhat.col(i) + theta * arma::inv(Lambda) * aux;

    sigmaHat = F * sigmaHat * F.t() + Qt - theta * arma::inv(Lambda) * theta.t();
  }
  
  xhat.row(0) = xhat.row(0)*std::sqrt(sigmaY)+mu;
  arma::mat yhat = (G * xhat);

  output["fitted"] = yhat.row(0);
  output["xhat"] = xhat;
  output["Sighat"] = sigmaHat;
  output["Theta"] = theta;
  output["Lambda"] = Lambda;
  output["Qt"] = Qt;

  return output;
}
