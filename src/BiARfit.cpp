#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

//' Fitted Values of BiAR model
//'
//' Fit a BiAR model to a bivariate irregularly observed time series.
//'
//' @param coef An array with the parameters of the BiAR model. The elements of the array are, in order, the autocorrelation and the cross correlation parameter of the BiAR model.
//' @param series1 Array with the observations of the first time series of the BiAR process.
//' @param series2 Array with the observations of the second time series of the BiAR process.
//' @param times Array with the irregular observational times.
//' @param series_esd1 Array with the measurements error standard deviations of the first time series of the BiAR process.
//' @param series_esd2 Array with the measurements error standard deviations of the second time series of the BiAR process.
//' @param zero_mean logical; if TRUE, the array y has zero mean; if FALSE, y has a mean different from zero.
//'
//' @return A list with the following components:
//' \itemize{
//' \item{rho}{ Estimated value of the contemporary correlation coefficient.}
//' \item{innov.var}{ Estimated value of the innovation variance.}
//' \item{fitted}{ Fitted values of the BiAR model.}
//' \item{fitted.state}{ Fitted state values of the BiAR model.}
//' \item{Lambda}{ Lambda value estimated by the BiAR model at the last time point.}
//' \item{Theta}{ Theta array estimated by the BiAR model at the last time point.}
//' \item{Sighat}{ Covariance matrix estimated by the BiAR model at the last time point.}
//' \item{Qt}{ Covariance matrix of the state equation estimated by the BiAR model at the last time point.}
//' }
//' @references
//' \insertRef{Elorrieta_2021}{iAR}
//' @seealso
//' \code{\link{gentime}}, \code{\link{BiARsample}}, \code{\link{BiARphikalman}}, \code{\link{BiARkalman}}
//'
//' @keywords internal
//' @examples
//' \dontshow{
//' n=80
//' set.seed(6714)
//' st<-gentime(n)
//' x=BiARsample(phiR=0.9,phiI=0.3,times=st,rho=0.9)
//' y=x$series
//' y1=y/apply(y,1,sd)
//' yerr1=rep(0,n)
//' yerr2=rep(0,n)
//' biar=BiARkalman(series1=y1[1,],series2=y1[2,],times=st,series_esd1 = yerr1,series_esd2=yerr2)
//' biar
//' predbiar=BiARfit(coef=c(biar$phiR,biar$phiI),series1=y1[1,],series2=y1[2,],times=st,series_esd1 = 
//' rep(0,length(y[1,])),series_esd2=rep(0,length(y[1,])))
//' rho=predbiar$rho
//' print(rho)
//' yhat=predbiar$fitted
//' }
//' @noRd
// [[Rcpp::export]]
List BiARfit(arma::vec coef, arma::vec series1, arma::vec series2, arma::vec times, arma::vec series_esd1, arma::vec series_esd2, bool zero_mean = true) {
  List output;
  
  arma::vec y1=series1;
  arma::vec y2=series2;
  arma::cx_double phi(coef[0], coef[1]);

  double phiMod = std::sqrt(std::pow(phi.real(), 2.0) + std::pow(phi.imag(), 2.0));
  if(phiMod >= 1.0) {
    return output;
  }

  if(zero_mean == FALSE) {
    y1 = y1 - arma::mean(y1);
    y2 = y2 - arma::mean(y2);
  }

  int n = y1.size();
  arma::mat sigmaY = arma::cov(arma::join_horiz(y1, y2));
  arma::mat sigmaHat = sigmaY * arma::eye(2,2);
  arma::mat xhat(2, n, fill::zeros);
  arma::vec delta = arma::diff(times);

  arma::mat Q = sigmaHat;
  arma::mat F(2, 2, fill::zeros);
  arma::mat G = arma::eye(2,2);

  double psi = (phi.imag() < 0 && phiMod < 1) ? -acos(phi.real()/phiMod) : acos(phi.real()/phiMod);

  arma::mat y = arma::join_vert(y1.t(), y2.t());

  arma::mat theta;
  arma::mat Lambda(2, 2, fill::zeros);
  arma::mat Qt;

  for (int i = 0; i < (n-1); ++i) {
    arma::mat R = {{std::pow(series_esd1[i+1], 2.0), 0}, {0, std::pow(series_esd2[i+1], 2.0)}};
    Lambda = G * sigmaHat * G.t() + R;

    if(det(Lambda) <= 0 || Lambda.has_nan()) {
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

    arma::mat M = phi2 * arma::eye(2,2);
    Qt = M * (Q - R);
    theta = F * sigmaHat * G.t();
    arma::mat aux = y.col(i) - (G * xhat.col(i));

    xhat.col(i + 1) = F * xhat.col(i) + theta * arma::inv(Lambda) * aux;

    if(arma::accu(R) == 0) {
      sigmaHat = Qt + F * (sigmaHat - sigmaHat.t()) * F.t();
    } else {
      sigmaHat = F * sigmaHat * F.t() + Qt - theta * arma::inv(Lambda) * theta.t();
    }
  }

  arma::mat yhat = G * xhat;

  arma::mat finalMat = y - G * xhat;
  arma::mat finalCor = arma::cor(finalMat.t());
  arma::mat finalCov = arma::cov(finalMat.t());

  output["rho"] = finalCor.at(0,1);
  output["innov.var"] = finalCov;
  output["fitted"] = yhat;
  output["fitted.state"] = xhat;
  output["Sighat"] = sigmaHat;
  output["Theta"] = theta;
  output["Lambda"] = Lambda;
  output["Qt"] = Qt;

  return output;
}
