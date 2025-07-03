#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

//' Forecast from CiAR model
//'
//' Forecast from models fitted by \code{\link{CiARkalman}}
//'
//' @param coef An array with the parameters of the CiAR model. The elements of the array are, in order, the real and the imaginary part of the phi parameter of the CiAR model.
//' @param series Array with the time series observations.
//' @param times Array with the observational times.
//' @param tAhead The time ahead for which the forecast is required.
//' @param standardized logical; if TRUE, the array series is standardized; if FALSE, series contains the raw time series
//' @param zero_mean logical; if TRUE, the array series has zero mean; if FALSE, series has a mean different from zero.
//'
//' @return A list with the following components:
//' \itemize{
//' \item{forecast}{ Point forecast in the time ahead required.}
//' }
//' @references
//' \insertRef{Elorrieta_2019}{iAR}
//'
//' @seealso
//'
//' \code{\link{CiARsample}}, \code{\link{CiARkalman}}, \code{\link{CiARfit}}
//'
//'
//' @keywords internal
//' @examples
//' \dontshow{
//' #Simulated Data
//' n=100
//' set.seed(6714)
//' st<-gentime(n)
//' x=CiARsample(phiR=0.9,phiI=0,times=st,c=1)
//' y=x$series
//' y1=y/sd(y)
//' n=length(y1)
//' p=trunc(n*0.99)
//' ytr=y1[1:p]
//' yte=y1[(p+1):n]
//' str=st[1:p]
//' ste=st[(p+1):n]
//' tahead=ste-str[p]
//'
//' ciar=CiARkalman(series=ytr,times=str)
//' forCIAR<-CiARforecast(coef=c(ciar$phiR,ciar$phiI),series=ytr,times=str,tAhead=tahead)
//' }
//' @noRd
// [[Rcpp::export]]
List CiARforecast(arma::vec coef, arma::vec series, arma::vec times, double tAhead, bool zero_mean=true, bool standardized=true) {
  List output;

  arma::vec y1=series;
  arma::cx_double phi(coef[0], coef[1]);
  double phiMod = std::sqrt(std::pow(phi.real(), 2.0) + std::pow(phi.imag(), 2.0));
  if(phiMod >= 1.0) {
    return output;
  }

  Function CiARfit("CiARfit", Environment::namespace_env("iAR"));
  List ciarFitOutput = CiARfit(coef, y1, times, zero_mean,standardized);
  
  double sigmaY = 1;
  double mu=0;
  if(standardized == false) {
    sigmaY = arma::var(y1);
  }
  if(zero_mean == false) {
    mu = arma::mean(y1);
    y1 = y1 - mu;
  }
  if(standardized == false) {
    y1 = y1/std::sqrt(sigmaY); 
  }
  arma::vec yhat = ciarFitOutput["fitted"];
  arma::mat xhat = ciarFitOutput["xhat"];
  arma::vec Theta = ciarFitOutput["Theta"];
  arma::vec Lambda = ciarFitOutput["Lambda"];
  arma::mat Sighat = ciarFitOutput["Sighat"];
  arma::mat Qt = ciarFitOutput["Qt"];

  int n = yhat.size();
  arma::mat G(1, 2, fill::zeros);
  G.row(0).col(0) = 1;

  double psi = -acos(phi.real()/phiMod);

  double phi2R = pow(phiMod, tAhead) * cos(tAhead * psi);
  double phi2I = pow(phiMod, tAhead) * sin(tAhead * psi);

  double yhat1 = 0;
  arma::mat xhat1(2, 1, fill::zeros);
  arma::mat sigHat2;

  arma::mat F = {{phi2R, -phi2I}, {phi2I, phi2R}};
  arma::vec temp = y1[n-1] - G * xhat.col(n-1);

  xhat1 = (F * xhat.col(n-1)) + Theta * arma::inv(Lambda) * temp;
  xhat1.row(0) = xhat1.row(0)*std::sqrt(sigmaY)+mu;
  
  yhat1 = arma::as_scalar(G * xhat1);

  output["forecast"] = yhat1;

  return output;
}
