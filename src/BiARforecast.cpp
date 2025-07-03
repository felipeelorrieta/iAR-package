#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

//' Forecast from BiAR model
//'
//' Forecast from models fitted by \code{\link{BiARkalman}}
//'
//' @param coef An array with the parameters of the BiAR model. The elements of the array are, in order, the autocorrelation and the cross-correlation coefficient of the BiAR model.
//' @param series1 Array with the observations of the first time series of the BiAR process.
//' @param series2 Array with the observations of the second time series of the BiAR process.
//' @param times Array with the observational times.
//' @param tAhead The time ahead for which the forecast is required.
//'
//' @return A list with the following components:
//' \itemize{
//' \item{forecast}{ Point forecast in the time ahead required.}
//' }
//' @references
//' \insertRef{Elorrieta_2021}{iAR}
//'
//' @seealso
//'
//' \code{\link{BiARsample}}, \code{\link{BiARkalman}}, \code{\link{BiARfit}}
//'
//'
//' @keywords internal
//' @examples
//' \dontshow{
//' #Simulated Data
//' n=100
//' set.seed(6714)
//' st<-gentime(n)
//' x=BiARsample(phiR=0.9,phiI=0.3,times=st)
//' biar=iAR::BiARkalman(series1=x$series[1,],series2=x$series[2,],times=st)
//' forBiAR<-BiARforecast(coef=c(biar$phiR,biar$phiI),series1=x$series[1,],series2=x$series[2,],
//' times=st,tAhead=c(1.3))
//' }
//' @noRd
// [[Rcpp::export]]
List BiARforecast(arma::vec coef, arma::vec series1, arma::vec series2, arma::vec times, double tAhead) {
  List output;

  arma::vec y1=series1;
  arma::vec y2=series2;
  arma::cx_double phi(coef[0], coef[1]);
  double phiMod = std::sqrt(std::pow(phi.real(), 2.0) + std::pow(phi.imag(), 2.0));
  if(phiMod >= 1.0) {
    return output;
  }

  int n = y1.size();
  arma::mat series_esd1(1, n, fill::zeros);
  arma::mat series_esd2(1, n, fill::zeros);

  Function BiARfit("BiARfit", Environment::namespace_env("iAR"));
  List biarFitOutput = BiARfit(coef, y1, y2, times, series_esd1, series_esd2);

  arma::mat yhat = biarFitOutput["fitted"];
  arma::mat xhat = biarFitOutput["fitted.state"];
  arma::mat theta = biarFitOutput["Theta"];
  arma::mat Lambda = biarFitOutput["Lambda"];
  arma::mat Sighat = biarFitOutput["Sighat"];
  arma::mat Qt = biarFitOutput["Qt"];


  arma::mat G = arma::eye(2,2);
  arma::mat y = arma::join_vert(y1.t(), y2.t());

  arma::mat xhat1(2, 1, fill::zeros);
  arma::mat yhat1(2, 1, fill::zeros);
  arma::mat sigHat2;
  arma::mat Lambda2;

  double psi = (phi.imag() < 0 && phiMod < 1) ? -acos(phi.real()/phiMod) : acos(phi.real()/phiMod);

  double phi2R = pow(phiMod, tAhead) * cos(tAhead * psi);
  double phi2I = pow(phiMod, tAhead) * sin(tAhead * psi);

  arma::mat F = {{phi2R, -phi2I}, {phi2I, phi2R}};
  arma::mat temp = y.col(n-1) - (G * xhat.col(n-1));
  xhat1 = F * xhat.col(n-1) + theta * arma::inv(Lambda) * temp;
  yhat1 = G * xhat1;

  sigHat2 = F * Sighat * F.t() + Qt - theta * arma::inv(Lambda) * theta.t();  
  Lambda2 = G * sigHat2 * G.t();

  output["forecast"] = yhat1;

  return output;
}
