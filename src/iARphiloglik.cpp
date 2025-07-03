#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;
using namespace arma;

//' Minus Log Likelihood of the iAR Model
//'
//' This function return the negative log likelihood of the iAR Model for a specific value of phi.
//'
//' @param coef A given phi coefficient of the iAR model.
//' @param series Array with the time series observations.
//' @param times Array with the irregular observational times.
//' @param series_esd Array with the measurements error standard deviations.
//' @param zero_mean logical; if TRUE, the array series has zero mean; if FALSE, series has a mean different from zero.
//' @param standardized logical; if TRUE, the array series was standardized; if FALSE, series contains the raw data
//'
//' @return Value of the negative log likelihood evaluated in phi.
//' @references
//' \insertRef{Eyheramendy_2018}{iAR}
//'
//' @seealso
//'
//' \code{\link{gentime}}, \code{\link{iARsample}}
//'
//' @keywords internal
//' @examples
//' \dontshow{
//' set.seed(6714)
//' st<-gentime(n=100)
//' y<-iARsample(coef=0.99,times=st)
//' y<-y$series
//' iARphiloglik(coef=0.8,series=y,times=st,series_esd=0)
//' }
//' @noRd
// [[Rcpp::export]]
double iARphiloglik(double coef, arma::vec series, arma::vec times, arma::vec series_esd, bool zero_mean = true, bool standardized = true) {
  arma::vec y = series;
  double x = coef;
  double sigma = 1;
  int mu = 0;

  int n = y.size();
  arma::vec delta(n - 1, fill::zeros);

  if(arma::sum(series_esd) != 0) {
    delta = series_esd.rows(1, n-1);
  }

  if(standardized==false) {
    sigma = arma::var(y);
  }

  if(zero_mean == false) {
    mu = arma::mean(y);
  }

  arma::vec d = arma::diff(times);
  arma::vec phi(n-1, fill::zeros);
  for(int i = 0; i < n-1; ++i) {
    phi[i] = std::pow(x, d[i]);
  }

  arma::vec yhat(n-1, fill::zeros);
  for(int i = 0; i < n-1; ++i) {
    yhat[i] = mu + phi[i] * (y[i] - mu);
  }

  double cte = (n/2) * std::log(2*datum::pi);

  arma::vec temp0 = sigma * (1 - arma::pow(phi, 2)) + arma::pow(delta, 2);
  arma::vec temp1 = arma::log(temp0);
  arma::vec temp2 = (y.rows(1, n-1) - yhat);
  arma::vec temp3 = arma::pow(temp2, 2);

  return cte + 0.5 * arma::sum(temp1 + temp3/temp0);
}
