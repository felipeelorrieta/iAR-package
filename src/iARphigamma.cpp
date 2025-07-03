#include <RcppArmadillo.h>
#include <string.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

//' Minus Log Likelihood iAR-Gamma Model
//'
//' This function return the negative log likelihood of the iAR-Gamma given specific values of phi, mean and sigma.
//'
//' @param coef An array with the parameters of the iAR-Gamma model. The first element of the array corresponding to the phi parameter, the second to the mean parameter, and the last one to the scale parameter sigma.
//' @param series Array with the time series observations.
//' @param times Array with the irregular observational times.
//' @param yest The estimate of a missing value in the time series. This function recognizes a missing value with a NA. If the time series does not have a missing value, this value does not affect the computation of the likelihood.
//'
//' @return Value of the negative log likelihood evaluated in phi, mean and sigma.
//' @references
//' \insertRef{Eyheramendy_2018}{iAR}
//'
//' @seealso
//'
//' \code{\link{gentime}}, \code{\link{iARgsample}}
//'
//' @keywords internal
//' @examples
//' \dontshow{
//' n=100
//' set.seed(6714)
//' st<-gentime(n)
//' y<-iARgsample(coef=0.9,times=st,sigma=1,mean=1)
//' iARphigamma(coef=c(0.9,1,1),series=y$series,times=st,yest=0)
//' }
//' @noRd
// [[Rcpp::export]]
double iARphigamma(arma::vec yest,arma::vec coef, arma::vec series, arma::vec times) {
  arma::vec y = series;
  if(y.has_nan() == true) {
     y.elem(arma::find_nonfinite(y))=yest;
  }
  double mu = arma::as_scalar(coef.row(1));
  double sigma = arma::as_scalar(coef.row(2));

  int n = y.size();
  arma::vec d = arma::diff(times);

  arma::vec xd(n-1, fill::zeros);
  for(int i = 0; i < n-1; ++i) {
    xd[i] = std::pow(arma::as_scalar(coef.row(0)), d[i]);
  }

  arma::vec yhat(n-1, fill::zeros);
  yhat = mu + xd % y.rows(0, n-2);

  arma::vec gL(n-1, fill::zeros);
  gL = sigma * (1 - arma::pow(xd,2));

  arma::vec beta(n-1, fill::zeros);
  beta=gL/yhat;

  arma::vec alpha(n-1, fill::zeros);
  for(int i = 0; i < (n-1); ++i) {
    alpha[i] = std::pow(yhat[i],2)/gL[i];
  }

  arma::vec temp0 = (-alpha)%arma::log(beta);
  arma::vec temp1 = (alpha - 1) % arma::log(y.rows(1,n-1));

  double out = arma::sum(temp0 - arma::lgamma(alpha) - (y.rows(1,n-1)/beta) + temp1) - arma::as_scalar(y.row(0));
  return -out;
}
