#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;
using namespace arma;

//' Minus Log Likelihood iAR-T Model
//'
//' This function return the negative log likelihood of the iAR-T given specific values of phi and sigma.
//'
//' @param coef An array with the parameters of the iAR-T model. The first element of the array corresponding to the phi parameter and the second element to the scale parameter sigma
//' @param series Array with the time series observations
//' @param times Array with the irregular observational times
//' @param df degrees of freedom
//' @param yest The estimate of a missing value in the time series. This function recognizes a missing value with a NA. If the time series does not have a missing value, this value does not affect the computation of the likelihood.
//'
//' @return Value of the negative log likelihood evaluated in phi,sigma and df.
//' @references
//' \insertRef{Eyheramendy_2018}{iAR}
//'
//' @seealso
//'
//' \code{\link{gentime}}, \code{\link{iARtsample}}
//'
//' @keywords internal
//' @examples
//' \dontshow{
//' n=300
//' set.seed(6714)
//' st<-gentime(n)
//' y<-iARtsample(coef=0.9,times=st,sigma=1,df=3)
//' iARphit(coef=c(0.9,1),series=y$series,times=st,yest=0)
//' }
//' @noRd
// [[Rcpp::export]]
double iARphit(arma::vec yest,arma::vec coef, arma::vec series, arma::vec times, double df = 3) {
  arma::vec y = series;
  if(y.has_nan() == true) {
     y.elem(arma::find_nonfinite(y))=yest;
  }
  double sigma = arma::as_scalar(coef.row(1));
  int n = y.size();
  arma::vec d = arma::diff(times);

  arma::vec xd(n-1, fill::zeros);
  for(int i = 0; i < n-1; ++i) {
    xd[i] = std::pow(arma::as_scalar(coef.row(0)), d[i]);
  }

  arma::vec yhat(n-1, fill::zeros);
  yhat = xd % y.rows(0, n-2);

  arma::vec gL(n-1, fill::zeros);
  gL = sigma * (1 - (arma::pow(xd, 2))) * ((df-2)/df);

  double temp1 = std::tgamma((df+1)/2);
  double temp2 = std::tgamma(df/2)*(std::sqrt(df*datum::pi));
  double cte = (n-1) * std::log(temp1/temp2);

  arma::vec stand(n-1, fill::zeros);
  stand = (y.rows(1, n-1) - yhat)/(arma::sqrt(gL));
  stand = arma::pow(stand, 2);

  double s1 = arma::sum(0.5 * arma::log(gL));
  double s2 = arma::sum(arma::log(1 + (1/df)*stand));

  double out = cte - s1 - ((df+1)/2) * s2 - 0.5 * (std::log(2*datum::pi) + std::pow(arma::as_scalar(y.row(0)), 2));

  return -out;
}
