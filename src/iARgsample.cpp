#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;
using namespace arma;

//' Simulate from an iAR-Gamma Model
//'
//' Simulates an iAR-Gamma Time Series Model.
//'
//' @param coef A coefficient of iAR-Gamma model. A value between 0 and 1.
//' @param times Array with observational times.
//' @param sigma Scale parameter of the IAR-Gamma process. A positive value.
//' @param mean Mean parameter of the IAR-Gamma process. A positive value.
//'
//' @return  A list with the following components:
//' \itemize{
//' \item{series}{ Array with simulated IAR-Gamma process.}
//' \item{times}{ Array with observation times.}
//' }
//' @references
//' \insertRef{Eyheramendy_2018}{iAR}
//'
//' @seealso
//'
//' \code{\link{gentime}}
//'
//' @keywords internal
//' @examples
//' \dontshow{
//' n=100
//' set.seed(6714)
//' st<-gentime(n)
//' y<-iARgsample(coef=0.9,times=st,sigma=1,mean=1)
//' plot(st,y$series,type='l')
//' hist(y$series,breaks=20)
//' }
//' @noRd
// [[Rcpp::export]]
List iARgsample(double coef, arma::vec times, int sigma = 1, int mean = 1) {
  List output;

  int n = times.size();
  arma::vec y(n, fill::zeros);
  arma::vec delta = arma::diff(times);

  y.row(0) = arma::randg(1, distr_param(1,1));

  arma::vec shape(n, fill::zeros);
  arma::vec scale(n, fill::zeros);
  arma::vec yhat(n, fill::zeros);

  for(int i = 1; i < n; ++i) {
    double phid = std::pow(coef, arma::as_scalar(delta.row(i-1)));
    yhat.row(i) = mean + phid * y.row(i-1);
    double gL = sigma * (1 - std::pow(phid, 2));
    shape.row(i) = std::pow(arma::as_scalar(yhat.row(i)), 2)/gL;
    scale.row(i) = (gL/yhat.row(i));

    double a = arma::as_scalar(shape.row(i));
    double b = arma::as_scalar(scale.row(i));

    y.row(i) = arma::randg(1, distr_param(a,b));
  }

  output["series"] = y.rows(0, n-1);
  
  return output;
}
