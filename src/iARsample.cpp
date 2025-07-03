#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;
using namespace arma;

//' Simulate from an iAR Model
//'
//' Simulates an iAR Time Series Model.
//'
//' @param coef A coefficient of iAR model. A value between 0 and 1
//' @param times Array with observational times.
//'
//' @return A list with the following components:
//' \itemize{
//' \item{series}{ Array with simulated iAR data.}
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
//' set.seed(6714)
//' st<-gentime(n=100)
//' y<-iARsample(coef=0.99,times=st)
//' y<-y$series
//' plot(st,y,type='l')
//' }
//' @noRd
// [[Rcpp::export]]
List iARsample(double coef, arma::vec times) {
  List output;
  
  int n = times.size();
  arma::mat Sigma(n, n, fill::zeros);

  for(int i = 0; i < n; ++i) {
    arma::vec d(n-i, fill::zeros);
    d = times[i] - times.rows(i, n-1);
    d = arma::abs(d);

    arma::vec temp(n-i, fill::zeros);
    for(int j = 0; j < (n-i); ++j) {
      temp[j] = std::pow(coef, d[j]);
    }

    Sigma.rows(i, (n-1)).col(i) = temp;
    Sigma.cols(i, (n-1)).row(i) = temp.t();

  }

  Function eigenRfunction("eigen");
  NumericMatrix sigmaRcpp = wrap(Sigma);
  Rcpp::List eigen_results = eigenRfunction(Named("x") = sigmaRcpp, Named("symmetric") = true);
  arma::vec eigenValues = Rcpp::as<arma::vec>(eigen_results[0]);
  arma::mat eigenVectors = Rcpp::as<arma::mat>(eigen_results[1]);

  arma::mat temp = arma::eye(n, n);
  temp.diag() = arma::sqrt(eigenValues);

  arma::mat A = eigenVectors * temp * eigenVectors.t();
  arma::vec e = arma::randn(n);
  arma::vec y = arma::vectorise(A * e);

  output["series"] = y.rows(0, n-1);

  return output;
}
