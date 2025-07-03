#include <RcppArmadillo.h>
#include <string.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

//' Minus Log Likelihood of the iAR Model estimated via Kalman Recursions
//'
//' This function return the negative log likelihood of the iAR process given a specific value of phi.
//'
//' @param coef A given phi coefficient of the iAR model.
//' @param series Array with the time series observations.
//' @param series_esd Array with the measurements error standard deviations.
//' @param times Array with the irregular observational times.
//' @param zero_mean logical; if TRUE, the array series has zero mean; if FALSE, series has a mean different from zero.
//' @param standardized logical; if TRUE, the array series is standardized; if FALSE, series contains the raw time series.
//' @param yest The estimate of a missing value in the time series. This function recognizes a missing value with a NA. If the time series does not have a missing value, this value does not affect the computation of the likelihood.
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
//' yerr=rep(0,100)
//' iARphikalman(coef=0.8,series=y,series_esd=yerr,times=st,yest=0)
//' }
//' @noRd
// [[Rcpp::export]]
double iARphikalman(arma::vec yest,double coef, arma::vec series, arma::vec series_esd, arma::vec times, bool zero_mean = true, bool standardized = true) {
  double x = coef;
  arma::vec y = series;
  arma::vec y_copy = series;

  //Removing NA elements before join them
  uvec y_NA_indexes = arma::find_nonfinite(y);

  if(y_NA_indexes.size() > 0) {
    y_copy.shed_rows(y_NA_indexes);
  }

  double sigmay = 1.0;

  if(zero_mean == false) {
    y = y - arma::mean(y_copy);
  }


  if(standardized == false) {
    sigmay = arma::var(y_copy);
  }

  int n = y.size();

  arma::vec delta(n, fill::zeros);
  arma::mat xhat(1, n, fill::zeros);
  arma::vec Sighat(1, 1, fill::zeros);

  Sighat[0] = sigmay;
  delta = arma::diff(times);

  arma::vec Q(1, 1, fill::zeros);
  Q[0] = Sighat[0];

  arma::vec F(1,1, fill::zeros);
  arma::vec G(1,1, fill::ones);

  double sumLambda = 0;
  double sumError = 0;

  //Check if x is NaN
  if(x != x) {
    x = 1.1;
  }

  if(std::abs(x) < 1.0) {
    for (int i = 0; i < (n-1); ++i) {
      arma::vec Lambda(1,1, fill::ones);
      arma::vec Theta(1,1, fill::ones);

      Lambda = G * Sighat * G.t() + arma::pow(series_esd.row(i+1), 2);

      if(Lambda[0] <= 0 || Lambda.has_nan()) {
        sumLambda = n * 1e10;
        break;
      }

      double phi2 = std::pow(x, arma::as_scalar(delta.row(i)));
      F[0] = phi2;

      phi2 = 1 - (std::pow(x, (2 * arma::as_scalar(delta.row(i)))));
      auto Qt = phi2 * Q;

      sumLambda = sumLambda + arma::as_scalar(arma::log(Lambda));
      Theta = F * Sighat * G.t();

      arma::vec yaux = y.row(i);
      arma::mat innov;

      innov = yaux - G * xhat.col(i);

      if(yaux.has_nan() == true) {
        arma::mat temp(1, 1, fill::zeros);
        temp.row(0) = yest;
        innov = temp - G * xhat.col(i);
      }


      sumError = sumError + arma::as_scalar((arma::pow(innov, 2)/Lambda));
      xhat.col(i + 1) = F * xhat.col(i) + Theta * arma::inv(Lambda) * innov;
      Sighat = F * Sighat * F.t() + Qt - Theta * arma::inv(Lambda) * Theta.t();
    }

    if(sumLambda != sumLambda) {
      return 1e10;
    } else {
      return ((sumLambda + sumError)/n);
    }

  } else {
    return 1e10;
  }
}
