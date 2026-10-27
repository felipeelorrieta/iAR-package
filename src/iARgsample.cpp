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
List iARgsample(double coef, arma::vec times, double sigma = 1.0, double mean = 1.0) {
  List output;

  int n = times.size();
  arma::vec y(n, fill::zeros);
  arma::vec delta = arma::diff(times);

  y.row(0) = arma::randg(1, distr_param(1, 1));

  arma::vec shape(n, fill::zeros);
  arma::vec scale(n, fill::zeros);
  arma::vec yhat(n, fill::zeros);

  for(int i = 1; i < n; ++i) {
    double phid   = std::pow(coef, arma::as_scalar(delta.row(i-1)));
    double yhat_i = mean + phid * arma::as_scalar(y.row(i-1));
    double gL     = sigma * (1.0 - std::pow(phid, 2.0));

    // Guard: gL can collapse to 0 when phi ≈ 1 and delta ≈ 0.
    // Clamp to a small positive floor to avoid shape/scale blowing up.
    if (gL <= 0.0) {
      Rcpp::warning(
        "iAR-Gamma: gL <= 0 at step %d (phi^(2*delta) >= 1). "
        "This can occur when consecutive times are very close together "
        "combined with phi near 1. Results at this step may be unreliable.",
        i
      );
      gL = 1e-10;
    }

    // Guard: yhat must be positive for a valid Gamma parameterisation.
    if (yhat_i <= 0.0) {
      Rcpp::warning(
        "iAR-Gamma: conditional mean (yhat) <= 0 at step %d. "
        "Ensure that the 'mean' parameter is strictly positive. "
        "Results at this step may be unreliable.",
        i
      );
      yhat_i = 1e-10;
    }

    shape.row(i) = std::pow(yhat_i, 2.0) / gL;
    scale.row(i) = gL / yhat_i;

    double a = arma::as_scalar(shape.row(i));
    double b = arma::as_scalar(scale.row(i));

    y.row(i) = arma::randg(1, distr_param(a, b));
  }

  output["series"] = y.rows(0, n-1);
  return output;
}
