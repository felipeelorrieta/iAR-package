#include <RcppArmadillo.h>

#include <string.h>
#include <iostream>
#include <cstdio>
#include <complex>
#include <iomanip>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;
using namespace arma;

void uvec_push(arma::uvec & v, unsigned int value) {
  arma::uvec av(1);
  av.at(0) = value;
  v.insert_rows(v.n_rows, av.row(0));
}

//' Minus Log Likelihood of the BiAR Model
//'
//' This function return the negative log likelihood of the BiAR process given specific values of phiR and phiI
//'
//'
//' @param coef An array with the parameters of the BiAR model. The elements of the array are, in order, the autocorrelation and the cross-correlation coefficient of the BiAR model.
//' @param series1 Array with the observations of the first time series of the BiAR process.
//' @param series2 Array with the observations of the second time series of the BiAR process.
//' @param times Array with the irregular observational times.
//' @param series_esd1 Array with the measurements error standard deviations of the first time series of the BiAR process.
//' @param series_esd2 Array with the measurements error standard deviations of the second time series of the BiAR process.
//' @param zero_mean logical; if TRUE, the array y has zero mean; if FALSE, y has a mean different from zero.
//' @param yest An array with the estimate of a missing value in one or both time series of the bivariate process. This function recognizes a missing value with a NA. If the bivariate time series does not have a missing value, this value does not affect the computation of the likelihood.
//'
//' @return Value of the negative log likelihood evaluated in phiR and phiI.
//'
//' @references
//' \insertRef{Elorrieta_2021}{iAR}
//' @seealso
//' \code{\link{gentime}}, \code{\link{BiARsample}}
//'
//'
//' @keywords internal
//' @examples
//' \dontshow{
//' n=300
//' set.seed(6714)
//' st<-gentime(n)
//' x=BiARsample(phiR=0.9,phiI=0.3,times=st)
//' y=x$series
//' y1=y[1,]
//' y2=y[2,]
//' yerr1=rep(0,n)
//' yerr2=rep(0,n)
//' BiARphikalman(coef=c(0.8,0.2),series1=y1,series2=y2,times=st,
//' series_esd1=yerr1,series_esd2=yerr2,yest=c(0,0))
//' }
//' @noRd
// [[Rcpp::export]]
double BiARphikalman(arma::vec yest,arma::vec coef, arma::vec series1, arma::vec series2, arma::vec times, arma::vec series_esd1, arma::vec series_esd2, bool zero_mean = true) {
  arma::cx_double phi(coef[0], coef[1]);
  arma::vec y1 = series1;
  arma::vec y2 = series2;
  double phiMod = std::sqrt(std::pow(phi.real(), 2.0) + std::pow(phi.imag(), 2.0));
  if(phiMod >= 1.0) {
    return 1e10;
  }

  arma::vec y1_copy = y1;
  arma::vec y2_copy = y2;

  //Removing NA elements before join them
  uvec y1_NA_indexes = arma::find_nonfinite(y1);
  uvec y2_NA_indexes = arma::find_nonfinite(y2);

  if(y1_NA_indexes.size() > 0 && y2_NA_indexes.size() == 0) {
    y1_copy.shed_rows(y1_NA_indexes);
    y2_copy.shed_rows(y1_NA_indexes);
  }

  if(y1_NA_indexes.size() == 0 && y2_NA_indexes.size() > 0) {
    y1_copy.shed_rows(y2_NA_indexes);
    y2_copy.shed_rows(y2_NA_indexes);
  }

  if(y1_NA_indexes.size() > 0 && y2_NA_indexes.size() > 0) {
    uvec combination = join_cols(y1_NA_indexes, y2_NA_indexes);
    uvec all_joined;
    std::set<int> uniqueSet( combination.begin(), combination.end() );

    for(auto elem: uniqueSet) {
      uvec_push(all_joined, elem);
    }

    y1_copy.shed_rows(all_joined);
    y2_copy.shed_rows(all_joined);
  }

  arma::mat y0 = arma::join_horiz(y1_copy, y2_copy);
  arma::mat sigmaY = arma::cov(y0);

  if(zero_mean == false) {
    y1 = y1 - arma::mean(y1_copy);
    y2 = y2 - arma::mean(y2_copy);
  }

  int n = y1.size();
  arma::mat sigmaHat = sigmaY * arma::eye(2,2);
  arma::mat xhat(2, n, fill::zeros);
  arma::vec delta = arma::diff(times);

  arma::mat Q = sigmaHat;
  arma::mat F(2, 2, fill::zeros);
  arma::mat G = arma::eye(2,2);

  double psi = (phi.imag() < 0 && phiMod < 1) ? -acos(phi.real()/phiMod) : acos(phi.real()/phiMod);

  arma::mat y = arma::join_vert(y1.t(), y2.t());

  double sumError = 0.0;
  double sumLambda = 0.0;

  for (int i = 0; i < (n-1); ++i) {
    arma::mat R = {{std::pow(series_esd1[i+1], 2.0), 0}, {0, std::pow(series_esd2[i+1], 2.0)}};
    arma::mat Lambda(2, 2, fill::zeros);
    Lambda = G * sigmaHat * G.t() + R;

    if(det(Lambda) <= 0 || Lambda.has_nan()) {
      sumLambda = n * 1e10;
      break;
    }

    double deltaPsi = delta[i] * psi;
    double phiModDelta = pow(phiMod, delta[i]);

    double phi2Real = phiModDelta * cos(deltaPsi);
    double phi2Imag = phiModDelta * sin(deltaPsi);

    arma::mat F = {{phi2Real, -phi2Imag}, {phi2Imag, phi2Real}};

    arma::cx_double phi2Inter = pow(phi, delta[i]);
    double phi2Mod = std::pow(std::sqrt(std::pow(phi2Inter.real(), 2.0) + std::pow(phi2Inter.imag(), 2.0)), 2.0);

    auto Qt = (1 - phi2Mod) * Q;

    sumLambda = sumLambda + log(arma::det(Lambda));
    arma::mat theta = F * sigmaHat * G.t();
    arma::vec yaux = y.col(i);

    arma::mat innov;

    if(yaux.row(0).has_nan() == true && yaux.row(1).has_nan() == true) {
        arma::mat temp(2, 1, fill::zeros);
        temp.row(0) = yest[0];
        temp.row(1) = yest[1];

        innov = temp - G * xhat.col(i);
    } else if(yaux.row(0).has_nan() == true && yaux.row(1).has_nan() == false) {
      arma::mat temp(2, 1, fill::zeros);
      temp.row(0) = yest[0];
      temp.row(1) = yaux.row(1);

      innov = temp - G * xhat.col(i);
    } else if(yaux.row(0).has_nan() == false && yaux.row(1).has_nan() == true) {
      arma::mat temp(2, 1, fill::zeros);
      temp.row(0) = yaux.row(0);
      temp.row(1) = yest[0];
      innov = temp - G * xhat.col(i);
    } else {
      innov = yaux - G * xhat.col(i);
    }

    arma::mat aux = innov;
    sumError = sumError + as_scalar(aux.t() * arma::inv(Lambda) * aux);

    xhat.col(i + 1) = F * xhat.col(i) + theta * arma::inv(Lambda) * aux;

    if(arma::accu(R) == 0) {
      sigmaHat = Qt + F * (sigmaHat - sigmaHat.t()) * F.t();
    } else {
      sigmaHat = F * sigmaHat * F.t() + Qt - theta * arma::inv(Lambda) * theta.t();
    }
  }

  if(sumLambda != sumLambda) {
    return 1e10;
  } else {
    return (sumLambda + sumError)/2;
  }
}
