#' Simulate from a BiAR Model
#'
#' Simulates a BiAR Time Series Model
#'
#' @param times Array with observational times.
#' @param phiR Autocorrelation coefficient of BiAR model. A value between -1 and 1.
#' @param phiI Crosscorrelation coefficient of BiAR model. A value between -1 and 1.
#' @param series_esd1 Array with the measurements error standard deviations of the first time series of the bivariate process.
#' @param series_esd2 Array with the measurements error standard deviations of the second time series of the bivariate process.
#' @param rho Contemporary correlation coefficient of BiAR model. A value between -1 and 1.
#'
#' @details The chosen phiR and phiI values must satisfy the condition $|phiR + i phiI| < 1$.
#'
#' @return A list with the following components:
#' \itemize{
#' \item{series}{ Matrix with the simulated BiAR process.}
#' }
#'
#'
#' @references
#' \insertRef{Elorrieta_2021}{iAR}
#' @seealso
#' \code{\link{gentime}}
#'
#' @keywords internal
#' @examples
#' \dontshow{
#' n=300
#' set.seed(6714)
#' st<-gentime(n)
#' x=BiARsample(phiR=0.9,phiI=0.3,times=st)
#' plot(st,x$series[1,],type='l')
#' plot(st,x$series[2,],type='l')
#' x=BiARsample(phiR=-0.9,phiI=-0.3,times=st)
#' plot(st,x$series[1,],type='l')
#' plot(st,x$series[2,],type='l')
#' }
#' @noRd
BiARsample<-function (times, phiR, phiI, series_esd1=0,series_esd2=0,rho = 0)
{
  n <- length(times)
  delta <- diff(times)
  x = matrix(0, nrow = 2, ncol = n)
  F = matrix(0, nrow = 2, ncol = 2)
  phi = complex(1, real = phiR, imaginary = phiI)
  if (Mod(phi) >= 1)
    stop("Mod of Phi must be less than one")
  Phi = Mod(phi)
  psi <- acos(phiR/Phi)
  if (phiI < 0)
    psi = -acos(phiR/Phi)
  #psi2 <- asin(phiI/Phi)
  e.R = rnorm(n)
  e.I = rnorm(n)
  state.error = rbind(e.R, e.I)
  Sigma = matrix(1, nrow = 2, ncol = 2)
  Sigma[1, 1] = 1
  Sigma[2, 2] = 1
  Sigma[1, 2] = rho * sqrt(Sigma[1, 1]) * sqrt(Sigma[2, 2])
  Sigma[2, 1] = Sigma[1, 2]
  B = svd(Sigma)
  A = matrix(0, nrow = 2, ncol = 2)
  diag(A) = sqrt(B$d)
  Sigma.root = (B$u) %*% A %*% B$u
  state.error = Sigma.root %*% state.error
  #measurement error
  w.R = rnorm(n)
  w.I = rnorm(n)
  observation.error = rbind(w.R, w.I)
  SigmaO = matrix(0, nrow = 2, ncol = 2)
  SigmaO[1, 1] = series_esd1**2
  SigmaO[2, 2] = series_esd2**2
  BO = svd(SigmaO)
  AO = matrix(0, nrow = 2, ncol = 2)
  diag(AO) = sqrt(BO$d)
  SigmaO.root = (BO$u) %*% AO %*% BO$u
  observation.error = SigmaO.root %*% observation.error
  G = diag(2)
  y = observation.error
  x[, 1] = state.error[, 1]
  for (i in 1:(n - 1)) {
    phi2.R <- (Phi^delta[i]) * cos(delta[i] * psi)
    phi2.I <- (Phi^delta[i]) * sin(delta[i] * psi)
    phi2 <- 1 - Mod(phi^delta[i])^2
    F[1, 1] = phi2.R
    F[1, 2] = -phi2.I
    F[2, 1] = phi2.I
    F[2, 2] = phi2.R
    x[, i + 1] = F %*% x[, i] + sqrt(phi2) * state.error[,i]
    y[, i ] = G %*% x[, i] + observation.error[,i]
  }
  y[, n ] = G %*% x[, n] + observation.error[,n]
  return(list(series = y))
}
