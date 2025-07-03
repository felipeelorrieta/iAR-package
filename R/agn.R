#' @name agn
#' @title Active Galactic Nuclei
#' @description Time series of the AGN MCG-6-30-15 measured in the K-band between 2006 August and 2011 July with the ANDICAM camera mounted on the 1.3 m telescope at Cerro Tololo Inter-American Observatory (CTIO)
#' @docType data
#' @usage agn
#' @format A data frame with 237 observations on the following 3 variables:
#' \describe{
#'   \item{t}{heliocentric Julian Day - 2450000}
#'   \item{m}{Flux $(10^(-15) ergs/s/cm^2 /A)$}
#'   \item{merr}{measurement error standard deviations.}
#' }
#' @references
#' \insertRef{Lira_2015}{iAR}
#' @examples
#' data(agn)
#' plot(agn$t,agn$m,type="l",ylab="",xlab="")
"agn"
