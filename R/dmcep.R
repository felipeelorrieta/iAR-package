#' Double Mode Cepheid.
#'
#' Time series of a double mode cepheid variable star obtained from OGLE.
#'
#' @format A data frame with 191 observations on the following 3 variables:
#' \describe{
#'   \item{t}{heliocentric Julian Day}
#'   \item{m}{magnitude}
#'   \item{merr}{measurement error of the magnitude (in mag).}
#' }
#' @details The dominant frequency computed by GLS for this light curve is 0.7410152.
#' The second frequency computed by GLS for this light curve is 0.5433353.
#' OGLE-ID:175210
#' @examples
#' data(dmcep)
#' f1=0.7410152
#' foldlc(dmcep,f1)
#' #fit=harmonicfit(dmcep,f1)
#' #f2=0.5433353
#' #foldlc(cbind(dmcep$t,fit$res,dmcep$merr),f2)
"dmcep"
