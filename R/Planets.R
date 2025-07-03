#' Transit of an extrasolar planet
#'
#' Time series corresponding to the residuals of the parametric model fitted by Jordan et al (2013) for a transit of an extrasolar planet.
#'
#' @format A data frame with 91 observations on the following 2 variables:
#' \describe{
#'   \item{t}{Time from mid-transit (hours).}
#'   \item{r}{Residuals of the parametric model fitted by Jordan et al (2013).}
#' }
#' @references
#' \insertRef{Jordan_etal13}{iAR}
#' @examples
#' data(Planets)
#' #plot(Planets[,1],Planets[,2],xlab='Time from mid-transit (hours)',ylab='Noise',pch=20)
"Planets"
