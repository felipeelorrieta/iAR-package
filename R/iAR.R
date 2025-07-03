#' iAR: Irregularly Observed Autoregressive Models
#'
#' Description: Data sets, functions and scripts with examples to implement autoregressive
#' models for irregularly observed time series. The models available in this package
#' are the irregular autoregressive model (Eyheramendy et al.(2018) <doi:10.1093/mnras/sty2487>),
#' the complex irregular autoregressive model (Elorrieta et al.(2019) <doi:10.1051/0004-6361/201935560>)
#' and the bivariate irregular autoregressive model (Elorrieta et al.(2021) <doi:10.1093/mnras/stab1216>)
#'
#'
#' "_PACKAGE"
#' @name iAR-package
#' @import S7 Rcpp ggplot2 zoo
#' @importFrom Rdpack reprompt
#' @importFrom stats nlminb rexp
#' @importFrom stats rnorm runif qnorm
#' @importFrom stats rt anova
#' @importFrom stats lm residuals
#' @importFrom stats optimize optim var
#' @importFrom stats pnorm rgamma
#' @importFrom stats sd na.omit
#' @useDynLib iAR, .registration = TRUE
NULL
